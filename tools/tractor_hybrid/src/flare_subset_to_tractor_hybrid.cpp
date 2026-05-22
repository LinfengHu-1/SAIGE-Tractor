// flare_subset_to_tractor_hybrid.cpp
// designed by Kai, implemented by codex
//
// Convert phased genotype VCF/BCF plus FLARE local ancestry VCF/BCF into
// ancestry-aware packed files for a SAIGE-TRACTOR genotype backend.

#include <htslib/hts.h>
#include <htslib/vcf.h>

#include <algorithm>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <unordered_set>
#include <vector>

struct RareCarrierPacked {
    uint32_t pos_index; // split biallelic variant ordinal
    uint32_t anc_hap;   // high 5 bits ancestry, low 27 bits hap_id
};

struct CarrierHap {
    uint32_t hap_id;
    uint8_t ancestry;
};

struct AncestryState {
    std::vector<std::vector<uint64_t>> masks;
    std::vector<int8_t> hap_ancestry;
};

struct OpenAncestryBlock {
    bool active = false;
    uint32_t block_id = 0;
    int geno_rid = -1;
    std::string chr;
    int64_t start_pos = -1;
    int64_t end_pos = -1;
    AncestryState state;
};

struct LaiRecord {
    bool valid = false;
    int geno_rid = -1;
    std::string chr;
    int64_t pos = 0;
    AncestryState state;
};

static constexpr uint32_t kMaxPackedHapId = (1u << 27) - 1u;

[[noreturn]] static void die(const char* fmt, ...) {
    va_list args;
    va_start(args, fmt);
    std::fputs("ERROR: ", stderr);
    std::vfprintf(stderr, fmt, args);
    std::fputc('\n', stderr);
    va_end(args);
    std::exit(1);
}

static inline uint32_t pack_anc_hap(uint32_t ancestry, uint32_t hap_id) {
    return ((ancestry & 31u) << 27) | (hap_id & 0x07ffffffu);
}

static inline void set_bit(std::vector<uint64_t>& bits, uint32_t h) {
    bits[h >> 6] |= (1ULL << (h & 63));
}

static inline bool is_gt_vector_end(int32_t value) {
    return value == bcf_int32_vector_end;
}

static bool same_masks(
    const std::vector<std::vector<uint64_t>>& a,
    const std::vector<std::vector<uint64_t>>& b
) {
    if (a.size() != b.size()) return false;

    for (size_t k = 0; k < a.size(); ++k) {
        if (a[k] != b[k]) return false;
    }

    return true;
}

static uint64_t tell_or_die(FILE* fp, const char* path) {
    off_t offset = ftello(fp);
    if (offset < 0) die("ftello failed for %s", path);
    return static_cast<uint64_t>(offset);
}

static FILE* open_output_or_die(const std::string& path, const char* mode) {
    FILE* fp = std::fopen(path.c_str(), mode);
    if (!fp) die("cannot open output file: %s", path.c_str());

    if (setvbuf(fp, nullptr, _IOFBF, 1 << 24) != 0) {
        die("setvbuf failed for output file: %s", path.c_str());
    }

    return fp;
}

static void write_exact(FILE* fp, const void* data, size_t bytes, const char* what) {
    if (bytes == 0) return;
    if (std::fwrite(data, 1, bytes, fp) != bytes) {
        die("failed writing %s", what);
    }
}

static void write_u32_le(FILE* fp, uint32_t value, const char* what) {
    uint8_t bytes[4] = {
        static_cast<uint8_t>(value),
        static_cast<uint8_t>(value >> 8),
        static_cast<uint8_t>(value >> 16),
        static_cast<uint8_t>(value >> 24)
    };
    write_exact(fp, bytes, sizeof(bytes), what);
}

static void write_u64_le(FILE* fp, uint64_t value, const char* what) {
    uint8_t bytes[8] = {
        static_cast<uint8_t>(value),
        static_cast<uint8_t>(value >> 8),
        static_cast<uint8_t>(value >> 16),
        static_cast<uint8_t>(value >> 24),
        static_cast<uint8_t>(value >> 32),
        static_cast<uint8_t>(value >> 40),
        static_cast<uint8_t>(value >> 48),
        static_cast<uint8_t>(value >> 56)
    };
    write_exact(fp, bytes, sizeof(bytes), what);
}

static void write_i64_le(FILE* fp, int64_t value, const char* what) {
    write_u64_le(fp, static_cast<uint64_t>(value), what);
}

static void write_string(FILE* fp, const std::string& value, const char* what) {
    if (value.size() > std::numeric_limits<uint32_t>::max()) {
        die("string too long while writing %s", what);
    }

    write_u32_le(fp, static_cast<uint32_t>(value.size()), what);
    write_exact(fp, value.data(), value.size(), what);
}

static void write_magic_header(FILE* fp, const char magic[8]) {
    write_exact(fp, magic, 8, "index header");
}

static void write_common_offset_idx_record(
    FILE* fp,
    uint64_t common_index,
    uint32_t global_variant_index,
    uint64_t mks_offset,
    uint64_t geno_offset
) {
    write_u64_le(fp, common_index, "common offset index");
    write_u32_le(fp, global_variant_index, "common offset index");
    write_u64_le(fp, mks_offset, "common offset index");
    write_u64_le(fp, geno_offset, "common offset index");
}

static void write_rare_offset_idx_record(
    FILE* fp,
    uint64_t rare_index,
    uint32_t global_variant_index,
    uint64_t mks_offset,
    uint64_t carrier_offset,
    uint32_t n_carriers
) {
    write_u64_le(fp, rare_index, "rare offset index");
    write_u32_le(fp, global_variant_index, "rare offset index");
    write_u64_le(fp, mks_offset, "rare offset index");
    write_u64_le(fp, carrier_offset, "rare offset index");
    write_u32_le(fp, n_carriers, "rare offset index");
}

static void write_anc_offset_idx_record(
    FILE* fp,
    uint32_t block_id,
    uint64_t mks_offset,
    uint64_t anc_offset
) {
    write_u32_le(fp, block_id, "ancestry offset index");
    write_u64_le(fp, mks_offset, "ancestry offset index");
    write_u64_le(fp, anc_offset, "ancestry offset index");
}

static void write_common_mks_record(
    FILE* fp,
    uint64_t common_index,
    uint32_t global_variant_index,
    const char* chr,
    int64_t pos,
    const std::string& id,
    const char* ref,
    const char* alt,
    uint32_t alt_index,
    uint32_t block_id,
    uint64_t geno_offset,
    uint32_t mac
) {
    write_u64_le(fp, common_index, "common mks");
    write_u32_le(fp, global_variant_index, "common mks");
    write_string(fp, chr, "common mks chr");
    write_i64_le(fp, pos, "common mks pos");
    write_string(fp, id, "common mks id");
    write_string(fp, ref, "common mks ref");
    write_string(fp, alt, "common mks alt");
    write_u32_le(fp, alt_index, "common mks");
    write_u32_le(fp, block_id, "common mks");
    write_u64_le(fp, geno_offset, "common mks");
    write_u32_le(fp, mac, "common mks");
}

static void write_rare_mks_record(
    FILE* fp,
    uint64_t rare_index,
    uint32_t global_variant_index,
    const char* chr,
    int64_t pos,
    const std::string& id,
    const char* ref,
    const char* alt,
    uint32_t alt_index,
    uint64_t carrier_offset,
    uint32_t n_carriers,
    uint32_t mac
) {
    write_u64_le(fp, rare_index, "rare mks");
    write_u32_le(fp, global_variant_index, "rare mks");
    write_string(fp, chr, "rare mks chr");
    write_i64_le(fp, pos, "rare mks pos");
    write_string(fp, id, "rare mks id");
    write_string(fp, ref, "rare mks ref");
    write_string(fp, alt, "rare mks alt");
    write_u32_le(fp, alt_index, "rare mks");
    write_u64_le(fp, carrier_offset, "rare mks");
    write_u32_le(fp, n_carriers, "rare mks");
    write_u32_le(fp, mac, "rare mks");
}

static void write_anc_mks_record(
    FILE* fp,
    uint32_t block_id,
    const char* chr,
    int64_t start_pos,
    int64_t end_pos,
    uint64_t anc_offset
) {
    write_u32_le(fp, block_id, "ancestry mks");
    write_string(fp, chr, "ancestry mks chr");
    write_i64_le(fp, start_pos, "ancestry mks start");
    write_i64_le(fp, end_pos, "ancestry mks end");
    write_u64_le(fp, anc_offset, "ancestry mks offset");
}

static void write_sidecars(
    bcf_hdr_t* ghdr,
    const std::string& samples_path,
    const std::string& meta_path,
    const char* geno_vcf,
    const char* flare_vcf,
    int n_samples,
    uint64_t n_haps,
    int n_words,
    int n_ancestries,
    int rare_threshold
) {
    FILE* samples_fp = open_output_or_die(samples_path, "w");
    for (int i = 0; i < n_samples; ++i) {
        std::fprintf(samples_fp, "%s\n", ghdr->samples[i]);
    }
    std::fclose(samples_fp);

    FILE* meta_fp = open_output_or_die(meta_path, "w");
    std::fprintf(meta_fp, "format_version\t1\n");
    std::fprintf(meta_fp, "n_samples\t%d\n", n_samples);
    std::fprintf(meta_fp, "n_haps\t%llu\n", static_cast<unsigned long long>(n_haps));
    std::fprintf(meta_fp, "n_words\t%d\n", n_words);
    std::fprintf(meta_fp, "n_ancestries\t%d\n", n_ancestries);
    std::fprintf(meta_fp, "rare_threshold\t%d\n", rare_threshold);
    std::fprintf(meta_fp, "source_genotype\t%s\n", geno_vcf);
    std::fprintf(meta_fp, "source_flare\t%s\n", flare_vcf);
    std::fclose(meta_fp);
}

static void write_anc_block(
    FILE* anc_fp,
    const char* anc_bin_path,
    FILE* anc_mks_fp,
    const char* anc_mks_path,
    FILE* anc_idx_fp,
    uint32_t block_id,
    const char* chr,
    int64_t start_pos,
    int64_t end_pos,
    const AncestryState& state,
    int n_ancestries,
    int n_words
) {
    uint64_t offset = tell_or_die(anc_fp, anc_bin_path);

    for (int k = 0; k < n_ancestries; ++k) {
        size_t wrote = std::fwrite(
            state.masks[k].data(),
            sizeof(uint64_t),
            static_cast<size_t>(n_words),
            anc_fp
        );

        if (wrote != static_cast<size_t>(n_words)) {
            die("failed writing ancestry block %u", block_id);
        }
    }

    uint64_t mks_offset = tell_or_die(anc_mks_fp, anc_mks_path);
    write_anc_mks_record(anc_mks_fp, block_id, chr, start_pos, end_pos, offset);
    write_anc_offset_idx_record(anc_idx_fp, block_id, mks_offset, offset);
}

static int check_sample_order(bcf_hdr_t* h1, bcf_hdr_t* h2) {
    int n1 = bcf_hdr_nsamples(h1);
    int n2 = bcf_hdr_nsamples(h2);
    if (n1 != n2) return 0;

    for (int i = 0; i < n1; ++i) {
        if (std::strcmp(h1->samples[i], h2->samples[i]) != 0) {
            std::fprintf(
                stderr,
                "Sample mismatch at %d: %s vs %s\n",
                i,
                h1->samples[i],
                h2->samples[i]
            );
            return 0;
        }
    }

    return 1;
}

static void build_state_from_flare(
    bcf_hdr_t* ahdr,
    bcf1_t* arec,
    int n_samples,
    int n_ancestries,
    int n_words,
    int32_t** an1_arr,
    int* nan1_arr,
    int32_t** an2_arr,
    int* nan2_arr,
    AncestryState& state
) {
    int n_an1 = bcf_get_format_int32(ahdr, arec, "AN1", an1_arr, nan1_arr);
    int n_an2 = bcf_get_format_int32(ahdr, arec, "AN2", an2_arr, nan2_arr);

    if (n_an1 <= 0 || n_an2 <= 0) {
        die("FLARE VCF must contain FORMAT/AN1 and FORMAT/AN2");
    }

    if (n_an1 != n_samples || n_an2 != n_samples) {
        die("FLARE VCF must contain scalar FORMAT/AN1 and FORMAT/AN2");
    }

    state.masks.assign(
        static_cast<size_t>(n_ancestries),
        std::vector<uint64_t>(static_cast<size_t>(n_words), 0)
    );
    state.hap_ancestry.assign(static_cast<size_t>(2 * n_samples), -1);

    for (int i = 0; i < n_samples; ++i) {
        int32_t raw_a1 = (*an1_arr)[i];
        int32_t raw_a2 = (*an2_arr)[i];

        if (raw_a1 == bcf_int32_missing || raw_a1 == bcf_int32_vector_end) {
            die("missing FORMAT/AN1 in FLARE record at sample index %d", i);
        }

        if (raw_a2 == bcf_int32_missing || raw_a2 == bcf_int32_vector_end) {
            die("missing FORMAT/AN2 in FLARE record at sample index %d", i);
        }

        int a1 = static_cast<int>(raw_a1);
        int a2 = static_cast<int>(raw_a2);

        if (a1 < 0 || a1 >= n_ancestries) {
            die("FORMAT/AN1 value out of range at sample index %d: %d", i, a1);
        }

        if (a2 < 0 || a2 >= n_ancestries) {
            die("FORMAT/AN2 value out of range at sample index %d: %d", i, a2);
        }

        uint32_t hap0 = static_cast<uint32_t>(2 * i);
        uint32_t hap1 = static_cast<uint32_t>(2 * i + 1);
        set_bit(state.masks[a1], hap0);
        set_bit(state.masks[a2], hap1);
        state.hap_ancestry[hap0] = static_cast<int8_t>(a1);
        state.hap_ancestry[hap1] = static_cast<int8_t>(a2);
    }
}

static void close_open_block(
    FILE* anc_fp,
    const char* anc_bin_path,
    FILE* anc_mks_fp,
    const char* anc_mks_path,
    FILE* anc_idx_fp,
    OpenAncestryBlock& block,
    uint32_t& n_blocks_written,
    int n_ancestries,
    int n_words
) {
    if (!block.active) return;

    write_anc_block(
        anc_fp,
        anc_bin_path,
        anc_mks_fp,
        anc_mks_path,
        anc_idx_fp,
        block.block_id,
        block.chr.c_str(),
        block.start_pos,
        block.end_pos,
        block.state,
        n_ancestries,
        n_words
    );

    ++n_blocks_written;
}

static void update_open_block_from_flare(
    FILE* anc_fp,
    const char* anc_bin_path,
    FILE* anc_mks_fp,
    const char* anc_mks_path,
    FILE* anc_idx_fp,
    OpenAncestryBlock& block,
    uint32_t& n_blocks_written,
    int n_ancestries,
    int n_words,
    int flare_geno_rid,
    const char* flare_chr,
    int64_t interval_start,
    int64_t interval_end,
    AncestryState& new_state
) {
    if (interval_start > interval_end) {
        die(
            "invalid LAI interval for %s: %lld > %lld",
            flare_chr,
            static_cast<long long>(interval_start),
            static_cast<long long>(interval_end)
        );
    }

    if (!block.active) {
        block.active = true;
        block.block_id = n_blocks_written;
        block.geno_rid = flare_geno_rid;
        block.chr = flare_chr;
        block.start_pos = interval_start;
        block.end_pos = interval_end;
        block.state = std::move(new_state);
        return;
    }

    bool same_chromosome = block.geno_rid == flare_geno_rid;
    bool unchanged = same_chromosome && same_masks(block.state.masks, new_state.masks);

    if (unchanged) {
        block.end_pos = std::max(block.end_pos, interval_end);
        return;
    }

    close_open_block(
        anc_fp,
        anc_bin_path,
        anc_mks_fp,
        anc_mks_path,
        anc_idx_fp,
        block,
        n_blocks_written,
        n_ancestries,
        n_words
    );

    block.active = true;
    block.block_id = n_blocks_written;
    block.geno_rid = flare_geno_rid;
    block.chr = flare_chr;
    block.start_pos = interval_start;
    block.end_pos = interval_end;
    block.state = std::move(new_state);
}

static std::string make_split_id(
    const char* raw_id,
    const char* chr,
    int64_t pos,
    const char* ref,
    const char* alt
) {
    if (raw_id == nullptr || std::strcmp(raw_id, ".") == 0) {
        return std::string(chr) + ":" + std::to_string(pos) + ":" + ref + ":" + alt;
    }

    return std::string(raw_id) + "_" + ref + "_" + alt;
}

static void add_alt_carrier(
    std::vector<std::vector<CarrierHap>>& carriers_by_alt,
    const AncestryState& state,
    int allele,
    uint32_t hap_id
) {
    if (allele <= 0 || allele >= static_cast<int>(carriers_by_alt.size())) return;

    int8_t ancestry = state.hap_ancestry[hap_id];
    if (ancestry < 0) return;

    carriers_by_alt[allele].push_back(
        CarrierHap{hap_id, static_cast<uint8_t>(ancestry)}
    );
}

static void validate_extra_ploidy(
    const int32_t* sample_gt,
    int ploidy,
    const char* chr,
    int64_t pos,
    int sample_index
) {
    for (int p = 2; p < ploidy; ++p) {
        int32_t gt = sample_gt[p];
        if (is_gt_vector_end(gt)) break;
        die(
            "non-diploid genotype at %s:%lld sample index %d",
            chr,
            static_cast<long long>(pos),
            sample_index
        );
    }
}

static void process_genotypes(
    bcf_hdr_t* ghdr,
    bcf1_t* grec,
    const OpenAncestryBlock& block,
    int n_samples,
    int32_t** gt_arr,
    int* ngt_arr,
    std::vector<std::vector<CarrierHap>>& carriers_by_alt
) {
    const char* chr = bcf_hdr_id2name(ghdr, grec->rid);
    int64_t pos = static_cast<int64_t>(grec->pos) + 1;

    int ngt = bcf_get_genotypes(ghdr, grec, gt_arr, ngt_arr);
    if (ngt <= 0) {
        die("missing FORMAT/GT at %s:%lld", chr, static_cast<long long>(pos));
    }

    if (ngt % n_samples != 0) {
        die("GT field length is not divisible by sample count at %s:%lld", chr, static_cast<long long>(pos));
    }

    int ploidy = ngt / n_samples;
    if (ploidy < 2) {
        die("expected diploid GT at %s:%lld", chr, static_cast<long long>(pos));
    }

    carriers_by_alt.assign(static_cast<size_t>(grec->n_allele), std::vector<CarrierHap>{});

    for (int i = 0; i < n_samples; ++i) {
        const int32_t* sample_gt = *gt_arr + static_cast<size_t>(i) * ploidy;
        int32_t g0 = sample_gt[0];
        int32_t g1 = sample_gt[1];

        validate_extra_ploidy(sample_gt, ploidy, chr, pos, i);

        if (is_gt_vector_end(g0) || is_gt_vector_end(g1) ||
            bcf_gt_is_missing(g0) || bcf_gt_is_missing(g1)) {
            die(
                "missing genotype at %s:%lld sample index %d",
                chr,
                static_cast<long long>(pos),
                i
            );
        }

        int allele0 = bcf_gt_allele(g0);
        int allele1 = bcf_gt_allele(g1);

        if (!bcf_gt_is_phased(g1)) {
            die(
                "unphased genotype at %s:%lld sample index %d",
                chr,
                static_cast<long long>(pos),
                i
            );
        }

        add_alt_carrier(
            carriers_by_alt,
            block.state,
            allele0,
            static_cast<uint32_t>(2 * i)
        );
        add_alt_carrier(
            carriers_by_alt,
            block.state,
            allele1,
            static_cast<uint32_t>(2 * i + 1)
        );
    }
}

static int read_next_record(htsFile* fp, bcf_hdr_t* hdr, bcf1_t* rec) {
    int ret = bcf_read(fp, hdr, rec);
    if (ret == 0) bcf_unpack(rec, BCF_UN_STR);
    return ret;
}

static bool read_next_lai_record(
    htsFile* fp,
    bcf_hdr_t* ahdr,
    bcf_hdr_t* ghdr,
    bcf1_t* rec,
    int n_samples,
    int n_ancestries,
    int n_words,
    int32_t** an1_arr,
    int* nan1_arr,
    int32_t** an2_arr,
    int* nan2_arr,
    std::unordered_set<std::string>& warned_missing_flare_contigs,
    LaiRecord& out
) {
    while (read_next_record(fp, ahdr, rec) == 0) {
        const char* a_chr = bcf_hdr_id2name(ahdr, rec->rid);
        int a_geno_rid = bcf_hdr_name2id(ghdr, a_chr);

        if (a_geno_rid < 0) {
            std::string contig(a_chr);
            if (warned_missing_flare_contigs.insert(contig).second) {
                std::fprintf(
                    stderr,
                    "WARNING: skipping FLARE contig absent from genotype header: %s\n",
                    a_chr
                );
            }
            continue;
        }

        out.valid = true;
        out.geno_rid = a_geno_rid;
        out.chr = a_chr;
        out.pos = static_cast<int64_t>(rec->pos) + 1;

        build_state_from_flare(
            ahdr,
            rec,
            n_samples,
            n_ancestries,
            n_words,
            an1_arr,
            nan1_arr,
            an2_arr,
            nan2_arr,
            out.state
        );

        return true;
    }

    out = LaiRecord{};
    return false;
}

static void print_usage(const char* prog) {
    std::fprintf(
        stderr,
        "Usage:\n"
        "  %s genotype.phased.vcf.gz flare.anc.vcf.gz n_ancestries rare_threshold out_prefix\n\n"
        "Example:\n"
        "  %s chr1.phased.vcf.gz chr1.flare.anc.vcf.gz 3 512 chr1\n",
        prog,
        prog
    );
}

int main(int argc, char** argv) {
    if (argc < 6) {
        print_usage(argv[0]);
        return 1;
    }

    const char* geno_vcf = argv[1];
    const char* flare_vcf = argv[2];
    int n_ancestries = std::atoi(argv[3]);
    int rare_threshold = std::atoi(argv[4]);
    const char* out_prefix = argv[5];

    if (n_ancestries <= 0 || n_ancestries > 32) {
        die("n_ancestries must be in [1, 32]");
    }

    if (rare_threshold < 0) {
        die("rare_threshold must be non-negative");
    }

    htsFile* gfp = bcf_open(geno_vcf, "r");
    htsFile* afp = bcf_open(flare_vcf, "r");

    if (!gfp || !afp) {
        die("cannot open input VCF/BCF");
    }

    bcf_hdr_t* ghdr = bcf_hdr_read(gfp);
    bcf_hdr_t* ahdr = bcf_hdr_read(afp);

    if (!ghdr || !ahdr) {
        die("cannot read input headers");
    }

    if (!check_sample_order(ghdr, ahdr)) {
        die("genotype and FLARE sample IDs must be identical and in the same order");
    }

    int n_samples = bcf_hdr_nsamples(ghdr);
    if (n_samples <= 0) {
        die("genotype VCF has no samples");
    }

    uint64_t n_haps = static_cast<uint64_t>(n_samples) * 2ULL;
    if (n_haps > static_cast<uint64_t>(kMaxPackedHapId) + 1ULL) {
        die("hap_id exceeds 27-bit packed limit");
    }

    int n_words = static_cast<int>((n_haps + 63ULL) / 64ULL);

    std::string common_bin = std::string(out_prefix) + ".common.geno.bin";
    std::string common_mks = std::string(out_prefix) + ".common.variant.mks";
    std::string common_idx = std::string(out_prefix) + ".common.variant.idx";
    std::string rare_bin = std::string(out_prefix) + ".rare.carrier.bin";
    std::string rare_mks = std::string(out_prefix) + ".rare.variant.mks";
    std::string rare_idx = std::string(out_prefix) + ".rare.variant.idx";
    std::string anc_bin = std::string(out_prefix) + ".ancblock.bin";
    std::string anc_mks = std::string(out_prefix) + ".ancblock.mks";
    std::string anc_idx = std::string(out_prefix) + ".ancblock.idx";
    std::string samples_path = std::string(out_prefix) + ".samples";
    std::string meta_path = std::string(out_prefix) + ".meta";

    FILE* common_fp = open_output_or_die(common_bin, "wb");
    FILE* common_mks_fp = open_output_or_die(common_mks, "wb");
    FILE* common_idx_fp = open_output_or_die(common_idx, "wb");
    FILE* rare_fp = open_output_or_die(rare_bin, "wb");
    FILE* rare_mks_fp = open_output_or_die(rare_mks, "wb");
    FILE* rare_idx_fp = open_output_or_die(rare_idx, "wb");
    FILE* anc_fp = open_output_or_die(anc_bin, "wb");
    FILE* anc_mks_fp = open_output_or_die(anc_mks, "wb");
    FILE* anc_idx_fp = open_output_or_die(anc_idx, "wb");

    const char common_mks_magic[8] = {'T', 'R', 'C', 'M', 'M', 'K', 'S', '1'};
    const char rare_mks_magic[8] = {'T', 'R', 'R', 'A', 'M', 'K', 'S', '1'};
    const char anc_mks_magic[8] = {'T', 'R', 'A', 'N', 'M', 'K', 'S', '1'};
    const char common_idx_magic[8] = {'T', 'R', 'C', 'M', 'I', 'D', 'X', '2'};
    const char rare_idx_magic[8] = {'T', 'R', 'R', 'A', 'I', 'D', 'X', '2'};
    const char anc_idx_magic[8] = {'T', 'R', 'A', 'N', 'I', 'D', 'X', '2'};
    write_magic_header(common_mks_fp, common_mks_magic);
    write_magic_header(rare_mks_fp, rare_mks_magic);
    write_magic_header(anc_mks_fp, anc_mks_magic);
    write_magic_header(common_idx_fp, common_idx_magic);
    write_magic_header(rare_idx_fp, rare_idx_magic);
    write_magic_header(anc_idx_fp, anc_idx_magic);

    write_sidecars(
        ghdr,
        samples_path,
        meta_path,
        geno_vcf,
        flare_vcf,
        n_samples,
        n_haps,
        n_words,
        n_ancestries,
        rare_threshold
    );

    bcf1_t* grec = bcf_init();
    bcf1_t* arec = bcf_init();

    int32_t* gt_arr = nullptr;
    int ngt_arr = 0;
    int32_t* an1_arr = nullptr;
    int nan1_arr = 0;
    int32_t* an2_arr = nullptr;
    int nan2_arr = 0;

    OpenAncestryBlock block;
    std::vector<std::vector<CarrierHap>> carriers_by_alt;
    std::unordered_set<std::string> warned_missing_flare_contigs;

    uint32_t n_blocks_written = 0;
    uint64_t common_index = 0;
    uint64_t rare_index = 0;
    uint32_t global_variant_index = 0;

    LaiRecord lai_record;
    bool has_lai_record = read_next_lai_record(
        afp,
        ahdr,
        ghdr,
        arec,
        n_samples,
        n_ancestries,
        n_words,
        &an1_arr,
        &nan1_arr,
        &an2_arr,
        &nan2_arr,
        warned_missing_flare_contigs,
        lai_record
    );
    int last_flare_geno_rid = -1;
    int64_t last_flare_pos = 0;

    while (read_next_record(gfp, ghdr, grec) == 0) {
        const char* g_chr = bcf_hdr_id2name(ghdr, grec->rid);
        int64_t g_pos = static_cast<int64_t>(grec->pos) + 1;

        while (has_lai_record) {
            bool current_block_covers =
                block.active &&
                block.geno_rid == grec->rid &&
                g_pos <= block.end_pos;

            if (lai_record.geno_rid > grec->rid ||
                (lai_record.geno_rid == grec->rid && current_block_covers)) {
                break;
            }

            LaiRecord interval_record = std::move(lai_record);
            has_lai_record = read_next_lai_record(
                afp,
                ahdr,
                ghdr,
                arec,
                n_samples,
                n_ancestries,
                n_words,
                &an1_arr,
                &nan1_arr,
                &an2_arr,
                &nan2_arr,
                warned_missing_flare_contigs,
                lai_record
            );

            while (has_lai_record &&
                   lai_record.geno_rid == interval_record.geno_rid &&
                   lai_record.pos == interval_record.pos) {
                interval_record = std::move(lai_record);
                has_lai_record = read_next_lai_record(
                    afp,
                    ahdr,
                    ghdr,
                    arec,
                    n_samples,
                    n_ancestries,
                    n_words,
                    &an1_arr,
                    &nan1_arr,
                    &an2_arr,
                    &nan2_arr,
                    warned_missing_flare_contigs,
                    lai_record
                );
            }

            if (last_flare_geno_rid > interval_record.geno_rid) {
                die("FLARE records are not sorted by genotype header contig order");
            }

            int64_t interval_start = 1;
            if (last_flare_geno_rid == interval_record.geno_rid) {
                if (interval_record.pos <= last_flare_pos) {
                    die(
                        "FLARE records are not strictly increasing at %s:%lld",
                        interval_record.chr.c_str(),
                        static_cast<long long>(interval_record.pos)
                    );
                }
                interval_start = last_flare_pos + 1;
            }

            update_open_block_from_flare(
                anc_fp,
                anc_bin.c_str(),
                anc_mks_fp,
                anc_mks.c_str(),
                anc_idx_fp,
                block,
                n_blocks_written,
                n_ancestries,
                n_words,
                interval_record.geno_rid,
                interval_record.chr.c_str(),
                interval_start,
                interval_record.pos,
                interval_record.state
            );

            last_flare_geno_rid = interval_record.geno_rid;
            last_flare_pos = interval_record.pos;
        }

        if (block.active && block.geno_rid == grec->rid && g_pos > block.end_pos &&
            (!has_lai_record || lai_record.geno_rid > grec->rid)) {
            block.end_pos = g_pos;
        }

        if (!block.active || block.geno_rid != grec->rid || g_pos > block.end_pos) {
            continue;
        }

        if (grec->n_allele < 2) {
            continue;
        }

        process_genotypes(
            ghdr,
            grec,
            block,
            n_samples,
            &gt_arr,
            &ngt_arr,
            carriers_by_alt
        );

        if (carriers_by_alt.empty()) {
            continue;
        }

        const char* raw_id = grec->d.id;
        const char* ref = grec->d.allele[0];

        for (int alt_idx = 1; alt_idx < grec->n_allele; ++alt_idx) {
            if (global_variant_index == std::numeric_limits<uint32_t>::max()) {
                die("global split variant index exceeds uint32_t limit");
            }

            const char* alt = grec->d.allele[alt_idx];
            const auto& carriers = carriers_by_alt[alt_idx];
            uint32_t mac = static_cast<uint32_t>(carriers.size());
            std::string split_id = make_split_id(raw_id, g_chr, g_pos, ref, alt);

            if (mac <= static_cast<uint32_t>(rare_threshold)) {
                uint64_t carrier_offset = tell_or_die(rare_fp, rare_bin.c_str());

                for (const CarrierHap& carrier : carriers) {
                    RareCarrierPacked packed{
                        global_variant_index,
                        pack_anc_hap(carrier.ancestry, carrier.hap_id)
                    };

                    size_t wrote = std::fwrite(
                        &packed,
                        sizeof(RareCarrierPacked),
                        1,
                        rare_fp
                    );

                    if (wrote != 1) {
                        die("failed writing rare carrier for split variant %u", global_variant_index);
                    }
                }

                uint64_t rare_mks_offset = tell_or_die(rare_mks_fp, rare_mks.c_str());
                write_rare_mks_record(
                    rare_mks_fp,
                    rare_index,
                    global_variant_index,
                    g_chr,
                    g_pos,
                    split_id,
                    ref,
                    alt,
                    static_cast<uint32_t>(alt_idx),
                    carrier_offset,
                    mac,
                    mac
                );
                write_rare_offset_idx_record(
                    rare_idx_fp,
                    rare_index,
                    global_variant_index,
                    rare_mks_offset,
                    carrier_offset,
                    mac
                );

                ++rare_index;
            } else {
                std::vector<uint64_t> bits(static_cast<size_t>(n_words), 0);

                for (const CarrierHap& carrier : carriers) {
                    set_bit(bits, carrier.hap_id);
                }

                uint64_t geno_offset = tell_or_die(common_fp, common_bin.c_str());
                size_t wrote = std::fwrite(
                    bits.data(),
                    sizeof(uint64_t),
                    static_cast<size_t>(n_words),
                    common_fp
                );

                if (wrote != static_cast<size_t>(n_words)) {
                    die("failed writing common genotype for split variant %u", global_variant_index);
                }

                uint64_t common_mks_offset = tell_or_die(common_mks_fp, common_mks.c_str());
                write_common_mks_record(
                    common_mks_fp,
                    common_index,
                    global_variant_index,
                    g_chr,
                    g_pos,
                    split_id,
                    ref,
                    alt,
                    static_cast<uint32_t>(alt_idx),
                    block.block_id,
                    geno_offset,
                    mac
                );
                write_common_offset_idx_record(
                    common_idx_fp,
                    common_index,
                    global_variant_index,
                    common_mks_offset,
                    geno_offset
                );

                ++common_index;
            }

            ++global_variant_index;
        }
    }

    close_open_block(
        anc_fp,
        anc_bin.c_str(),
        anc_mks_fp,
        anc_mks.c_str(),
        anc_idx_fp,
        block,
        n_blocks_written,
        n_ancestries,
        n_words
    );

    if (gt_arr) std::free(gt_arr);
    if (an1_arr) std::free(an1_arr);
    if (an2_arr) std::free(an2_arr);

    bcf_destroy(grec);
    bcf_destroy(arec);
    bcf_hdr_destroy(ghdr);
    bcf_hdr_destroy(ahdr);
    bcf_close(gfp);
    bcf_close(afp);

    std::fclose(common_fp);
    std::fclose(common_mks_fp);
    std::fclose(common_idx_fp);
    std::fclose(rare_fp);
    std::fclose(rare_mks_fp);
    std::fclose(rare_idx_fp);
    std::fclose(anc_fp);
    std::fclose(anc_mks_fp);
    std::fclose(anc_idx_fp);

    std::fprintf(stderr, "Finished.\n");
    std::fprintf(stderr, "Global split variants: %u\n", global_variant_index);
    std::fprintf(stderr, "Common variants:       %llu\n", static_cast<unsigned long long>(common_index));
    std::fprintf(stderr, "Rare variants:         %llu\n", static_cast<unsigned long long>(rare_index));
    std::fprintf(stderr, "Ancestry blocks:       %u\n", n_blocks_written);

    return 0;
}
