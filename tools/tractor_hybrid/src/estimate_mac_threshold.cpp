// estimate_mac_threshold.cpp
// designed by Kai, implemented by codex
//
// Estimate a sparse/dense MAC threshold for the ancestry-aware packed backend.

#include <htslib/hts.h>
#include <htslib/vcf.h>

#include <algorithm>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

struct Options {
    const char* genotype_path = nullptr;
    double storage_weight = 1.0;
    double query_weight = 0.0;
    double dense_word_cost = 1.0;
    double sparse_carrier_cost = 1.0;
    uint64_t max_records = 0;
};

struct ScanStats {
    int n_samples = 0;
    uint64_t n_haps = 0;
    uint64_t n_words = 0;
    uint64_t records_seen = 0;
    uint64_t split_variants = 0;
    uint64_t total_alt_carriers = 0;
    uint32_t max_mac = 0;
    std::map<uint32_t, uint64_t> mac_counts;
};

struct ThresholdScore {
    uint32_t threshold = 0;
    uint64_t rare_variants = 0;
    uint64_t common_variants = 0;
    uint64_t rare_carriers = 0;
    uint64_t common_dense_words = 0;
    long double payload_bytes = 0.0L;
    long double query_units = 0.0L;
    long double score = 0.0L;
};

[[noreturn]] static void die(const char* fmt, ...) {
    va_list args;
    va_start(args, fmt);
    std::fputs("ERROR: ", stderr);
    std::vfprintf(stderr, fmt, args);
    std::fputc('\n', stderr);
    va_end(args);
    std::exit(1);
}

static uint64_t parse_u64_arg(const char* value, const char* name) {
    char* end = nullptr;
    unsigned long long parsed = std::strtoull(value, &end, 10);
    if (!end || *end != '\0') die("invalid %s: %s", name, value);
    return static_cast<uint64_t>(parsed);
}

static double parse_double_arg(const char* value, const char* name) {
    char* end = nullptr;
    double parsed = std::strtod(value, &end);
    if (!end || *end != '\0' || parsed < 0.0) die("invalid %s: %s", name, value);
    return parsed;
}

static void print_usage(const char* prog) {
    std::fprintf(
        stderr,
        "Usage:\n"
        "  %s genotype.phased.vcf.gz [options]\n\n"
        "Options:\n"
        "  --storage-weight FLOAT       Weight for payload bytes (default: 1.0)\n"
        "  --query-weight FLOAT         Weight for per-query work units (default: 0.0)\n"
        "  --dense-word-cost FLOAT      Cost per dense uint64 word scanned (default: 1.0)\n"
        "  --sparse-carrier-cost FLOAT  Cost per sparse carrier scanned (default: 1.0)\n"
        "  --max-records N              Scan only the first N VCF records (default: all)\n",
        prog
    );
}

static Options parse_options(int argc, char** argv) {
    Options opt;
    if (argc < 2) {
        print_usage(argv[0]);
        std::exit(1);
    }

    opt.genotype_path = argv[1];

    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        auto require_value = [&](const char* name) -> const char* {
            if (i + 1 >= argc) die("%s requires a value", name);
            return argv[++i];
        };

        if (arg == "--storage-weight") {
            opt.storage_weight = parse_double_arg(require_value("--storage-weight"), "--storage-weight");
        } else if (arg == "--query-weight") {
            opt.query_weight = parse_double_arg(require_value("--query-weight"), "--query-weight");
        } else if (arg == "--dense-word-cost") {
            opt.dense_word_cost = parse_double_arg(require_value("--dense-word-cost"), "--dense-word-cost");
        } else if (arg == "--sparse-carrier-cost") {
            opt.sparse_carrier_cost = parse_double_arg(require_value("--sparse-carrier-cost"), "--sparse-carrier-cost");
        } else if (arg == "--max-records") {
            opt.max_records = parse_u64_arg(require_value("--max-records"), "--max-records");
        } else if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            std::exit(0);
        } else {
            die("unknown option: %s", arg.c_str());
        }
    }

    if (opt.storage_weight == 0.0 && opt.query_weight == 0.0) {
        die("at least one of --storage-weight or --query-weight must be positive");
    }

    return opt;
}

static void validate_extra_ploidy(const int32_t* sample_gt, int ploidy, const char* chr, int64_t pos, int sample_index) {
    for (int p = 2; p < ploidy; ++p) {
        int32_t gt = sample_gt[p];
        if (gt == bcf_int32_vector_end) break;
        die("non-diploid genotype at %s:%lld sample index %d", chr, static_cast<long long>(pos), sample_index);
    }
}

static ScanStats scan_mac_distribution(const char* path, uint64_t max_records) {
    htsFile* fp = bcf_open(path, "r");
    if (!fp) die("cannot open genotype VCF/BCF: %s", path);

    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    if (!hdr) die("cannot read genotype header: %s", path);

    ScanStats stats;
    stats.n_samples = bcf_hdr_nsamples(hdr);
    if (stats.n_samples <= 0) die("genotype VCF has no samples");
    stats.n_haps = static_cast<uint64_t>(stats.n_samples) * 2ULL;
    stats.n_words = (stats.n_haps + 63ULL) / 64ULL;

    bcf1_t* rec = bcf_init();
    int32_t* gt_arr = nullptr;
    int ngt_arr = 0;

    while (bcf_read(fp, hdr, rec) == 0) {
        if (max_records != 0 && stats.records_seen >= max_records) break;
        ++stats.records_seen;

        bcf_unpack(rec, BCF_UN_STR);
        if (rec->n_allele < 2) continue;

        const char* chr = bcf_hdr_id2name(hdr, rec->rid);
        int64_t pos = static_cast<int64_t>(rec->pos) + 1;

        int ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
        if (ngt <= 0) die("missing FORMAT/GT at %s:%lld", chr, static_cast<long long>(pos));
        if (ngt % stats.n_samples != 0) {
            die("GT field length is not divisible by sample count at %s:%lld", chr, static_cast<long long>(pos));
        }

        int ploidy = ngt / stats.n_samples;
        if (ploidy < 2) die("expected diploid GT at %s:%lld", chr, static_cast<long long>(pos));

        std::vector<uint32_t> mac_by_alt(static_cast<size_t>(rec->n_allele), 0);

        for (int i = 0; i < stats.n_samples; ++i) {
            const int32_t* sample_gt = gt_arr + static_cast<size_t>(i) * ploidy;
            validate_extra_ploidy(sample_gt, ploidy, chr, pos, i);

            int32_t g0 = sample_gt[0];
            int32_t g1 = sample_gt[1];
            if (g0 == bcf_int32_vector_end || g1 == bcf_int32_vector_end ||
                bcf_gt_is_missing(g0) || bcf_gt_is_missing(g1)) {
                die("missing genotype at %s:%lld sample index %d", chr, static_cast<long long>(pos), i);
            }
            if (!bcf_gt_is_phased(g1)) {
                die("unphased genotype at %s:%lld sample index %d", chr, static_cast<long long>(pos), i);
            }

            int allele0 = bcf_gt_allele(g0);
            int allele1 = bcf_gt_allele(g1);
            if (allele0 > 0 && allele0 < rec->n_allele) ++mac_by_alt[allele0];
            if (allele1 > 0 && allele1 < rec->n_allele) ++mac_by_alt[allele1];
        }

        for (int alt = 1; alt < rec->n_allele; ++alt) {
            uint32_t mac = mac_by_alt[alt];
            ++stats.mac_counts[mac];
            ++stats.split_variants;
            stats.total_alt_carriers += mac;
            stats.max_mac = std::max(stats.max_mac, mac);
        }
    }

    if (gt_arr) std::free(gt_arr);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(fp);

    return stats;
}

static ThresholdScore score_threshold(const ScanStats& stats, const Options& opt, uint32_t threshold) {
    ThresholdScore result;
    result.threshold = threshold;

    for (const auto& [mac, count] : stats.mac_counts) {
        if (mac <= threshold) {
            result.rare_variants += count;
            result.rare_carriers += static_cast<uint64_t>(mac) * count;
        } else {
            result.common_variants += count;
        }
    }

    result.common_dense_words = result.common_variants * stats.n_words;
    result.payload_bytes =
        static_cast<long double>(result.rare_carriers) * 8.0L +
        static_cast<long double>(result.common_dense_words) * 8.0L;
    result.query_units =
        static_cast<long double>(result.rare_carriers) * opt.sparse_carrier_cost +
        static_cast<long double>(result.common_dense_words) * opt.dense_word_cost;
    result.score =
        static_cast<long double>(opt.storage_weight) * result.payload_bytes +
        static_cast<long double>(opt.query_weight) * result.query_units;

    return result;
}

static ThresholdScore find_best_threshold(const ScanStats& stats, const Options& opt) {
    ThresholdScore best = score_threshold(stats, opt, 0);

    for (const auto& [mac, _count] : stats.mac_counts) {
        ThresholdScore candidate = score_threshold(stats, opt, mac);
        if (candidate.score < best.score ||
            (candidate.score == best.score && candidate.threshold > best.threshold)) {
            best = candidate;
        }
    }

    return best;
}

static std::string format_bytes(long double bytes) {
    const char* units[] = {"B", "KiB", "MiB", "GiB", "TiB", "PiB"};
    int unit = 0;
    while (bytes >= 1024.0L && unit < 5) {
        bytes /= 1024.0L;
        ++unit;
    }

    char buf[128];
    std::snprintf(buf, sizeof(buf), "%.2Lf %s", bytes, units[unit]);
    return buf;
}

static uint32_t model_break_even_threshold(const ScanStats& stats, const Options& opt) {
    long double dense_cost =
        static_cast<long double>(opt.storage_weight) * 8.0L * stats.n_words +
        static_cast<long double>(opt.query_weight) * opt.dense_word_cost * stats.n_words;
    long double sparse_per_mac =
        static_cast<long double>(opt.storage_weight) * 8.0L +
        static_cast<long double>(opt.query_weight) * opt.sparse_carrier_cost;

    if (sparse_per_mac <= 0.0L) return 0;
    long double threshold = dense_cost / sparse_per_mac;
    if (threshold > static_cast<long double>(UINT32_MAX)) return UINT32_MAX;
    return static_cast<uint32_t>(threshold);
}

static void print_score_row(const char* label, const ThresholdScore& score) {
    std::cout << std::left << std::setw(14) << label
              << std::right << std::setw(10) << score.threshold
              << std::setw(14) << score.rare_variants
              << std::setw(14) << score.common_variants
              << std::setw(18) << score.rare_carriers
              << std::setw(16) << format_bytes(score.payload_bytes)
              << std::setw(18) << std::fixed << std::setprecision(0) << score.query_units
              << "\n";
}

int main(int argc, char** argv) {
    Options opt = parse_options(argc, argv);
    ScanStats stats = scan_mac_distribution(opt.genotype_path, opt.max_records);

    if (stats.split_variants == 0) die("no split ALT variants found");

    uint32_t break_even = model_break_even_threshold(stats, opt);
    ThresholdScore best = find_best_threshold(stats, opt);

    std::cout << "MAC threshold estimator\n";
    std::cout << "Samples:              " << stats.n_samples << "\n";
    std::cout << "Haplotypes:           " << stats.n_haps << "\n";
    std::cout << "Dense words/variant:  " << stats.n_words << "\n";
    std::cout << "Dense bytes/variant:  " << stats.n_words * 8ULL << "\n";
    std::cout << "VCF records scanned:  " << stats.records_seen << "\n";
    std::cout << "Split ALT variants:   " << stats.split_variants << "\n";
    std::cout << "Max MAC observed:     " << stats.max_mac << "\n";
    std::cout << "Mean MAC observed:    "
              << std::fixed << std::setprecision(2)
              << static_cast<double>(stats.total_alt_carriers) / static_cast<double>(stats.split_variants)
              << "\n";
    std::cout << "Model break-even MAC: " << break_even << "\n";
    std::cout << "Recommended threshold:" << best.threshold << "\n\n";

    std::cout << std::left << std::setw(14) << "candidate"
              << std::right << std::setw(10) << "threshold"
              << std::setw(14) << "rare_vars"
              << std::setw(14) << "common_vars"
              << std::setw(18) << "rare_carriers"
              << std::setw(16) << "payload"
              << std::setw(18) << "query_units"
              << "\n";

    print_score_row("best", best);
    print_score_row("break_even", score_threshold(stats, opt, break_even));

    std::vector<uint32_t> candidates = {64, 128, 256, 512, 1024, 2048, 4096, 8192};
    for (uint32_t candidate : candidates) {
        if (candidate <= stats.max_mac || candidate <= break_even) {
            print_score_row(("T=" + std::to_string(candidate)).c_str(), score_threshold(stats, opt, candidate));
        }
    }

    return 0;
}
