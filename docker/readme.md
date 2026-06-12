# SAIGE-Tractor Docker image

Builds a single image with **file conversion** (Kai's `tractor_hybrid` tools), **GWAS**
(step1/step2/step3 + createSparseGRM), and **cross-ancestry rg/h2** (no-GRM null + RHE) entry points.
See `docker/Dockerfile` for the layer-by-layer rationale.

## Build

Run from the package dir (`SAIGE-Tractor/`); context is `.`, `.dockerignore` excludes `.pixi`, `.git`,
the mac `plink2_includes.a`, etc.

```bash
cd SAIGE-Tractor
DOCKER_BUILDKIT=1 docker buildx build \
  --platform linux/amd64 \
  -f docker/Dockerfile \
  -t saigetractor:rg-h2 \
  --load .
```

- `--platform linux/amd64` — target x86_64 clusters (runs locally on Apple Silicon via Colima+Rosetta).
- `docker buildx` (not plain `docker build`) — **required**: the Dockerfile uses `RUN --mount=type=cache`
  (ccache) and `# syntax=docker/dockerfile:1`, which need BuildKit/buildx.
- `--load` — load the finished image into the local docker image store (so `docker images` /
  `docker save` see it).

The first build is slow (emulated amd64 + a full pixi conda env). Re-running is fast: layers are
ordered so the expensive converter + pixi env + plink2 layers are cached, and a package-code edit only
re-runs `COPY → R CMD INSTALL` (with **ccache**, so unchanged `.cpp` are cache hits).

### One-time host setup (macOS, done once / per new machine)
```bash
brew install colima docker docker-buildx
mkdir -p ~/.docker/cli-plugins
ln -sfn "$(brew --prefix)/opt/docker-buildx/bin/docker-buildx" ~/.docker/cli-plugins/docker-buildx
colima start --cpu 4 --memory 8 --disk 60 --vm-type vz --vz-rosetta   # also: once per reboot
```

## Entry points (all on PATH in the image)
- **File conversion:** `flare_subset_to_tractor_hybrid`, `estimate_mac_threshold`
- **GWAS:** `step1_fitNULLGLMM.R`, `step2_SPAtests.R` (also integrated `--estimate_cross_anc_rg`),
  `step3_LDmat.R`, `createSparseGRM.R`, `step0_extractUnrelatedFromGRM.R`
- **cross-ancestry rg/h2:** `step1_fitNULL_noGRM.R`, `step2_rgRHEonly.R`, `step3_crossAncestryRgRHE.R`

```bash
docker run --rm -v "$PWD:/work" -w /work saigetractor:rg-h2 step1_fitNULLGLMM.R --help
docker run --rm -v "$PWD:/work" -w /work saigetractor:rg-h2 \
  flare_subset_to_tractor_hybrid geno.phased.vcf.gz flare.anc.vcf.gz 2 10 out_prefix
```

## Ship to the cluster as Singularity (no Docker needed on the cluster)
```bash
docker save saigetractor:rg-h2 | gzip > st.tar.gz        # on a machine with docker
# scp st.tar.gz to the cluster, then on the cluster (Apptainer/Singularity present):
gunzip st.tar.gz && apptainer build st.sif docker-archive://st.tar
apptainer exec --bind "$PWD:/work" st.sif step2_rgRHEonly.R --help
```

## Publish to Docker Hub (optional, legacy flow)
```bash
docker tag  saigetractor:rg-h2 wzhou88/saigetractor:rg-h2
docker push wzhou88/saigetractor:rg-h2
```
