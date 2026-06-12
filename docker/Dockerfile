# syntax=docker/dockerfile:1
FROM ubuntu:20.04

LABEL maintainer="wzhou@broadinstitute.org"

ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8
ENV OMP_NUM_THREADS=1
# Force step_2 to use 1 single thread. More threads are ineffective
ENV PIXI_NO_PATH_UPDATE=1

# Build (linux/amd64 for the cluster). Context = this package dir (SAIGE-Tractor); see .dockerignore.
# REQUIRES BuildKit (for the ccache `RUN --mount=type=cache`); prefix DOCKER_BUILDKIT=1:
#   DOCKER_BUILDKIT=1 docker build --platform linux/amd64 -f docker/Dockerfile -t saigetractor:rg-h2 .
# Convert to Singularity ON THE CLUSTER (no docker needed there):
#   docker save saigetractor:rg-h2 | gzip > st.tar.gz          # on a machine with docker
#   # scp st.tar.gz to cluster, then:
#   gunzip st.tar.gz && apptainer build st.sif docker-archive://st.tar
#
# LAYER ORDERING IS DELIBERATE (fast incremental rebuilds): the expensive, code-INDEPENDENT layers
# (converter, pixi conda env, plink2 archive) come BEFORE `COPY . .`, so editing package code only
# re-runs the cheap tail (R CMD INSTALL + script exposure). The multi-GB pixi env is NOT re-downloaded
# unless pixi.toml / pixi.lock change. So a typical `R/`-or-`src/`-edit rebuild is minutes, not ~hour.

RUN apt-get update -qq && \
    DEBIAN_FRONTEND=noninteractive apt-get -y install --no-install-recommends \
      curl \
      ca-certificates \
      ccache \
      g++ \
      libhts-dev \
      make \
      pkg-config \
      time \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# --- Layer A: Kai's tractor_hybrid converter (file conversion). Independent of R/pixi, built against
#     system htslib so it's a standalone binary. Cached unless tools/tractor_hybrid changes. ---
COPY tools/tractor_hybrid /app/tools/tractor_hybrid
RUN mkdir -p /scripts && \
   cp tools/tractor_hybrid/Makefile /scripts/Makefile && \
   cp -R tools/tractor_hybrid/src /scripts/src && \
   make -C /scripts CXX=g++ CXXFLAGS="-O3 -std=c++17 -Wall -Wextra -Wpedantic" BIN_DIR=/scripts/bin && \
   /scripts/bin/flare_subset_to_tractor_hybrid 2>&1 | grep -q "Usage" && \
   /scripts/bin/estimate_mac_threshold 2>&1 | grep -q "Usage"

# --- Layer B: pixi conda env (R + RcppArmadillo + MetaSKAT + openblas + htslib + gxx_linux-64 ...) +
#     CRAN-only lintools. THE EXPENSIVE LAYER. Copy ONLY the manifest first so it is reused across
#     package-code changes; re-runs ONLY when pixi.toml / pixi.lock / activate.sh change. ---
COPY pixi.toml pixi.lock activate.sh ./
# pixi downloads are resumable/cached, but a single `tls handshake eof` from conda.anaconda.org (common
# over emulated amd64 networking) aborts the whole solve. Retry the install: each attempt resumes from
# the cache and makes progress; the trailing `pixi install` outside the loop fails the build for real if
# all retries are exhausted. Lower concurrency also eases the TLS pressure.
RUN curl -fsSL https://pixi.sh/install.sh | sh && \
    mv /root/.pixi/bin/pixi /bin && \
    (pixi config set --global concurrent-downloads 2 || true) && \
    for i in 1 2 3 4 5 6 7 8; do pixi install && break || { echo "== pixi install attempt $i failed; retrying =="; sleep 15; }; done && \
    pixi install && \
    rm -rf /root/.cache && \
    pixi run install-github-deps

# --- Layer C: bundled plink2 STATIC archive at the package root (/app/plink2_includes.a). Depends
#     only on the pixi env (the task downloads plink-ng itself), not on package source. Built with the
#     pixi $CXX (gxx_linux-64) + `ar rcs` exactly as src/Makevars expects (`PKG_LIBS += ../plink2_includes.a`).
#     Leaves /app/plink-ng (needed by Makevars -I) and /app/plink2_includes.a. ---
RUN pixi run install-plink2

# --- Layer D: NOW copy the package source. A change to R/, src/, DESCRIPTION, etc. invalidates from
#     HERE down only; layers A-C stay cached. .dockerignore keeps .pixi, plink-ng, and the mac
#     plink2_includes.a out of the context so this COPY doesn't clobber the artifacts built above. ---
COPY . .

# compile + install the R package (GWAS *and* cross-ancestry rg/h2 RHE both live here).
# This is the main layer that re-runs on a code edit. ccache (persisted via the BuildKit cache mount)
# turns unchanged .cpp into cache hits -> editing one file recompiles ONLY that file (seconds), not all
# of src/. We point R's CC/CXX at `ccache <conda compiler>` via a personal ~/.R/Makevars, which R reads
# AFTER (and so overrides) the Makevars.site that activate.sh writes. Absolute /usr/bin/ccache avoids
# any PATH surprises under `pixi run`.
RUN --mount=type=cache,target=/ccache \
    pixi run bash -c '\
      CCBIN=$(ls $CONDA_PREFIX/bin/*-gcc $CONDA_PREFIX/bin/*-cc 2>/dev/null | head -1); \
      mkdir -p /root/.R; \
      { echo "CC = /usr/bin/ccache $CCBIN"; \
        echo "CXX = /usr/bin/ccache $CXX"; \
        echo "CXX11 = /usr/bin/ccache $CXX"; \
        echo "CXX14 = /usr/bin/ccache $CXX"; \
        echo "CXX17 = /usr/bin/ccache $CXX"; } > /root/.R/Makevars; \
      export CCACHE_DIR=/ccache CCACHE_MAXSIZE=3G; \
      R CMD INSTALL . && ccache -s'

ENV PATH=/scripts/bin:/usr/local/bin:/usr/local/sbin:/usr/sbin:/usr/bin:/sbin:/bin

# Expose ALL CLI entry points on PATH. Every script runs through the pixi-env Rscript, so normalize
# the shebang to the pixi-run form (idempotent for scripts that already have it) before installing.
#   GWAS:                 step1_fitNULLGLMM, step2_SPAtests (also does integrated --estimate_cross_anc_rg),
#                         step3_LDmat, createSparseGRM, step0_extractUnrelatedFromGRM
#   cross-ancestry rg/h2: step1_fitNULL_noGRM, step2_rgRHEonly, step3_crossAncestryRgRHE
RUN set -e; \
   SHEBANG='#!/usr/bin/env -S pixi run --manifest-path /app/pixi.toml Rscript'; \
   for s in step1_fitNULLGLMM step2_SPAtests step3_LDmat createSparseGRM \
            step0_extractUnrelatedFromGRM \
            step1_fitNULL_noGRM step2_rgRHEonly step3_crossAncestryRgRHE; do \
     cp "extdata/$s.R" "/usr/local/bin/$s.R"; \
     sed -i "1s|^#!.*|$SHEBANG|" "/usr/local/bin/$s.R"; \
     chmod a+x "/usr/local/bin/$s.R"; \
   done; \
   for s in step1_fitNULLGLMM step2_SPAtests step3_LDmat createSparseGRM \
            step1_fitNULL_noGRM step2_rgRHEonly step3_crossAncestryRgRHE; do \
     echo "=== $s.R --help ==="; \
     "$s.R" --help >/dev/null 2>&1 || echo "WARNING: $s.R --help returned nonzero"; \
   done

# Final assertion: BOTH tractor_hybrid converter binaries are callable by BARE NAME via PATH
# (not just by absolute path), so the build fails here if /scripts/bin isn't wired up.
RUN set -e; \
   command -v flare_subset_to_tractor_hybrid; \
   command -v estimate_mac_threshold; \
   flare_subset_to_tractor_hybrid 2>&1 | grep -q "Usage"; \
   estimate_mac_threshold 2>&1 | grep -q "Usage"; \
   echo "OK: both tractor_hybrid converter tools are on PATH and runnable"

CMD ["bash"]
