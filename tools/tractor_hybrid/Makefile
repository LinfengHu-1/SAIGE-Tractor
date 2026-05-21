# designed by Kai, implemented by codex

CXX ?= g++
PKG_CONFIG ?= pkg-config

BIN_DIR ?= bin
FLARE_TARGET := $(BIN_DIR)/flare_subset_to_tractor_hybrid
ESTIMATE_TARGET := $(BIN_DIR)/estimate_mac_threshold

PKG_HTSLIB_CFLAGS := $(shell $(PKG_CONFIG) --cflags htslib 2>/dev/null)
PKG_HTSLIB_LIBS := $(shell $(PKG_CONFIG) --libs htslib 2>/dev/null)

HTSLIB_CFLAGS ?= $(PKG_HTSLIB_CFLAGS)
HTSLIB_LIBS ?= $(PKG_HTSLIB_LIBS)
ifeq ($(strip $(HTSLIB_LIBS)),)
HTSLIB_LIBS := -lhts
endif

CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
CPPFLAGS += $(HTSLIB_CFLAGS)
LDLIBS += $(HTSLIB_LIBS)

.PHONY: all clean check-deps

all: check-deps $(FLARE_TARGET) $(ESTIMATE_TARGET)

check-deps:
	@printf '#include <htslib/hts.h>\n#include <htslib/vcf.h>\n' | \
		$(CXX) $(CPPFLAGS) -x c++ -E - >/dev/null 2>&1 || \
		( echo "ERROR: htslib headers not found. Install htslib or set HTSLIB_CFLAGS/HTSLIB_LIBS." >&2; exit 1 )

$(FLARE_TARGET): src/flare_subset_to_tractor_hybrid.cpp
	mkdir -p $(BIN_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

$(ESTIMATE_TARGET): src/estimate_mac_threshold.cpp
	mkdir -p $(BIN_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

clean:
	rm -rf $(BIN_DIR)
