# Makefile -- riordan-graph-enumerator
#
# Targets
#   make            build every tool into build/ (C++20 compiler, C compiler for the bundled nauty)
#   make check      self-tests, about 15 s: labeled enumeration, recognizer, u(1..8) by methods
#                   A, B and C, the listing tool and the sampler
#   make clean      remove build/
#
# Variables (override on the command line)
#   CXX, CC         compilers, e.g. make CXX=clang++ CC=clang on macOS
#   CXXFLAGS        optimisation flags (default -O3), e.g. CXXFLAGS="-O3 -march=native"
#
# nauty is compiled with MAXN=WORDSIZE=64 (one 64-bit word per adjacency row) and thread-local
# storage (USE_TLS), so that canonical labelling can run in several threads. geng, used only by
# method C, is compiled separately with WORDSIZE=32.
CXX      ?= g++
CC       ?= gcc
CXXFLAGS ?= -O3
NAUTYDIR := third_party/nauty
NF       := -DMAXN=WORDSIZE -DWORDSIZE=64
B        := build
NAUTY_SRC := nauty nautil naugraph schreier naurng nautinv
NAUTY_OBJ := $(NAUTY_SRC:%=$(B)/nauty/%.o)
CXXALL   := -std=c++20 $(CXXFLAGS) $(NF) -Iinclude -I$(NAUTYDIR) -pthread
HDR      := $(wildcard include/rgraph/*.hpp)

TOOLS := rgcount rglist rgsample tree_count recognize_all bell structure complement_labeled complement_unlabeled labelstats labelstats2 \
         quick_check recognize_test geng

all: $(TOOLS:%=$(B)/%)

$(B)/nauty/%.o: $(NAUTYDIR)/%.c
	@mkdir -p $(B)/nauty
	$(CC) -O3 $(NF) -DUSE_TLS -c $< -o $@

$(B)/rgcount: tools/rgcount.cpp $(HDR) $(NAUTY_OBJ)
	$(CXX) $(CXXALL) $< $(NAUTY_OBJ) -o $@
$(B)/rglist: tools/rglist.cpp $(HDR) $(NAUTY_OBJ)
	$(CXX) $(CXXALL) $< $(NAUTY_OBJ) -o $@
$(B)/rgsample: tools/rgsample.cpp $(HDR) $(NAUTY_OBJ)
	$(CXX) $(CXXALL) $< $(NAUTY_OBJ) -o $@
$(B)/tree_count: verify/tree_count.cpp $(NAUTY_OBJ)
	$(CXX) $(CXXALL) $< $(NAUTY_OBJ) -o $@
$(B)/recognize_all: verify/recognize_all.cpp
	$(CXX) -std=c++20 $(CXXFLAGS) $< -o $@
$(B)/bell: analysis/bell.cpp
	$(CXX) -std=c++20 $(CXXFLAGS) -pthread $< -o $@
$(B)/structure: analysis/structure.cpp $(HDR) $(NAUTY_OBJ)
	$(CXX) $(CXXALL) $< $(NAUTY_OBJ) -o $@
$(B)/complement_labeled: analysis/complement_labeled.cpp $(HDR) $(NAUTY_OBJ)
	$(CXX) $(CXXALL) $< $(NAUTY_OBJ) -o $@
$(B)/complement_unlabeled: analysis/complement_unlabeled.cpp $(HDR) $(NAUTY_OBJ)
	$(CXX) $(CXXALL) $< $(NAUTY_OBJ) -o $@
$(B)/labelstats: analysis/labelstats.cpp $(HDR) $(NAUTY_OBJ)
	$(CXX) $(CXXALL) $< $(NAUTY_OBJ) -o $@
$(B)/labelstats2: analysis/labelstats2.cpp $(HDR) $(NAUTY_OBJ)
	$(CXX) $(CXXALL) $< $(NAUTY_OBJ) -o $@
$(B)/quick_check: tests/quick_check.cpp $(HDR) $(NAUTY_OBJ)
	$(CXX) $(CXXALL) $< $(NAUTY_OBJ) -o $@
$(B)/recognize_test: tests/recognize_test.cpp $(HDR)
	$(CXX) -std=c++20 $(CXXFLAGS) -Iinclude $< -o $@
# nauty's geng (all graphs of order n), used by method C
$(B)/geng: $(NAUTYDIR)/geng.c $(NAUTYDIR)/gtools.c
	@mkdir -p $(B)
	$(CC) -O3 -DMAXN=WORDSIZE -DWORDSIZE=32 -o $@ $(NAUTYDIR)/geng.c $(NAUTYDIR)/gtools.c \
	  $(NAUTYDIR)/nauty.c $(NAUTYDIR)/nautil.c $(NAUTYDIR)/naugraph.c $(NAUTYDIR)/schreier.c $(NAUTYDIR)/naurng.c

check: all
	./$(B)/quick_check
	./$(B)/recognize_test
	sh tests/small_counts.sh
	sh tests/rglist_test.sh
	sh tests/rgsample_test.sh

clean:
	rm -rf $(B)

.PHONY: all check clean
