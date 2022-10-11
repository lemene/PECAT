
BUILD_BIN_DIR = ../build/bin
BUILD_OBJ_DIR = ../build/obj
BUILD_DIR += $(BUILD_BIN_DIR) $(BUILD_OBJ_DIR)

BUILD_DIR += $(BUILD_OBJ_DIR)/prog $(BUILD_OBJ_DIR)/../thirdparty/edlib/edlib/src $(BUILD_OBJ_DIR)/correct \
             $(BUILD_OBJ_DIR)/bridge $(BUILD_OBJ_DIR)/assemble $(BUILD_OBJ_DIR)/correct/align $(BUILD_OBJ_DIR)/kmer \
			 $(BUILD_OBJ_DIR)/utils $(BUILD_OBJ_DIR)/phase $(BUILD_OBJ_DIR)/assemble/simplify

simplify_src = assemble/simplify/simplifier.cpp ./assemble/simplify/spur.cpp ./assemble/simplify/bridge.cpp \
	./assemble/simplify/transitive.cpp ./assemble/simplify/cross.cpp ./assemble/simplify/low_quality.cpp \
	./assemble/simplify/unreliable.cpp ./assemble/simplify/loop.cpp ./assemble/simplify/bubble.cpp \
	./assemble/simplify/best.cpp ./assemble/simplify/semi_bubble.cpp ./assemble/simplify/repeat.cpp \
	./assemble/simplify/phase.cpp ./assemble/simplify/pbridge.cpp ./assemble/simplify/pspur.cpp

libfsa_src = ./utils/argument_parser.cpp ./utils/logger.cpp overlap.cpp read_store.cpp sequence.cpp\
             utility.cpp file_io.cpp seq_io.cpp overlap_store.cpp graph.cpp\
			 assemble/overlap_assemble.cpp simple_align.cpp assemble/string_graph.cpp \
			 assemble/read_variants.cpp assemble/string_node.cpp assemble/string_edge.cpp \
			 assemble/asm_dataset.cpp assemble/asm_options.cpp assemble/stat_readinfo.cpp\
			 ../thirdparty/edlib/edlib/src/edlib.cpp correct/aligner.cpp ./correct/dw.cpp  ./correct/diff.cpp  \
			 ./correct/align/edlib_aligner.cpp ./correct/align/ksw2_aligner.cpp \
			 ./utils/string_pool.cpp ./utils/pair_file.cpp ./utils/misc.cpp \
			 ./phase/phase_info.cpp ./phase/contig_phaser.cpp $(simplify_src)

libfsa_src += $(wildcard ./correct/align/*.cpp)

libfsa_obj:=$(patsubst %.cpp, $(BUILD_OBJ_DIR)/%.o, $(libfsa_src))

fsa_ol_assemble_src = ./prog/fsa_ol_assemble.cpp
fsa_ol_assemble_obj:=$(patsubst %.cpp, $(BUILD_OBJ_DIR)/%.o, $(fsa_ol_assemble_src))

fsa_ctg_bridge_src = ./prog/fsa_ctg_bridge.cpp bridge/contig_bridge.cpp bridge/contig_graph.cpp bridge/contig_link.cpp bridge/contig_link_store.cpp
fsa_ctg_bridge_obj:=$(patsubst %.cpp, $(BUILD_OBJ_DIR)/%.o, $(fsa_ctg_bridge_src))

fsa_rd_haplotype_src = ./prog/fsa_rd_haplotype.cpp ./phase/read_haplotype.cpp ./phase/local_phaser.cpp ./phase/utility.cpp  ./phase/phs_options.cpp ./phase/phs_dataset.cpp
fsa_rd_haplotype_obj:=$(patsubst %.cpp, $(BUILD_OBJ_DIR)/%.o, $(fsa_rd_haplotype_src))

fsa_rd_tools_src = ./prog/fsa_rd_tools.cpp read_tools.cpp
fsa_rd_tools_obj:=$(patsubst %.cpp, $(BUILD_OBJ_DIR)/%.o, $(fsa_rd_tools_src))

fsa_ol_tools_src = ./prog/fsa_ol_tools.cpp overlap_tools.cpp 
fsa_ol_tools_obj:=$(patsubst %.cpp, $(BUILD_OBJ_DIR)/%.o, $(fsa_ol_tools_src))

fsa_misc_tools_src = ./prog/fsa_misc_tools.cpp misc_tools.cpp
fsa_misc_tools_obj:=$(patsubst %.cpp, $(BUILD_OBJ_DIR)/%.o, $(fsa_misc_tools_src))

#./correct/alignment_graph.cpp ./correct/alignment_graph2.cpp
fsa_rd_correct_src = ./prog/fsa_rd_correct.cpp  \
                     ./correct/read_correct.cpp  ./correct/alignment_graph.cpp
fsa_rd_correct_obj:=$(patsubst %.cpp, $(BUILD_OBJ_DIR)/%.o, $(fsa_rd_correct_src))

# ./correct/alignment_graph.cpp ./correct/alignment_graph2.cpp
fsa_ctg_correct_src = ./prog/fsa_ctg_correct.cpp \
                     ./correct/contig_correct.cpp ./correct/alignment_graph.cpp
fsa_ctg_correct_obj := $(patsubst %.cpp, $(BUILD_OBJ_DIR)/%.o, $(fsa_ctg_correct_src))

fsa_ol_purge_src = ./prog/fsa_ol_purge.cpp ./overlap_purge.cpp
fsa_ol_purge_obj:=$(patsubst %.cpp, $(BUILD_OBJ_DIR)/%.o, $(fsa_ol_purge_src))

fsa_ol_refine_src = ./prog/fsa_ol_refine.cpp ./overlap_refine.cpp
fsa_ol_refine_obj:=$(patsubst %.cpp, $(BUILD_OBJ_DIR)/%.o, $(fsa_ol_refine_src))

fsa_kmer_bin_src = ./prog/fsa_kmer_bin.cpp ./kmer/kmer_bin.cpp
fsa_kmer_bin_obj:=$(patsubst %.cpp, $(BUILD_OBJ_DIR)/%.o, $(fsa_kmer_bin_src))

ifeq ($(BUILDDEBUG), 1)
	CFLAGS =  -Wall -O0  -D_FILE_OFFSET_BITS=64 -g -I.  -I ../thirdparty/edlib/edlib/include -I ../thirdparty/ksw2
	CXXFLAGS = -std=c++11 -Wall -O0  -D_FILE_OFFSET_BITS=64 -g -I.  -I ../thirdparty/edlib/edlib/include -I ../thirdparty/ksw2
else
	CFLAGS =  -Wall -O3  -D_FILE_OFFSET_BITS=64 -I. -I ../thirdparty/edlib/edlib/include -I ../thirdparty/ksw2
	CXXFLAGS = -std=c++11 -Wall -O3 -D_FILE_OFFSET_BITS=64 -I.  -I ../thirdparty/edlib/edlib/include -I ../thirdparty/ksw2
endif

LDFLAGS = -pthread -lz $(BUILD_BIN_DIR)/libksw2.a $


obj:=$(libfsa_obj) $(fsa_ol_assemble_obj) \
	$(fsa_ctg_bridge_obj) $(fsa_rd_correct_obj) \
	$(fsa_rd_haplotype_obj) \
	$(fsa_ol_tools_obj) $(fsa_rd_tools_obj) $(fsa_misc_tools_obj) \
	$(fsa_ol_purge_obj) $(fsa_ol_refine_obj) $(fsa_ctg_correct_obj) $(fsa_kmer_bin_obj)

dep:=$(patsubst %.o,%.d, $(obj))


bin = libfsa.a fsa_ol_assemble fsa_ctg_bridge \
      fsa_rd_haplotype fsa_ol_tools fsa_rd_tools \
	  fsa_misc_tools fsa_rd_correct fsa_ol_purge fsa_ol_refine \
	  fsa_ctg_correct fsa_kmer_bin

bin:=$(patsubst %,  $(BUILD_BIN_DIR)/%, $(bin))


all: $(BUILD_DIR) $(bin) 

ifneq ($(MAKECMDGOALS),clean)
-include $(dep)
endif

$(BUILD_DIR) :
	mkdir -p $(BUILD_DIR) 


$(BUILD_BIN_DIR)/libfsa.a: $(libfsa_obj)
	ar -r $@ $^

$(BUILD_BIN_DIR)/fsa_ol_assemble:$(fsa_ol_assemble_obj) $(BUILD_BIN_DIR)/libfsa.a
	g++ $^ -o $@ $(LDFLAGS)

$(BUILD_BIN_DIR)/fsa_ctg_bridge:$(fsa_ctg_bridge_obj) $(BUILD_BIN_DIR)/libfsa.a
	g++ $^ -o $@ $(LDFLAGS)

$(BUILD_BIN_DIR)/fsa_rd_haplotype:$(fsa_rd_haplotype_obj) $(BUILD_BIN_DIR)/libfsa.a
	g++ $^ -o $@ $(LDFLAGS)

$(BUILD_BIN_DIR)/fsa_ol_tools:$(fsa_ol_tools_obj) $(BUILD_BIN_DIR)/libfsa.a
	g++ $^ -o $@ $(LDFLAGS)

$(BUILD_BIN_DIR)/fsa_rd_tools:$(fsa_rd_tools_obj) $(BUILD_BIN_DIR)/libfsa.a
	g++ $^ -o $@ $(LDFLAGS)

$(BUILD_BIN_DIR)/fsa_misc_tools:$(fsa_misc_tools_obj) $(BUILD_BIN_DIR)/libfsa.a
	g++ $^ -o $@ $(LDFLAGS)

$(BUILD_BIN_DIR)/fsa_rd_correct: $(fsa_rd_correct_obj) $(BUILD_BIN_DIR)/libfsa.a
	g++ $^ -o $@ $(LDFLAGS)

$(BUILD_BIN_DIR)/fsa_ctg_correct: $(fsa_ctg_correct_obj) $(BUILD_BIN_DIR)/libfsa.a
	g++ $^ -o $@ $(LDFLAGS)

$(BUILD_BIN_DIR)/fsa_ol_purge: $(fsa_ol_purge_obj) $(BUILD_BIN_DIR)/libfsa.a
	g++ $^ -o $@ $(LDFLAGS)

	
$(BUILD_BIN_DIR)/fsa_ol_refine: $(fsa_ol_refine_obj) $(BUILD_BIN_DIR)/libfsa.a
	g++ $^ -o $@ $(LDFLAGS)

$(BUILD_BIN_DIR)/fsa_kmer_bin: $(fsa_kmer_bin_obj) $(BUILD_BIN_DIR)/libfsa.a
	g++ $^ -o $@ $(LDFLAGS)


%.o:%.c
	gcc $(CFLAGS) -c -MMD -MP -MF "$(@:%.o=%.d)" -MT "$(@:%.o=%.o)" -o "$@" "$<"

%.o:%.cpp
	g++ $(CXXFLAGS) -c -MMD -MP -MF "$(@:%.o=%.d)" -MT "$(@:%.o=%.o)" -o "$@" "$<"

$(BUILD_OBJ_DIR)/%.o:%.cpp
	g++ $(CXXFLAGS) -c -MMD -MP -MF "$(@:%.o=%.d)" -MT "$(@:%.o=%.o)" -o "$@" "$<"
	
#$(BUILD_OBJ_DIR)/edlib.o: ../thirdparty/edlib/edlib/src/edlib.cpp
#	g++ $(CXXFLAGS) -c -o "$@" "$<"

#$(BUILD_OBJ_DIR)/ksw2_extz.o: ../thirdparty/ksw2/ksw2_extz.c
#	gcc $(CFLAGS) -c -o "$@" "$<"


clean:
	rm $(obj) $(dep) $(bin) -f

