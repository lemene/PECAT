
BUILD_BIN_DIR = ../build/bin
BUILD_OBJ_DIR = ../build/obj
BUILD_DIR += $(BUILD_BIN_DIR) $(BUILD_OBJ_DIR)

BUILD_DIR += $(BUILD_OBJ_DIR)/prog $(BUILD_OBJ_DIR)/../thirdparty/edlib/edlib/src $(BUILD_OBJ_DIR)/correct \
             $(BUILD_OBJ_DIR)/bridge $(BUILD_OBJ_DIR)/assemble $(BUILD_OBJ_DIR)/correct/align $(BUILD_OBJ_DIR)/kmer \
			 $(BUILD_OBJ_DIR)/utils $(BUILD_OBJ_DIR)/phase $(BUILD_OBJ_DIR)/assemble/simplify $(BUILD_OBJ_DIR)/assemble/graph \
			 $(BUILD_OBJ_DIR)/overlap $(BUILD_OBJ_DIR)/polish $(BUILD_OBJ_DIR)/align $(BUILD_OBJ_DIR)/base \
			 $(BUILD_OBJ_DIR)/mapping

simplify_src = assemble/simplify/simplifier.cpp ./assemble/simplify/spur.cpp ./assemble/simplify/bridge.cpp \
	./assemble/simplify/transitive.cpp ./assemble/simplify/cross.cpp ./assemble/simplify/low_quality.cpp \
	./assemble/simplify/unreliable.cpp ./assemble/simplify/loop.cpp ./assemble/simplify/bubble.cpp \
	./assemble/simplify/best.cpp ./assemble/simplify/semi_bubble.cpp ./assemble/simplify/repeat.cpp \
	./assemble/simplify/phase.cpp ./assemble/simplify/pbridge.cpp ./assemble/simplify/pspur.cpp \
	./assemble/simplify/ptransitive.cpp ./assemble/simplify/extend.cpp

align_src = align/match_info.cpp

graph_src = assemble/graph/cross_node.cpp

correct_src = ./correct/crr_options.cpp ./correct/crr_dataset.cpp ./correct/crr_options.cpp 

phase_src = ./phase/hic_read_infos.cpp

assemble_src = ./assemble/contig_generator.cpp

overlap_src = ./overlap/mapping.cpp

kmer_src = ./kmer/ranked_kmers.cpp ./kmer/kmer.cpp ./kmer/minimizer_counter.cpp ./kmer/kmer_counter.cpp \
           ./kmer/minimizer_graph.cpp ./kmer/minimizer_store.cpp ./kmer/mkseq_store.cpp


libfsa_src = ./utils/argument_parser.cpp ./utils/logger.cpp overlap.cpp read_store.cpp sequence.cpp\
             utility.cpp file_io.cpp seq_io.cpp overlap_store.cpp graph.cpp sequence_store.cpp\
			 assemble/overlap_assemble.cpp simple_align.cpp assemble/string_graph.cpp \
			 assemble/read_variants.cpp assemble/string_node.cpp assemble/string_edge.cpp \
			 assemble/asm_dataset.cpp assemble/asm_options.cpp assemble/stat_readinfo.cpp\
			 ../thirdparty/edlib/edlib/src/edlib.cpp correct/aligner.cpp ./correct/dw.cpp  ./correct/diff.cpp  \
			 ./correct/align/edlib_aligner.cpp ./correct/align/ksw2_aligner.cpp  \
			 ./utils/string_pool.cpp ./utils/pair_file.cpp ./utils/misc.cpp ./utils/project_file.cpp\
			 ./phase/phase_info.cpp ./phase/contig_phaser.cpp ./base/string_store.cpp \
			 $(phase_src) $(simplify_src) $(correct_src) $(graph_src) \
			 $(assemble_src) $(overlap_src) $(align_src) $(kmer_src)

libfsa_src += $(wildcard ./correct/align/*.cpp)

libfsa_obj:=$(patsubst %.cpp, $(BUILD_OBJ_DIR)/%.o, $(libfsa_src))

fsa_ol_assemble_src = ./prog/fsa_ol_assemble.cpp
fsa_ol_assemble_obj:=$(patsubst %.cpp, $(BUILD_OBJ_DIR)/%.o, $(fsa_ol_assemble_src))

fsa_ctg_bridge_src = ./prog/fsa_ctg_bridge.cpp bridge/contig_bridge.cpp bridge/contig_graph.cpp bridge/contig_link.cpp bridge/contig_link_store.cpp
fsa_ctg_bridge_obj:=$(patsubst %.cpp, $(BUILD_OBJ_DIR)/%.o, $(fsa_ctg_bridge_src))

fsa_rd_haplotype_src = ./prog/fsa_rd_haplotype.cpp ./phase/read_haplotype.cpp ./phase/local_phaser.cpp ./phase/utility.cpp  ./phase/phs_options.cpp ./phase/phs_dataset.cpp
fsa_rd_haplotype_obj:=$(patsubst %.cpp, $(BUILD_OBJ_DIR)/%.o, $(fsa_rd_haplotype_src))

fsa_ctg_call_src = ./prog/fsa_ctg_call.cpp ./phase/contig_call.cpp ./phase/local_phaser.cpp ./phase/utility.cpp  ./phase/phs_options.cpp ./phase/phs_dataset.cpp
fsa_ctg_call_obj:=$(patsubst %.cpp, $(BUILD_OBJ_DIR)/%.o, $(fsa_ctg_call_src))

fsa_rd_consistent_src = ./prog/fsa_rd_consistent.cpp ./phase/read_consistent.cpp ./phase/local_phaser.cpp ./phase/utility.cpp  ./phase/phs_options.cpp ./phase/phs_dataset.cpp
fsa_rd_consistent_obj:=$(patsubst %.cpp, $(BUILD_OBJ_DIR)/%.o, $(fsa_rd_consistent_src))

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

fsa_ctg_polish_src = ./prog/fsa_ctg_polish.cpp ./correct/alignment_graph.cpp\
                     ./polish/contig_polish.cpp ./polish/pol_dataset.cpp ./polish/pol_options.cpp\
					 ./polish/contig_analyzer.cpp ./polish/coverage_info.cpp ./polish/window_slider.cpp \
					 ./polish/contig_fragment.cpp ./polish/contig_graph.cpp ./polish/multi_coverage.cpp
					 
fsa_ctg_polish_obj := $(patsubst %.cpp, $(BUILD_OBJ_DIR)/%.o, $(fsa_ctg_polish_src))

fsa_ol_purge_src = ./prog/fsa_ol_purge.cpp ./overlap_purge.cpp
fsa_ol_purge_obj:=$(patsubst %.cpp, $(BUILD_OBJ_DIR)/%.o, $(fsa_ol_purge_src))

fsa_ol_refine_src = ./prog/fsa_ol_refine.cpp ./overlap_refine.cpp
fsa_ol_refine_obj:=$(patsubst %.cpp, $(BUILD_OBJ_DIR)/%.o, $(fsa_ol_refine_src))

fsa_kmer_tools_src = ./prog/fsa_kmer_tools.cpp ./kmer/mprog_kmer_tools.cpp
fsa_kmer_tools_obj:=$(patsubst %.cpp, $(BUILD_OBJ_DIR)/%.o, $(fsa_kmer_tools_src))

fsa_sam_tools_src = ./prog/fsa_sam_tools.cpp ./sam_tools.cpp
fsa_sam_tools_obj:=$(patsubst %.cpp, $(BUILD_OBJ_DIR)/%.o, $(fsa_sam_tools_src))

fsa_rd_call_src = ./prog/fsa_rd_call.cpp ./phase/read_call.cpp
fsa_rd_call_obj:=$(patsubst %.cpp, $(BUILD_OBJ_DIR)/%.o, $(fsa_rd_call_src))

fsa_map_tools_src = ./prog/fsa_map_tools.cpp ./mapping/mprog_map_tools.cpp \
					./mapping/prog_coverage_profile.cpp ./polish/multi_coverage.cpp
fsa_map_tools_obj:=$(patsubst %.cpp, $(BUILD_OBJ_DIR)/%.o, $(fsa_map_tools_src))

ifeq ($(BUILDDEBUG), 1)
	CFLAGS =  -Wall -O0  -D_FILE_OFFSET_BITS=64 -g -I.  -I ../thirdparty/edlib/edlib/include -I ../thirdparty/ksw2  -I ../thirdparty/htslib
	CXXFLAGS = -std=c++14 -Wall -O0  -D_FILE_OFFSET_BITS=64 -g -I.  -I ../thirdparty/edlib/edlib/include -I ../thirdparty/ksw2  -I ../thirdparty/htslib
else
	CFLAGS =  -Wall -O3  -D_FILE_OFFSET_BITS=64 -I. -I ../thirdparty/edlib/edlib/include -I ../thirdparty/ksw2 -I ../thirdparty/htslib
	CXXFLAGS = -std=c++14 -Wall -O3 -D_FILE_OFFSET_BITS=64 -I.  -I ../thirdparty/edlib/edlib/include -I ../thirdparty/ksw2  -I ../thirdparty/htslib
endif

LDFLAGS = -pthread -lz $(BUILD_BIN_DIR)/libksw2.a ../thirdparty/htslib/libhts.a  


obj:=$(libfsa_obj) $(fsa_ol_assemble_obj) \
	$(fsa_ctg_bridge_obj) $(fsa_rd_correct_obj) \
	$(fsa_rd_haplotype_obj) $(fsa_rd_consistent_obj) $(fsa_ctg_call_obj) \
	$(fsa_ol_tools_obj) $(fsa_rd_tools_obj) $(fsa_misc_tools_obj) \
	$(fsa_ol_purge_obj) $(fsa_ol_refine_obj) $(fsa_kmer_bin_obj) \
	$(fsa_ctg_polish_obj) $(fsa_sam_tools_obj) $(fsa_kmer_tools_obj) \
	$(fsa_map_tools_obj) $(fsa_rd_call_obj)


dep:=$(patsubst %.o,%.d, $(obj))

bin = libfsa.a fsa_ol_assemble fsa_ctg_bridge \
      fsa_rd_haplotype fsa_ol_tools fsa_rd_tools fsa_sam_tools\
	  fsa_misc_tools fsa_rd_correct fsa_ol_purge fsa_ol_refine \
	  fsa_ctg_call fsa_rd_consistent fsa_ctg_polish \
	  fsa_kmer_tools fsa_map_tools #fsa_rd_call

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

$(BUILD_BIN_DIR)/fsa_ctg_call:$(fsa_ctg_call_obj) $(BUILD_BIN_DIR)/libfsa.a
	g++ $^ -o $@ $(LDFLAGS)

$(BUILD_BIN_DIR)/fsa_rd_consistent:$(fsa_rd_consistent_obj) $(BUILD_BIN_DIR)/libfsa.a
	g++ $^ -o $@ $(LDFLAGS)

$(BUILD_BIN_DIR)/fsa_ol_tools:$(fsa_ol_tools_obj) $(BUILD_BIN_DIR)/libfsa.a
	g++ $^ -o $@ $(LDFLAGS)

$(BUILD_BIN_DIR)/fsa_rd_tools:$(fsa_rd_tools_obj) $(BUILD_BIN_DIR)/libfsa.a
	g++ $^ -o $@ $(LDFLAGS)

$(BUILD_BIN_DIR)/fsa_misc_tools:$(fsa_misc_tools_obj) $(BUILD_BIN_DIR)/libfsa.a
	g++ $^ -o $@ $(LDFLAGS)

$(BUILD_BIN_DIR)/fsa_rd_correct: $(fsa_rd_correct_obj) $(BUILD_BIN_DIR)/libfsa.a
	g++ $^ -o $@ $(LDFLAGS)

$(BUILD_BIN_DIR)/fsa_ctg_polish: $(fsa_ctg_polish_obj) $(BUILD_BIN_DIR)/libfsa.a
	g++ $^ -o $@ $(LDFLAGS)

$(BUILD_BIN_DIR)/fsa_ol_purge: $(fsa_ol_purge_obj) $(BUILD_BIN_DIR)/libfsa.a
	g++ $^ -o $@ $(LDFLAGS)

	
$(BUILD_BIN_DIR)/fsa_ol_refine: $(fsa_ol_refine_obj) $(BUILD_BIN_DIR)/libfsa.a
	g++ $^ -o $@ $(LDFLAGS)

$(BUILD_BIN_DIR)/fsa_kmer_bin: $(fsa_kmer_bin_obj) $(BUILD_BIN_DIR)/libfsa.a
	g++ $^ -o $@ $(LDFLAGS)

$(BUILD_BIN_DIR)/fsa_kmer_tools: $(fsa_kmer_tools_obj) $(BUILD_BIN_DIR)/libfsa.a
	g++ $^ -o $@ $(LDFLAGS)

$(BUILD_BIN_DIR)/fsa_sam_tools: $(fsa_sam_tools_obj) $(BUILD_BIN_DIR)/libfsa.a
	g++ $^ -o $@ $(LDFLAGS)


# $(BUILD_BIN_DIR)/fsa_rd_call: $(fsa_rd_call_obj) $(BUILD_BIN_DIR)/libfsa.a
# 	echo $(fsa_rd_call_obj) $(BUILD_BIN_DIR)/libfsa.a
# 	g++ $^ -o $@ $(LDFLAGS)


$(BUILD_BIN_DIR)/fsa_map_tools: $(fsa_map_tools_obj) $(BUILD_BIN_DIR)/libfsa.a
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

