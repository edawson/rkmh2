CXX:=g++
CXXFLAGS:= -std=c++14 -O3 -fopenmp -mtune=native -march=native 
PROFILEFLAGS:= -std=c++14 -O2 -g -pg -fopenmp -mtune=native -march=native 
COMPATFLAGS:= -std=c++14 -O3 -fopenmp -mtune=ivybridge -march=ivybridge 
DEBUGFLAGS:= -std=c++14 -O0 -g -ggdb -fopenmp
LD_INC_FLAGS:= -Isparsepp/sparsepp -Ipliib -Ikseqpp/include/kseq++ -Imkmh
LD_LIB_FLAGS:= -lz

rkmh2: rkmh2.cpp mkmh/mkmh.hpp kseqpp/include/kseq++/config.hpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LD_INC_FLAGS) $(LD_LIB_FLAGS) 

profile:
	$(CXX) $(PROFILEFLAGS) -o $@ $^ $(LD_INC_FLAGS) $(LD_LIB_FLAGS)

debug-rkmh2: rkmh2.cpp mkmh/mkmh.hpp kseqpp/include/kseq++/config.hpp
	$(CXX) $(DEBUGFLAGS) -o $@ $^ $(LD_INC_FLAGS) $(LD_LIB_FLAGS) 

compat: rkmh2.cpp mkmh/mkmh.hpp kseqpp/include/kseq++/config.hpp
	$(CXX) $(COMPATFLAGS) -o rkmh2 $^ $(LD_INC_FLAGS) $(LD_LIB_FLAGS) 


kseqpp/include/kseq++/config.hpp:
	+cd kseqpp && cmake .

test: rkmh2
	./rkmh2 filter -r data/zika.refs.fa -f data/z1.fq -t 1

clean:
	$(RM) rkmh2
	$(RM) debug-rkmh2

.PHONY: test clean compat profile
