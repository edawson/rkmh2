CXX:=g++
CXXFLAGS:= -std=c++14 -O3 -fopenmp -mtune=native -march=native 
LD_INC_FLAGS:= -Isparsepp/sparsepp -Ipliib -Ikseqpp/include/kseq++ -Imkmh
LD_LIB_FLAGS:= -lz

kramer: krmr.cpp mkmh/mkmh.hpp kseqpp/include/kseq++/kseq++.hpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LD_INC_FLAGS) $(LD_LIB_FLAGS) 

kseqpp/include/kseq++/kseq++.hpp:
	+cd kseqpp && cmake .

test: kramer
	./kramer filter -r data/zika.refs.fa -f data/z1.fq -t 1

clean:
	$(RM) kramer

.PHONY: test clean
