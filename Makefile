CXX:=g++
CXXFLAGS:= -std=c++14 -O3 -fopenmp -mtune=native -march=native 
LD_INC_FLAGS:= -Isparsepp/sparsepp -Ipliib -Ikseqpp/include

krmr: krmr.cpp 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LD_INC_FLAGS) $(LD_LIB_FLAGS) 
