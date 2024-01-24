MAKE := make

SRCS := DPMM.c DensityOutput.c PolyaUrnDPMM.c PosteriorDensityOutput.c PredictiveDensityOutput.c WalkerDPMM.c csv_helpers.c helpers.c js_output.c log.c plain_text.c read_data.c sort.c work.c
OBJS := $(patsubst %.c,obj/%.o,$(SRCS))
HEADERS := $(patsubst %.c,src/%.h,$(SRCS)) src/carbondate_internal.h

CXX := g++
CXXFLAGS := -std=c++11 -O2 -DOXCAL_RELEASE -pthread
INC := -IRmath/include

.PHONY: clean all Rmath

all: Rmath ex/carbondate

Rmath libRmath.a:
	cd Rmath; $(MAKE)

ex/carbondate: obj/main.o $(OBJS)
	mkdir -p ex/
	cp Oxcal.dat ex/
	cp curves/intcal* ex/
	cp curves/shcal* ex/
	$(CXX) -o ex/carbondate $(CXXFLAGS) obj/main.o $(OBJS) libRmath.a

obj/main.o: main.cpp $(HEADERS) inc/carbondate.h
	mkdir -p obj/
	$(CXX) -c $(CXXFLAGS) main.cpp $(INC) -o obj/main.o

obj/%.o: src/%.cpp $(HEADERS)
	mkdir -p obj/
	$(CXX) -c $(CXXFLAGS) $< $(INC) -o $@

clean:
	-rm -rf obj ex libRmath.a
