MAKE := make
MKFLAGS =

CC := g++
CCFLAGS := -std=c++11 -O2 -DOXCAL_RELEASE -pthread

OBJ = DPMM.o DensityOutput.o PolyaUrnDPMM.o PosteriorDensityOutput.o PredictiveDensityOutput.o WalkerDPMM.o csv_helpers.o helpers.o log.o plain_text.o read_data.o sort.o work.o
OBJ := $(patsubst %,obj/%,$(OBJ))

.PHONY: clean all Rmath

all: Rmath ex/carbondate

Rmath libRmath.a:
	cd Rmath; $(MAKE) $(MKFLAGS)

ex/carbondate: obj/main.o $(OBJ)
	mkdir -p ex/
	cp Oxcal.dat ex/
	cp curves/intcal* ex/
	$(CC) -o ex/carbondate $(CCFLAGS) obj/main.o $(OBJ) libRmath.a

obj/main.o:
	mkdir -p obj/
	$(CC) -c $(CCFLAGS) main.cpp -IRmath/include -o obj/main.o

obj/%.o: src/%.cpp src/%.h
	mkdir -p obj/
	$(CC) -c $(CCFLAGS) $< -Iinclude -IRmath/include -o $@

clean:
	-rm -rf obj ex libRmath.a
