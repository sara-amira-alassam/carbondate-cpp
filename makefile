MAKE := make
MKFLAGS =

CC := g++
CCFLAGS := -std=c++11 -O2

OBJ = DPMM.o DensityOutput.o PolyaUrnDPMM.o PosteriorDensityOutput.o PredictiveDensityOutput.o WalkerDPMM.o csv_helpers.o helpers.o log.o plain_text.o read_data.o sort.o work.o
OBJ := $(patsubst %,obj/%,$(OBJ))

.PHONY: clean all Rmath carbondate

all: Rmath carbondate

Rmath Rmath.a:
	cd Rmath; $(MAKE) $(MKFLAGS)

carbondate: obj/main.o $(OBJ)
	g++ -o carbondate $(CCFLAGS) obj/main.o $(OBJ) Rmath.a

obj/main.o:
	g++ -c $(CCFLAGS) main.cpp -IRmath/include -o obj/main.o

obj/%.o: src/%.cpp
	g++ -c $(CCFLAGS) $< -Iinclude -IRmath/include -o $@

clean:
	-rm -rf obj ex Rmath.a Rmath/*.o
