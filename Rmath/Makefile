SRCS := $(wildcard *.c)
OBJS := $(patsubst %.c,../obj/%.o,$(SRCS))

CC := gcc
CCFLAGS := -O2 -DMATHLIB_STANDALONE -pthread
INC := -Iinclude -Iinclude/R_ext

.PHONY: clean all

all: ../libRmath.a

.NOTPARALLEL ../libRmath.a: $(OBJS)
	ar rc $@ $(OBJS)

../obj/%.o: %.c include/Rmath.h
	mkdir -p ../obj/
	$(CC) -c $< $(INC) $(CCFLAGS) -o $@

clean:
	-rm -rf $(OBJS) ../libRmath.a
