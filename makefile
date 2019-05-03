VER=0.1
DIST=qubic$(VER)
PROGS=qubic
SRCS=struct.c read_array.c make_graph.c get_options.c fib.c write_block.c main.c expand.c
OBJS=$(SRCS:.c=.o) 
CC=g++

CUDA_INCLUDEPATH=/usr/local/cuda/include

NVCC=nvcc
NVCC_OPTS=-arch=sm_30
GCC_OPTS=-std=c++11 -g -O0 -Wall
CUDA_HOME=/usr/local/cuda/include
CUDA_LD_FLAGS=-L$(CUDA_HOME) -lcuda -lcudart 

LDFLAGS=-static  -lm 
CFLAGS=-O0 -Wall -ansi -I.  -DVER=$(VER)

all: $(PROGS)

cluster.o:cluster.cu cluster.h utils.h cluster_parallel.cu
	$(NVCC) -c cluster.cu $(NVCC_OPTS)

${PROGS}: $(OBJS) cluster.o
	$(CC) $(OBJS) -o $@ $(CUDA_LD_FLAGS) #$(LDFLAGS)

.o:
	$(CC) $(GCC_OPTS) $< -o $@ 

clean:
	rm -f $(PROGS)
	rm -f *.o
	rm -f *.rules
	rm -f *.chars
	rm -f *.blocks
	rm -f *.expansion

dist:
	$(MAKE) clean
	cd .. && tar czvf $(DIST).tar.gz $(DIST)/

Ecoli.chars:
	./${PROGS} -i Ecoli

CRP.blocks:
	./${PROGS} -i CRP 

test: Ecoli.chars CRP.blocks
	./${PROGS} -i Ecoli.chars -b CRP.blocks -s
