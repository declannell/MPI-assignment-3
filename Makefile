CC = mpicc

CFLAGS = -g -Wall #-D_CC_OVERLAP

LDFLAGS = -lm

POISSOBJS = decomp1d.o jacobi.o
#POISSOBJS = decomp1d.o jacobi.o gfunc.o

EXECS = rma_blocking_2d_decomp rma_pscw_2d_decomp

all: $(EXECS)

rma_blocking_2d_decomp: rma_blocking_2d_decomp.o  $(POISSOBJS) 
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

rma_pscw_2d_decomp: rma_pscw_2d_decomp.o  $(POISSOBJS) 
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)


tags:
	etags *.c *.h

.PHONY: clean tags tests

clean:
	$(RM) *.o $(EXECS) $(TESTS) TAGS tags
