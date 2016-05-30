FC=gfortran
LFLAGS=-fopenmp -g
CFLAGS=-O3 -fopenmp -g

all: test

test: sort.o mrgrnk.o test.o
	$(FC) $(LFLAGS) -o $@ $^

%.o: %.f90
	$(FC) -c $(CFLAGS) -o $@ $^

clean:
	rm -f *.o sort
