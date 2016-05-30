FC=gfortran
LFLAGS=-fopenmp -g
CFLAGS=-O3 -fopenmp -g

sort: sort.o mrgrnk.o
	$(FC) $(LFLAGS) -o $@ $^

%.o: %.f90
	$(FC) -c $(CFLAGS) -o $@ $^

clean:
	rm -f *.o sort
