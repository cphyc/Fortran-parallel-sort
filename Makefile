FC=gfortran
LFLAGS=-O3 -fopenmp
CFLAGS=-O3 -fopenmp

all: test

test: m_mrgrnk.o mod_sort.o test.o
	$(FC) $(LFLAGS) -o $@ $^

%.o: %.f90
	$(FC) -c $(CFLAGS) -o $@ $^

clean:
	rm -f *.o *.mod test