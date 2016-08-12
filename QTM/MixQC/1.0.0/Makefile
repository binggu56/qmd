objects = long_range.o reduce.o qm.o derivs.o input.o tstamp.o vinit.o fit.o gasdev.o ran1.o
objects2 = long_range.o reduce.o cont.o derivs.o input.o tstamp.o vinit.o qpot.o gasdev.o ran2.o
FC = ifort
switch = -O3
EXEC = qm
EXEC2 = cont
LIB = -mkl=sequential
#LIB = /usr/lib64/atlas/liblapack.so.3.0
#LIB = -llapack  
$(EXEC): $(objects)
	$(FC) -o $(EXEC) $(LIB) $(switch) $(objects)
$(EXEC2): $(objects2)
	$(FC) -o $(EXEC2) $(LIB) $(switch) $(objects2)


%.o: %.f
	$(FC) -c  $<
%.o: %.f90
	$(FC) -c $<

clean:
	rm $(objects)
veryclean:
	rm *.o
