objects = qm.o ran2.o mod.o qpot.o gasdev.o
FC = gfortran
switch = -O3
EXEC = qm
#LIB = /usr/lib64/atlas/liblapack.so.3.0
LIB = -llapack  
$(EXEC): $(objects)
	$(FC) -o $(EXEC) $(LIB) $(switch) $(objects)

#cdat.o:mod.f
#	$(FC) -c mod.f
cdat.mod: mod.f
	$(FC) -c mod.f
%.o: %.f cdat.mod
	$(FC) -c  $<

clean:
	rm $(objects) *.mod
veryclean:
	rm *.o
