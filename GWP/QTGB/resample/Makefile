objects = mod.o qm.o ham.o ran2.o gasdev.o eigen.o
FC = gfortran 
switch = -O3 -g 
EXEC = qm.x
LIB = /usr/lib64/atlas/liblapack.so.3.0
#LIB = -llapack  
$(EXEC): $(objects)
	$(FC) -o $(EXEC) $(LIB) $(switch) $(objects)

cdat.mod : mod.f90
	$(FC) -c mod.f90

#cdat.mod: mod.f90
#	$(FC) -c mod.f90

%.o: %.f 
	$(FC) -c  $<

%.o: %.f90
	$(FC) -c  $<

clean:
	rm $(objects) *.mod
veryclean:
	rm *.o
