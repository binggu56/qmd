objects = mod.o cor.o init.o qm.o ran2.o gasdev.o eigen.o
FC = gfortran 
switch = -O3 -g 
EXEC = qm.x
LIB = /usr/lib64/atlas/liblapack.so.3.0
#LIB = -llapack  
$(EXEC): $(objects)
	$(FC) -o $(EXEC) $(LIB) $(switch) $(objects)

#cdat.o:mod.f
#	$(FC) -c mod.f
cdat.mod: mod.f90
	$(FC) -c mod.f90
%.o: %.f cdat.mod
	$(FC) -c  $<

%.o: %.f90 cdat.mod
	$(FC) -c  $<
clean:
	rm $(objects) *.mod
veryclean:
	rm *.o
