objects = init.o qm.o ran2.o  gasdev.o eigen.o
FC =  gfortran 
switch = -O3 -g 
EXEC = qm.x
#LIB = /usr/lib64/atlas/liblapack.so.3.0
LIB = -llapack  
$(EXEC): $(objects)
	$(FC) -o $(EXEC) $(LIB) $(switch) $(objects)

#cdat.o:mod.f
#	$(FC) -c mod.f

%.o: %.f 
	$(FC) -c  $<

%.o: %.f90
	$(FC) -c  $<
clean:
	rm $(objects) *.mod
veryclean:
	rm *.o
