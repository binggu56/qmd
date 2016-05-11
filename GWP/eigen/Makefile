objects = qm.o 
FC = gfortran
switch = -O3
EXEC = qm
#LIB = /usr/lib64/atlas/liblapack.so.3.0
LIB = -llapack  
$(EXEC): $(objects)
	$(FC) -o $(EXEC) $(LIB) $(switch) $(objects)

#cdat.o:mod.f
#	$(FC) -c mod.f
#cdat.mod: mod.f
#	$(FC) -c mod.f
%.o: %.f 
	$(FC) -c  $<

clean:
	rm $(objects) *.mod
veryclean:
	rm *.o
