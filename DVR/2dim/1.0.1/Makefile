objects  = qm_dvr.o
FC = gfortran
switch = -O3 
EXEC = qm.x 
LIB = /usr/lib64/atlas/liblapack.so.3.0
#LIB = -llapack  
FLAGS = -w
$(EXEC): $(objects)
	$(FC) -o $(EXEC) $(LIB) $(FLAGS) $(switch) $(objects)
%.o: %.f
	$(FC) -c  $<
#%.O: %.f
#	$(FC) -c $<

clean:
	rm $(objects)
veryclean:
	rm *.o
