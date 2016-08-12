objects = main.o
FC = gfortran
switch = -O3
EXEC = qm
EXEC2 = cont
LIB = /usr/lib64/atlas/liblapack.so.3.0
#LIB = -llapack  
$(EXEC): $(objects)
	$(FC) -o $(EXEC) $(LIB) $(switch) $(objects)

%.o: %.f
	$(FC) -c  $<
#%.O: %.f
#	$(FC) -c $<

clean:
	rm $(objects)
veryclean:
	rm *.o
