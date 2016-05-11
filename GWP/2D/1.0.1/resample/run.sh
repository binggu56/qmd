
#if [ -e $1 ]
#  then
#	rm -r $1
#fi
#mkdir $1




for i in `seq 1 $1`;
do
	cp record.dat temp.dat 
	./qm.x | tee log
done    


