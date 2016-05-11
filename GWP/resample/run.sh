
#if [ -e $1 ]
#  then
#	rm -r $1
#fi
#mkdir $1
cd .. 
make 
./qm.x 
cd resample
rm *.dat
cp ../temp.dat . 
cp ../cor.dat .

for i in `seq 1 $1`;
do
#	cp record.dat temp.dat 
	./qm.x | tee log
done    


