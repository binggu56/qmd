if [ -e $1 ]
  then
    rm -r $1
fi

mkdir $1
cp qm.x IN qm.f $1
cp en.py $1
cp derivs.f90 $1
#cp pople $1
cd $1
#qsub pople
./qm.x & 
