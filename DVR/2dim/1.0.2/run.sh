if [ -e $1 ]
  then
    rm -r $1
fi

mkdir $1
cp qm.x INdvr qm_dvr.f $1
#cp *.py $1
#cp pople $1
cd $1
#qsub pople
./qm.x 
