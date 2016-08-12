if [ -e $1 ]
  then
	rm -r $1
fi
mkdir $1
cp HigherOrder_1D.py $1
cp plt.py $1 
cd $1
~/anaconda3/bin/python3 HigherOrder_1D.py 

