mkdir $4
counter=1
until [ $counter -gt $1 ]
do
echo $counter

./SVM_multi.R -s $1 -p $counter -m $2 -u $3 -f $4 &

((counter++))
done
echo All done
