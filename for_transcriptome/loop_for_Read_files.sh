counter=1
until [ $counter -gt $1 ]
do
echo $counter

tmp_all_fil.R -s $1 -p $counter &

((counter++))
done
echo All done
