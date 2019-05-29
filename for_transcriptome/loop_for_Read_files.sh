mkdir $3
counter=1
until [ $counter -gt $1 ]
do
echo $counter

./Read_events_multi.R -s $1 -p $counter -i $2 -o $3 &

((counter++))
done
echo All done
