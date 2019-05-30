ls -1 */*.tmp | awk -F '/' '{print $2}' | sort | uniq > names

mkdir combined
cd combined
for i in `cat ../names` ; do
    cat ../header ../*/$i > $i
done
echo done
