mkdir $1
cd $1

awk '{print $1 "\t" $2 "\t"  $3 "\t" $4 "\t" $6 "\t"  $7 "\t" $8 "\t" $9>$1".tmp"}' ../$2.event
rm contig.tmp
head -n1 ../$2.event | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $6 "\t" $7 "\t" $8 "\t" $9}' > header
echo "done"

