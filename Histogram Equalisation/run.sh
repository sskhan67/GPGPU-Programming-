rm  histrogram.pgm 


for ((j=0;j<30;j++))
do

./histo Lenna.pgm  $k >> histo.csv
done
