cd $2

rm *.dat
cp $1/*.dat .

echo See run detail at $2/outInfo.dat

./scarf

cp outInfo.dat $1/
cp *_field.dat $1/