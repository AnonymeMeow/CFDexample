del "%2\*.dat"
copy "%1\*.dat" "%2\"

echo See run detail at %2/outInfo.dat

%2\grid.exe
mpiexec -n 8 %2\scarf.exe

copy "%2\outInfo.dat" "%1\"
copy "%2\*_field.dat" "%1\"