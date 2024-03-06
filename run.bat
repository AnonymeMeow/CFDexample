cd %3

del "*.dat"
copy "%1\*.dat" "."

echo See run detail at %3/outInfo.dat

%2\grid.exe
mpiexec -n 8 %3\scarf.exe

copy "outInfo.dat" "%1\"
copy "*_field.dat" "%1\"