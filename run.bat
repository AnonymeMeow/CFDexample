cd %2

del "*.dat"
copy "%1\*.dat" "."

echo See run detail at %2/outInfo.dat

.\grid.exe
mpiexec -n 8 %2\scarf.exe

copy "outInfo.dat" "%1\"
copy "*_field.dat" "%1\"