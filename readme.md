This repository presents a program to simulate reactive flow with three examples attached ([Cavity](./Cavity/), [Plate](./Plate/) and [Tube](./Tube/)).

This program is currently changing its structure from multi-process running (mpi) to multi-thread running.

# How to run

The [CMakeLists](./CMakeLists.txt) file defines three custom targets, each builds and runs the corresponding example.

To run the example you want to run, you just need to run

`cmake --build ./build --config Release --target ${target name}_build_and_run -j 10`

in your command line in the root directory of this repository. Remember to replace `${target name}` with the name of the example (`Cavity`, `Plate` or `Tube`) you want to run.