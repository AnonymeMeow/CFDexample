This branch is the original version of the program, whose performance will serve as the baseline for future optimizations.

Note: The original version of the program has a bug that skips some of the steps in the solving process, which makes the program 'faster'. But I still choose this version as the baseline, because the original program is from a published paper, by doing so I can better compare my results with the existing ones. (And I'm confident that I can make the program run faster ᕙ( •̀ ᗜ •́ )ᕗ)

This program can be built and run on Linux platforms, you just need to change the `MPI_DIR` in [CMakeLists file](./CMakeLists.txt) to the MPI installation directory on your computer and run

`cmake --build ./build --config Release --target ${target name}_build_and_run`

in your shell. The `${target name}` is the name of the target(`Cavity`, `Plate` or `Tube`).