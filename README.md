# Libigl-Inline-Expansion

This project uses inline-expansion, which generates code for a given structure of (sparse) matrix and operation, for some of the functions implemented in [libigl](https://libigl.github.io/). The examples we pick normally fufil the following requirement:

1. Once the structure and how to produce the result is set, the structure of the final matrix remain the same, only the values will differ based on input.
2. The number of operations is finite, it will not be a random operation on each entry, but some pattern is formed.

## Requirement

We use [Intel MKL](https://software.intel.com/en-us/mkl) not only for comparison for some of the examples, but also use the pardiso solver implemented in MKL for solving a sparse linear equation. Hence, downloading and installing is necessary. 

After you installed MKL, please set systemv ariable <em>MKLROOT</em> to where you installed MKL, normally the line of code works:

```bash
export MKLROOT=/opt/intel/mkl/
```

The project is built using cmake, so please make sure cmake is installed on your machine.

We also used clang as our compiler during the testing due to part of the code generated cannot be compiled using gcc, specifically, it is the direct initialization of a __mm128d variable like this:

```c++
__mm128d m1 = {a, b};
```

So please change the compiler to clang by doing the following:

```bash
export CC=/usr/bin/clang
export CXX=/usr/bin/clang++
```

After this, cmake will recognize clang as its default compiler.

## Build

After you have downloaded this project, please run the following code to build the project:

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j
```

## Run an example

To run an example, please do:

```bash
cd build && ./tutorial/xxxx_bin
```

The GUI is removed from all the example unlike the original ones in libigl since we only need to profile the speed.

The output vary from test to test, if you run SLIM or Cotangent Flow example, then the final object will be dumped, you will find them in build folder.

If you run optical flow example, the final flow in U and V direction will be dumped as two text files which you can later read.

If you run Cotangent Matrix example, only the speed result and the error will be shown.

The numeric tests are a set of tests for profiling our performance of some operations on random matrix comparing to MKL, the result will be written to a text file in build folder.

## Generate speed result for all tests

A bash script is provided for doing all of the tests, another bash script is provided to produce all the graphs. To reproduce the result in the paper, please do:

```bash
./all_tests.sh
./all_plot.sh
```

You can also just pick the test to run by just going into the script and pick the ones you want to do.

The numeric tests will take longer and is recommended to be commented out if you just want to see the results for real world applications.
