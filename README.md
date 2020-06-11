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

## To use EGGS:

EGGS can be plugged into current code with the following steps, suppose you have declared two Eigen sparse matrix:

```c++
Eigen::SparseMatrix<double> A;
Eigen::SparseMatrix<double> B;
// setup A and B
```

Then first, you will need to convert the structure using the implemented function in [utils.h](https://github.com/txstc55/EGGS/blob/master/include/igl/inline_expansion/utils.hpp) under include/igl/inline_expansion

```c++
Eigen::SparseMatrix<ie::NumericType> A_Numeric = ie::to_sparse_numeric<double, Eigen::ColMajor>(A, 0);
Eigen::SparseMatrix<ie::NumericType> B_Numeric = ie::to_sparse_numeric<double, Eigen::ColMajor>(B, 1);
```

The first input of this matrix will be the sparse matrix itself, the second one is the matrix id. The matrix id is the order of how you will pass in the datas when doing numeric computation. Note that, you will eventually group the input values in a vector, so the matrix ids should start from 0 and no gaps between ids.

At this step, you can perform symbolic executions, for example, you can do:

```c++
Eigen::SparseMatrix<ie::NumericType> Result_Numeric = A_Numeric * B_Numeric;
```

And the SparseMatrix Result_Numeric will have structual information of what is done to each entries in A and B to get the result. 

In order to assemble the final matrix, you should save the outerIndex and innerIndex from Result_Numeric since we only compute the values.

Now you will want to pass in the structure information to an executor by the following:

```c++
ie::NumericExecutor ex = ie::NumericExecutor(Result_Numeric, 0);
```

The 0 was a selection for method which was in a very old version, for now just leave it there.

The executor will generate code and compile it and link it against the application. 

After the compilation is done, you can use the executor to compute the numerical value. A typical way to do this when the value if A and B keeps changing is the following:

```c++
// The following are typically only done once
std::vector<vector<double>> DATAS;          // declare a vector for input data
DATAS.resize(2);                            // in this case, we have two input datas
DATAS[0].resize(A.nonZeros());              // reserve space for values of A
DATAS[1].resize(B.nonZeros());              // reserve space for values of B
std::vector<double> RESULTS;                // declare a vector for output
RESULTS.resize(Result_Numeric.nonZeros());  // reserve spcae for outputs

// The following are done everytime inputs are changed
DATAS[0].assign(A.valuePtr(), A.valuePtr() + A.nonZeros()); // assign the values of A
DATAS[1].assign(B.valuePtr(), B.valuePtr() + B.nonZeros()); // assign the values of B
ex.ExecuteMulti(DATAS, RESULTS);                            // execute the multi threaded version
```

After the execution, the numeric values will be in the RESULTS vector. You can pre-store the structure information of Result_Numeric (the outer index and inner index) to reassemble an Eigen Sparse Matrix. An implementation is done at [test_utils.hpp](https://github.com/txstc55/EGGS/blob/master/tutorial/803_Numeric/test_utils.hpp) callsed ConstructSparseMatrix.
