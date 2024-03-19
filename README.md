### Beware: this is alpha code at best

This is an attempt to implement the core functionality of [The Virtual Brain](https://www.thevirtualbrain.org/tvb/zwei/brainsimulator-software) brain simulator [root](https://github.com/the-virtual-brain/tvb-root) project in C++.

### Dependencies:

* [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) version 3.4+
* Zlib
* Intel [Math Kernel Library](https://software.intel.com/content/www/us/en/develop/documentation/get-started-with-mkl-for-dpcpp/top.html) 

### Compilers:

So far, it has been tested with:

* VS C++ 2019 on Windows 10
* gcc 9.3 on Ubuntu 20.04
* Intel DPC++ 2021.2.0 on Windows 10

### Installation on Ubuntu Linux:

The dependencies can be found in the following system packages:

```libboost-all-dev libeigen3-dev```

First build the library.

```shell
$ cd src
$ mkdir build
$ cd build
$ cmake .. -DCMAKE_BUILD_TYPE=Release
$ make
```

Then proceed with the examples. 

```shell
$ cd examples
$ mkdir build
$ cd build
$ cmake .. -DCMAKE_BUILD_TYPE=Release
$ make
```


### Example of usage

Use the `test_simulation` to run from the command line.

```shell
examples/build/test_simulation -h
Options:
  -h [ --help ]                    Help screen
  --params arg                     Model parameters
  --noise arg                      Vector with noise sigmas for each state 
                                   variable
  --sc-matrix arg                  Structural connectivity matrix
  --model arg                      Whole brain model
  --length-matrix arg              Connection lengths matrix matrix
  --speed arg (=1000000)           Signal speed
  --time-start arg (=0)            Start of simulation (ms)
  --time-end arg (=10000)          End of simulation (ms)
  --dt arg (=0.100000001)          Integration step (ms)
  --params-file arg                NPZ file with simulation parameters
  --out-file-prefix arg (=out_sim) Output file prefix
```

Following line runs the Zerlaut model on a connectome of 66 nodes for 10s of simulated time. 

```shell
$ examples/build/test_simulation --model ZerlautAdaptationSecondOrder --sc-matrix Data_Raw/Human_66.npz --noise 0.0 0.0 0.0 0.0  0.0 0.0 0.0 0.1
```
