### Beware: this is aplha code at best

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

```libboost-all libeigen3-dev```

First build the library.

```shell
$ cd src
$ mkdir build
$ cd build
$ cmake ..
$ make
```

Then proceed with the examples. Make sure, that the `CPATH` environment variable points to directories with `Python.h` and `numpy` headers (e.g. `env/lib/python3.10/site-packages/numpy/core/include`).

```shell
$ cd examples
$ mkdir build
$ cd build
$ cmake ..
$ make
```
