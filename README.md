Int-DWTs
========

Repository to store the source code of the Int-DWTs library, being developed at LUPS, in UFPEL, RS, Brazil.

Maintainer: Vinícius Rodrigues dos Santos (`fireapache`)

Prerequisites
=============

You'll need `make` and `clang` and `cxsc` in order to compile the application.

`$ sudo apt-get install make`

`$ sudo apt-get install clang`

To get `cxsc` from terminal:

`$ wget http://www2.math.uni-wuppertal.de/~xsc/xsc/cxsc/cxsc-2-5-4.tar.gz`

`$ tar -zxvf cxsc-2-5-4.tar.gz`

`$ cd cxsc-2-5-4/`

`$ ./install_cxsc`

IMPORTANT: Use `/home/<you>/cxsc/` as installation path! (which is the default option)

How To Build
============

Pre-requisits:

`$ git clone https://www.github.com/fireapache/Int-DWTs/`

`$ cd Int-DWTs/`

`$ make tests`


How To Use
==========

And you can run the test envioriment by: `./test.exe`

* `-l` lists all tests.
* `-t` run a test
* `-ft` run a fundamental test

