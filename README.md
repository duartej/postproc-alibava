# Postproc-Alibava   [![Build Status](https://travis-ci.org/duartej/postproc-alibava.svg?branch=master)](https://travis-ci.org/duartej/postproc-alibava)
Some utilities for the post-processing of data acquired with
the alibava systems [alibava DAQ](https://www.alibavasystems.com)

 * *fortythieves*: a utility to convert the raw binary data 
 format obtained from the Alibava DAQ into a ROOT trees

 author: Jordi Duarte-Campderros (June.2017)

### Compilation
Create a *build* directory. Configure and compile the code using
cmake inside the build directory:
```bash
$ mkdir build
$ cd build
$ cmake ..
$ make install
```
Per default it will create the executable *fortythieves* in the
```$HOME/.local/bin/``` directory, so you should have the environment
variable ```PATH``` pointing to that folder:
```bash
$ export PATH=$PATH:$HOME/.local/bin
```
As well, it will create some shared libraries under the 
```$HOME/.local/lib```, so you need to update also the environment
variable ```LD_LIBRARY-PATH```:
```bash
$ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/.local/lib
```
You can change the install directory by using the cmake variable
```CMAKE_INSTALL_PREFIX```
```bash
$ cmake -DCMAKE_INSTALL_PREFIX=your_favorite_path ..
```
and then don't forget to export the ```PATH``` and the ```LD_LIBRARY_PATH``` to
include both folders.

#### Dependencies
 * ROOT >= 6.0 (Note that some problems has been spotted when using 5.34)
 * CMAKE >= 2.8

### Usage
After succesful compilation and the exportation of the environment
variable, you are ready to use the executables of this package. 
Please take a look to the ```help``` option to use them:
```bash
$ fortythieves -h
usage: fortythieves [OPTIONS] alibava_data.raw

Convert a raw binary data from the ALIBAVA DAQ into a ROOT file

[OPTIONS]
 -o name of the ROOT output file [fortythieves.root]
 -p flag to store pedestal and noise header [false]
 -r run number [-1]
 -h show this help
```   

