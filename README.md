# Postproc-Alibava 
Some utilities for the post-processing of data acquired with
the alibava systems [alibava DAQ](https://www.alibavasystems.com)

 * *fortythieves*: a utility to convert the raw binary data 
 format obtained from the Alibava DAQ into a ROOT trees

### Compilation
Create a build the directory of the source code:
```bash
mkdir build
cd build
cmake ..
make install
```
Per default it will create the executable *fortythieves* in the
```$HOME/.bin/``` directory, so you should have the environment
variable ```PATH``` pointing to that folder:
```bash
export PATH=$PATH:$HOME/.bin
```
You can change the install directory by using the cmake variable
```CMAKE_INSTALL_PREFIX```
```bash
cmake -DCMAKE_INSTALL_PREFIX=your_favorite_path ..
```
and then don't forget to export the ```PATH``` variable to include
that folder.

#### Dependencies
 * ROOT
 * CMAKE

### Usage
After succesful compilation and the exportation of the environment
variable, you are ready to use the executables of this package. 
Please take a look to the ```help``` option to use them:
```bash
usage: fortythieves [OPTIONS] alibava_data.raw

Convert a raw binary data from the ALIBAVA DAQ into a ROOT file

[OPTIONS]
 -o name of the ROOT output file [fortythieves.root]
 -p flag to store pedestal and noise header [false]
 -r run number [-1]
 -h show this help
```   

