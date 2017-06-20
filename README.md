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
 -p alibava raw-binary containing the pedestal run
 -c alibava raw-binary containing the calibration run
 -r run number [-1]
 -h show this help
```   
To obtain a root file with signal subtracted by pedestal and common noise, you need
to provide the ```-p``` option with the pedestal file obtained using a pedestal run
with the ALIBAVA-gui. You can provide a raw file obtained using a calibration run with 
the ```-c``` option as well. In that case, a new branch will contain the signal in 
number of electrons (pedestal and noise subtracted).
```bash
$ fortythieves -r <runNumber> -p <raw_pedestal_file> -c <calibration_file> -o outputfile.root <raw_beam_file> 
```
After the previous command is launched, the ```outputfile.root``` file will contain four TTrees:
* **runHeader**: the run header branches related extracted from the <raw_beam_file>. Note that the run header usually
contains data which is constant along the whole run.
* **Events**: the Event related branches (raw ADC counts per beetle, time per event, temperature, etc..)
* **postproc_runHeader**: this tree is only present whenever the ```-p``` or ```-c``` option is active, and it stores the pedestals and noise per channel calculated per chip, and the number of electrons per ADC count
   * ```pedestal_cmmd_beetle<chipnumber>```
   * ```noise_cmmd_beetle<chipnumber>```
   * ```electronADC_beetle<chipnumber>```
* **postproc_Events**: this tree is only present whenever the ```-p``` or ```-c``` option is active as well, and it stores the ADC counts per channel/event with the pedestals and noise subtracted, and the number of electrons per channel/event with the pedestals and noise subtracted (if ```-p``` option also)
   * ```postproc_data_beetle<chipnumber>```
   * ```postproc_cal_data_bettle<chipnumber>```

Note that in cases of ```-p``` and/or ```-c``` option present, an extra root file will be created 
with the result of the ROOT conversion performed at the pedestal and/or calibration raw files.
   
Use the ```AddFriend``` mechanism to relate and connect the original and the ```postproc``` versions of the Trees in order to share and use information between them.
```bash
# Add the Events tree as friend to postproc_Events (you can do it as they have the same
# number of events)
postproc_Events->AddFriend("Events");
# Now you can use the branches of Events (as eventTime) as if they belong to postproc_Events
postproc_Events->Draw("postproc_data_beetle1[13]","eventTime < 30 && eventTime > 3");
```
