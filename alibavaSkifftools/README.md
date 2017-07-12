# alibavaSkifftools
Configuration for the Alibava Marlin processors defined at
the [EUTelescope](https://github.com/duartej/eutelescope) package. 
Also it deals with the static data produced during the IFCA test-beam
campaign realized at CERN SPS during May-2017.
 * Scripts
   * *open_sesame*: multitool suitable to prepare and create Marlin jobs
   related with the Alibava and the Eutelescope data. The tool also extracts
   info related with the static data created at the May-2017 Test-beam.

 * Modules
   * *SPS2017TB_metadata*: May-2017 test-beam data-related
   * *steering_processing*: prepare and create alibava and Eutelescope Marlin
   processor jobs

author: Jordi Duarte-Campderros (June.2017)

### Compilation
See [README](https://github.com/duartej/postproc-alibava)

### Usage
#### open_sesame multitool
The tool has two main subcommands: 
 * *list_files*: list all the available ALIBAVA RAW data files and
 associate each **beam** file with their **pedestal** and **calibration** 
 file
 * *steering*  : build the needed steering files to run a given step
 of the marlin framework reconstruction for the ALIBAVA or TELESCOPE data.


#### The name of the package?
This package is intended to configure and prepare the whole set
of Alibava Marlin processors through the use of xml steering 
files: A pedantic reference to the skiff and the tools used by 
Santiago (the Hemingway's hero from *The old man and the sea*) 
to catch his *marlin*. 
