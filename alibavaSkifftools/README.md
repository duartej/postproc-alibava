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
 * *list_files*: available ALIBAVA RAW data files from CERN-SPS TB May-2017
 campaign.
   * List all the available data present in the folder <input_folder>
   ```bash
   $ open_sesame list_files /eos/user/d/duarte/alibavas_data
   ```
   * Get the **pedestal** and **calibration** file associated to a given **beam** file (_-b option_)
   ```bash
   $ open_sesame list_files /eos/user/d/duarte/alibavas_data -b 393_2017-05-21_17-48_gerva_MB2_N1-3_-50V_-76d0uA_-25C_lat132_beam.dat
   ```
 * *steering*  : build the needed steering files to run a given step
 of the marlin framework reconstruction for the ALIBAVA or TELESCOPE data.
   * Get the list of available reconstruction steps
   ```bash
   $ open_sesames steering_files -p
   ```
   * Get the full reconstruction chain for an alibava raw data
   ```bash
   $ ALIBAVA_DATA=/eos/user/d/duarte/alibavas_data
   $ open_sesame steering alibava_full_reco --alibava-input-filename ${ALIBAVA_DATA}/N1-3_0_b1/393_2017-05-21_17-48_gerva_MB2_N1-3_-50V_-76d0uA_-25C_lat132_beam.dat --pedestal-input-filename ${ALIBAVA_DATA}/N1-3_0_b1/2017-05-21_15-23_gerva_MB2_N1-3_-50V_-80d0uA_-25C_lat132_ped.dat
   # The script also creates a bash file to run the whole chain
   $ ./alibava_full_reconstruction.sh
   ```


##### The name of the package?
This package is intended to configure and prepare the whole set
of Alibava Marlin processors through the use of xml steering 
files: A pedantic reference to the skiff and the tools used by 
Santiago (the Hemingway's hero from *The old man and the sea*) 
to catch his *marlin*. 
