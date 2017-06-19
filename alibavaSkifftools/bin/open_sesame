#!/usr/bin/env python
"""
Extract the list of alibava raw data files from a parent folder.
The folder should contain subdirectories with the name
of the involved sensors, and the file names follows a pre-defined
naming convention.
"""
__author__ = "Jordi Duarte-Campderros"
__credits__ = ["Jordi Duarte-Campderros"]
__version__ = "v0.1"
__maintainer__ = "Jordi Duarte-Campderros"
__email__ = "jorge.duarte.campderros@cern.ch"
__status__ = "Development"

def main(parent_folder):
    """
    """
    import alibavaSkifftools.SPS2017TB_metadata import *

if __name__ == '__main__':
    from argparse import ArgumentParser
    from alibavaSkifftools.SPS2017TB_metadata import eospath

    usage  = "Extract the list of alibava raw data files from a parent folder"
    usage += "The folder should contain subdirectories with the name of the"
    usage += " involved sensors; and the file names follow a pre-defined naming"
    usage += "convention"

    parser = ArgumentParser(prog='open_sesame',description=usage)

    parser.add_argument('parent_folder',help="The parent folder to start"\
            " to search down [Default: see SPS2017TB_metadata.eospath]")
    parser.set_defaults(parent_folder=eospath)

    main(parser.parent_folder)
