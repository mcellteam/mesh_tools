"""
Copyright (C) 2020 by
The Salk Institute for Biological Studies

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

For the complete terms of the GNU General Public License, please see this URL:
http://www.gnu.org/licenses/gpl-2.0.html
"""

"""
example:
python scale_traces.py -a test/mito10_backchecked.amod -s 0.6457

creates file test/mito10_backchecked.amod.out
"""

import sys
import os
import argparse
import json
import numpy as np
from scipy import ndimage
from scipy import linalg

import transform_traces

class Options:
    def __init__(self):
        self.scale_factor = ''
        self.input_amod_file = ''
        self.output_amod_file = ''

    def __repr__(self):
        attrs = vars(self)
        return ", ".join("%s: %s" % item for item in attrs.items())
            
    @staticmethod
    def create_argparse():
        parser = argparse.ArgumentParser(description='Traces transformation tool')
        parser.add_argument(
            '-s', '--scale', type=str, 
            help='scale factor')
        parser.add_argument(
            '-a', '--amod', type=str, 
            help='amod file containing definition of traces done on top of the original images')
        parser.add_argument(
            '-o', '--output', type=str, 
            help='name of output .amod file')
        return parser

    def process_opts(self):
        
        parser = Options.create_argparse()
        args = parser.parse_args()
        
        if args.scale:
            self.scale_factor = float(args.scale)
        else:
            print("Input scale factor file must be set")
            return False
        
        if args.amod:
            self.input_amod_file = args.amod
        else:
            print("Input amod file must be set")
            return False
            
        if args.output:
            self.output_amod_file = args.output
        else:
            self.output_amod_file = self.input_amod_file + '.out'

        return True


def scale_contour_line(line, opts):
    # parse countour line into a 2D point with W coordinate
    data = line.split()
    assert len(data) == 3
    id = int(data[2])
    x = float(data[0])
    y = float(data[1])
 
    # simply multiply the coordinates
    x = x * opts.scale_factor
    y = y * opts.scale_factor

    return str(x) + ' ' + str(y) + ' ' + str(id) + '\n'        

def main():
    opts = Options()
    ok = opts.process_opts()
    if not ok:
        sys.exit(1)
    
    ok = transform_traces.check_amod_file(opts.input_amod_file)
    if not ok:
        sys.exit(1)
            
    transform_traces.process_amod_file(opts, None, None, scale_contour_line)
    
if __name__ == '__main__':
    main()            
        