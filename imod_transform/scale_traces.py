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
        self.scale_factor = None
        self.scale_file = None
        self.scale_factors_dict = None # set when scale_file is loaded 
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
            '-f', '--scale-file', type=str, 
            help='scale file in the format where first column is object id and the second column is scaling for this object')
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

        if args.scale_file:
            self.scale_file = args.scale_file
            
        if self.scale_file and self.scale_factor: 
            print("Only one of scale factor and scale file must be set")
            return False
            
        if not self.scale_file and not self.scale_factor: 
            print("Input scale factor or scale file must be set")
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


def scale_contour_line(index, line, opts):
    # parse countour line into a 2D point with W coordinate
    data = line.split()
    assert len(data) == 3
    id = int(data[2])
    x = float(data[0])
    y = float(data[1])
    
    if opts.scale_factor:
        scale = opts.scale_factor
    else:
        assert opts.scale_factors_dict
        if index not in opts.scale_factors_dict:
            print("Error: amod file uses object with index " + str(index) + 
                  " but scale for this object was not supplied.")
            sys.exit(1)
        scale = opts.scale_factors_dict[index]
 
    # simply multiply the coordinates
    x = x * scale
    y = y * scale

    return str(x) + ' ' + str(y) + ' ' + str(id) + '\n'        


def load_scale_file(opts):
    res = {}
    with open(opts.scale_file, 'r') as f:
        for line in f:
            if not line.strip():
                continue
            index_scale = line.split()
            if len(index_scale) != 2:
                print("Error in " + opts.scale_file + " - could not process line '" + line + "'.")
                sys.exit(1)
            res[int(index_scale[0])] = float(index_scale[1])
    return res


def main():
    opts = Options()
    ok = opts.process_opts()
    if not ok:
        sys.exit(1)
    
    if opts.scale_file: 
        opts.scale_factors_dict = load_scale_file(opts)
        print("Loaded scales:")
        print(opts.scale_factors_dict)
    
    ok = transform_traces.check_amod_file(opts.input_amod_file)
    if not ok:
        sys.exit(1)
            
    transform_traces.process_amod_file(opts, None, None, scale_contour_line)
    
    print("Created file " + opts.output_amod_file + ".")
    
if __name__ == '__main__':
    main()            
        