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

# https://towardsdatascience.com/image-geometric-transformation-in-numpy-and-opencv-936f5cd1d315
# https://docs.opencv.org/2.4/modules/imgproc/doc/geometric_transformations.html?highlight=warpaffine#void%20warpAffine(InputArray%20src,%20OutputArray%20dst,%20InputArray%20M,%20Size%20dsize,%20int%20flags,%20int%20borderMode,%20const%20Scalar&%20borderValue)
# 2:25 of https://robotacademy.net.au/lesson/describing-rotation-and-translation-in-2d/#:~:text=The%20homogeneous%20transformation%20matrix%20T,in%20this%20single%203x3%20matrix.

"""
display traces:

convert -compress none *.jpg picts.tif

tif2mrc -p 16.0 *.tif proj.mrc
alterheader proj.mrc -org "3872 3408 1824"
"""

import sys
import os
import argparse
import json
import numpy as np
from scipy import ndimage
from scipy import linalg


class NumpyArrayEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)
    
    
class Options:
    def __init__(self):
        self.xg_file = ''
        self.sections_file = ''
        self.output_json_file = ''

    def __repr__(self):
        attrs = vars(self)
        return ", ".join("%s: %s" % item for item in attrs.items())
            
    @staticmethod
    def create_argparse():
        parser = argparse.ArgumentParser(description='Traces transformation tool')
        parser.add_argument(
            '-x', '--xg', type=str, 
            help='xg file')
        parser.add_argument(
            '-s', '--sections', type=str, 
            help='sections file'
            )
        parser.add_argument(
            '-o', '--output', type=str, 
            help='name of output .json file')
        return parser

    def process_opts(self):
        
        parser = Options.create_argparse()
        args = parser.parse_args()
        
        if args.xg:
            self.xg_file = args.xg
        else:
            print("Input xg file must be set")
            return False
        
        if args.sections:
            self.sections_file = args.sections
        else:
            print("Input sections file must be set")
            return False  

        if args.output:
            self.output_json_file = args.output
        else:
            self.output_json_file = self.xg_file + '.json'

        return True
    
def load_xg_file(file_name):
    # returns a list 2x3 numpy arrays
    res = []
    with open(file_name, 'r') as f:
        for line in f:
            items = line.split()
            assert len(items) == 6
            mat = np.empty([2, 3])
            mat[0][0] = float(items[0])
            mat[1][0] = float(items[1])
            mat[0][1] = float(items[2])
            mat[1][1] = float(items[3])
            mat[0][2] = float(items[4])
            mat[1][2] = float(items[5])
            res.append(mat)
    
    return res

def load_sections_file(file_name):
    # returns a list of pairs
    res = []
    with open(file_name, 'r') as f:
        for line in f:
            items = line.split()
            assert len(items) == 2
            indices = items[1].split('-')
            assert len(indices) == 2
            res.append( (int(indices[0]), int(indices[1])) )
    
    return res


def generate_transforms_file(opts, transform_matrices, sections):
    # prepare json dict
    transforms = []
    for k in range(len(sections)):
        s = sections[k]
        t = transform_matrices[k]
        
        for id in range(s[0], s[1] + 1):
            item = {}
            item['id'] = id
            item['affine_matrix'] = t
            transforms.append(item)
    
    print(transforms)
    
    with open(opts.output_json_file, 'w') as f:
        json.dump(transforms, f, cls=NumpyArrayEncoder, indent=4, separators=(',', ': '))
            
def main():
    opts = Options()
    ok = opts.process_opts()
    if not ok:
        sys.exit(1)
        
    transform_matrices = load_xg_file(opts.xg_file)
    if not transform_matrices:
        print("Error while loading " + opts.xg_file)
        sys.exit(1)
    print("Loaded xg file transforms, got " + str(len(transform_matrices)) + " entries.")
    print(transform_matrices)
    
    sections = load_sections_file(opts.sections_file)
    if not sections:
        print("Error while loading " + opts.sections_file)
        sys.exit(1)
    print("Loaded sections file, got " + str(len(sections)) + " entries.")
    print(sections)
    
    if len(transform_matrices) != len(sections):
        print("Error: xf and sections file contain different counts of sections.")
        sys.exit(1)
    
    generate_transforms_file(
        opts, transform_matrices, sections)
        
            
if __name__ == '__main__':
    main()         
            