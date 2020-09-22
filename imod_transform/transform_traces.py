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

class Options:
    def __init__(self):
        self.swift_ir_file = ''
        self.rev_transform_file = ''
        self.input_amod_file = ''
        self.output_amod_file = ''
        self.y_crop = 0
        self.x_margin = 0
        self.y_margin = 0

    def __repr__(self):
        attrs = vars(self)
        return ", ".join("%s: %s" % item for item in attrs.items())
            
    @staticmethod
    def create_argparse():
        parser = argparse.ArgumentParser(description='Traces transformation tool')
        parser.add_argument(
            '-s', '--swift-ir', type=str, 
            help='json file produced by SWIFT-IR that contains information on transformations done during alignment')
        parser.add_argument(
            '-r', '--rev-transform', type=str, 
            help='json file containing information on transformation done on original images before passed onto SWIFT-IR'
                 ', it must be supplied as a reverse transformation, i.e. how to get from SWIFT-IR images to original images'
            )
        parser.add_argument(
            '-a', '--amod', type=str, 
            help='amod file containing definition of traces done on top of the original images')
        parser.add_argument(
            '-y', '--y-crop', type=int, 
            help='how many pixels were cropped from the bottom of the image before alignment with SWIFT-IR was done, default is 0')
        parser.add_argument(
            '-mx', '--x-margin', type=int, 
            help='x-margin used when generating new images with SWIFT-IR, default is 0')
        parser.add_argument(
            '-my', '--y-margin', type=int, 
            help='y-margin used when generating new images with SWIFT-IR, default is 0')
        parser.add_argument(
            '-o', '--output', type=str, 
            help='name of output .amod file')
        return parser

    def process_opts(self):
        
        parser = Options.create_argparse()
        args = parser.parse_args()
        
        if args.swift_ir:
            self.swift_ir_file = args.swift_ir
        else:
            print("Input swift-ir file must be set")
            return False
        
        # optional?
        if args.rev_transform:
            self.rev_transform_file = args.rev_transform
        else:
            print("Input reverse transforms file must be set")
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

        if args.y_crop:
            self.y_crop = args.y_crop 
    
        if args.x_margin:
            self.x_margin = args.x_margin 

        if args.y_margin:
            self.y_margin = args.y_margin 

        return True

        
def get_dict_value(d, keys, fname):
    res = d
    for k in keys:
        if k not in res:
            print("Error: did not find expected JSON key '" + k + "' in hierarchy " +
                  str(keys) + " in file '" + fname + "'.")
            sys.exit(1)
        res = res[k]
    return res
                  
        
def load_rev_transforms(rev_transform_file):
    # returns dictionary containing id -> 2x3 numpy array
    res = {}
    with open(rev_transform_file, 'r') as fin:
        json_transforms = json.load(fin)
        
        for t in json_transforms:
            id = get_dict_value(t, ['id'], rev_transform_file)
            atm = np.array(get_dict_value(t, ['affine_matrix'], rev_transform_file))
            assert atm.shape == (2, 3)
            res[id] = atm
        return res
    return None


def load_swift_transforms(swift_transforms_file):
    # returns dictionary containing id -> 2x3 numpy array
    res = {}
    with open(swift_transforms_file, 'r') as fin:
        json_swift = json.load(fin)
        
        alignment_stack = get_dict_value(
            json_swift, ['data', 'scales', 'scale_1', 'alignment_stack'], swift_transforms_file)
        
        for align in alignment_stack:

            atm = np.array(get_dict_value(
                align, ['align_to_ref_method', 'method_results', 'cumulative_afm'], swift_transforms_file))
            assert atm.shape == (2, 3)
            
            fname = get_dict_value(align, ['images', 'base', 'filename'], swift_transforms_file)
            try:
                id = int(fname.split('.')[-2])
            except:
                print("Error: could not get layer id from filename '" + fname + "'.")
                sys.exit(1)
                
            res[id] = atm
        return res
    return None



def transform_contour_line(line, rev_transforms, swift_transforms, height, opts):
    data = line.split()
    assert len(data) == 3
    id = int(data[2])
    point = np.array([[float(data[0])], [float(data[1])], [1.0]])

    if id not in rev_transforms:
        print("Error: reverse transform for layer id " + str(id) + " was not found.")
        sys.exit(1)
    if id not in swift_transforms:
        print("Error: SWIFT-IR transform for layer id " + str(id) + " was not found.")
        sys.exit(1)

    # transformation was applied on image with reflected y coordinate 
    point[1][0] = height - point[1][0] 
    
    # apply both transformations
    mrev = np.vstack((rev_transforms[id], [0.0, 0.0, 1.0]))
    mswift = np.vstack((swift_transforms[id], [0.0, 0.0, 1.0]))
    
    mrev_inv = np.linalg.inv(mrev)
    mswift_inv = np.linalg.inv(mswift)

    # apply inverse transformation that was applied on 
    # original image to create sources for SWIFT-IR alignment      
    point_rev = mrev_inv.dot(point)

    # and move by y offset - how much were the images cropped
    point_crop = point_rev 
    point_crop[1][0] += opts.y_crop
    
    # and apply inverse transformation of cummulative 
    # transformation from SWIFT-IR 
    point_swift = mswift_inv.dot(point_crop)

    point_margin = point_swift
    point_margin[0][0] += opts.x_margin
    point_margin[1][0] -= opts.y_margin

    # reflect back along Y axis
    point_margin[1][0] = height - point_margin[1][0]
    
    return str(point_margin[0][0]) + ' ' + str(point_margin[1][0]) + ' ' + str(id) + '\n'
    

def process_amod_file(opts, rev_transforms, swift_transforms):
    
    with open(opts.input_amod_file, 'r') as fin:
        with open(opts.output_amod_file, 'w') as fout:
            remaining_countours = 0
            height = -1
            for line in fin:
                if remaining_countours == 0:
                    if 'max' in line:
                        items = line.split()
                        assert len(items) == 4
                        height = int(items[2])                   
                        
                    if 'contour' in line:
                        items = line.split()
                        assert len(items) == 4
                        remaining_countours = int(items[3])
                    
                    # copy line verbatim
                    fout.write(line)
                else:
                    # transform line
                    res_line = transform_contour_line(
                        line, 
                        rev_transforms, swift_transforms,
                        height,
                        opts
                    )
                    fout.write(res_line)
                    remaining_countours -= 1
                 
            
def main():
    opts = Options()
    ok = opts.process_opts()
    if not ok:
        sys.exit(1)
        
    rev_transforms = load_rev_transforms(opts.rev_transform_file)
    if not rev_transforms:
        print("Error while loading " + opts.rev_transform_file)
        sys.exit(1)
    print("Loaded reverse transforms, got " + str(len(rev_transforms)) + " entries.")
    
    swift_transforms = load_swift_transforms(opts.swift_ir_file)
    if not swift_transforms:
        print("Error while loading " + opts.swift_ir_file)
        sys.exit(1)
    print("Loaded SWIFT-IR transforms, got " + str(len(swift_transforms)) + " entries.")
    
    if rev_transforms.keys() != swift_transforms.keys():
        print("Warning: Reverse and SWIFT-IR transforms contain different image ids!")
        print("Reverse transform ids:" + str(rev_transforms.keys()))
        print("SWIFT-IR transform ids:" + str(swift_transforms.keys()))
              
    process_amod_file(
        opts, rev_transforms, swift_transforms)
    
            
if __name__ == '__main__':
    main()            