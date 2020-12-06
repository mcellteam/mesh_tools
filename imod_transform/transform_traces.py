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

""" 
example:
python transform_traces.py -s test/G3-SA_original_alignem_project_v3_ident_phys_sec_3.json -r test/ref_interf.xg.json -a test/mito10_backchecked.amod -o test/mito10_backchecked_transformed.amod
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
        self.swift_ir_file = None
        self.rev_transform_file = None
        self.not_invert_rev_transform = False
        self.input_amod_file = ''
        self.output_amod_file = ''
        self.y_crop = 0
        self.x_margin = 0
        self.y_margin = 0
        self.margin_as_first = False

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
            '-n', '--not-invert-rev-transform', action='store_true',
            help='do not invert the reverse transformation suppoled through the -r/--rev-transform argument, default is false'
            )
        parser.add_argument(
            '-a', '--amod', type=str, 
            help='amod file containing definition of traces done on top of the original images')
        parser.add_argument(
            '-y', '--y-crop', type=int, 
            help='how many pixels were cropped from the bottom of the image before alignment with SWIFT-IR was done, default is 0')
        parser.add_argument(
            '-m1', '--margin-as-first', action='store_true', 
            help='x and y margin is applied first then the transformations are done, default is false')
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
            print("Input swift-ir file was not set, this transformation will be ignored.")
        
        if args.rev_transform:
            self.rev_transform_file = args.rev_transform
        else:
            print("Input reverse transforms file was not set, this transformation will be ignored.")
            
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
            
        if args.margin_as_first:
            self.margin_as_first = True

        if args.not_invert_rev_transform:
            self.not_invert_rev_transform = True 

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
    # returns dictionary containing id -> 2x3 numpy array or None
    if not rev_transform_file:
        return None
    
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
    # returns dictionary containing id -> 2x3 numpy array or None
    if not swift_transforms_file:
        return None
    
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

    if rev_transforms and id not in rev_transforms:
        print("Error: reverse transform for layer id " + str(id) + " was not found.")
        sys.exit(1)
    if swift_transforms and id not in swift_transforms:
        print("Error: SWIFT-IR transform for layer id " + str(id) + " was not found.")
        sys.exit(1)

    # transformation was applied on image with reflected y coordinate 
    point[1][0] = height - point[1][0] 
    
    point_margin1 = point
    if opts.margin_as_first:
        point_margin1[0][0] += opts.x_margin
        point_margin1[1][0] -= opts.y_margin
    
    # read selected transformations
    if rev_transforms:
        mrev = np.vstack((rev_transforms[id], [0.0, 0.0, 1.0]))
    else:
        mrev = np.identity(3)
        
    if swift_transforms:
        mswift = np.vstack((swift_transforms[id], [0.0, 0.0, 1.0]))
    else:
        mswift = np.identity(3)
    
    
    if opts.not_invert_rev_transform:
        mrev_inv = mrev
    else:
        mrev_inv = np.linalg.inv(mrev)
    
    mswift_inv = np.linalg.inv(mswift)

    # apply inverse transformation that was applied on 
    # original image to create sources for SWIFT-IR alignment      
    point_rev = mrev_inv.dot(point_margin1)

    # and move by y offset - how much were the images cropped
    point_crop = point_rev 
    point_crop[1][0] += opts.y_crop
    
    # and apply inverse transformation of cummulative 
    # transformation from SWIFT-IR 
    point_swift = mswift_inv.dot(point_crop)
    
    point_margin2 = point_swift
    if not opts.margin_as_first:
        point_margin2[0][0] += opts.x_margin
        point_margin2[1][0] -= opts.y_margin

    # reflect back along Y axis
    point_margin2[1][0] = height - point_margin2[1][0]
    
    return str(point_margin2[0][0]) + ' ' + str(point_margin2[1][0]) + ' ' + str(id) + '\n'
    

# signature for custom_processing_function - (line, opts), returns resultign line
def process_amod_file(opts, rev_transforms, swift_transforms, custom_processing_function=None):
    
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
                    if not custom_processing_function:
                        # transform line
                        res_line = transform_contour_line(
                            line, 
                            rev_transforms, swift_transforms,
                            height,
                            opts
                        )
                    else:
                        res_line = custom_processing_function(
                            line, 
                            opts
                        )
                    fout.write(res_line)
                    remaining_countours -= 1
                 

def check_amod_file(input_amod_file):
    with open(input_amod_file, 'rb') as fin:
        magic = b'IMODV1.2'
        d = fin.read(len(magic))
        if d == magic:
            fname = os.path.splitext(input_amod_file)[0]
            print("Error: input file '" + input_amod_file + 
                  "' seems to be binary .mod file because it starts with magic header " + str(magic)[1:] + ". "
                  "It can be converted with 'imodinfo -a " + input_amod_file + " > " + fname + ".ascii.amod'.")
            return False
    return True

            
def main():
    opts = Options()
    ok = opts.process_opts()
    if not ok:
        sys.exit(1)
        
    rev_transforms = load_rev_transforms(opts.rev_transform_file)
    if rev_transforms:
        print("Loaded reverse transforms, got " + str(len(rev_transforms)) + " entries.")
    
    swift_transforms = load_swift_transforms(opts.swift_ir_file)
    if swift_transforms:
        print("Loaded SWIFT-IR transforms, got " + str(len(swift_transforms)) + " entries.")
    
    if rev_transforms and swift_transforms and rev_transforms.keys() != swift_transforms.keys():
        print("Warning: Reverse and SWIFT-IR transforms contain different image ids!")
        print("Reverse transform ids:" + str(rev_transforms.keys()))
        print("SWIFT-IR transform ids:" + str(swift_transforms.keys()))
    
    ok = check_amod_file(opts.input_amod_file)
    if not ok:
        sys.exit(1)
            
    process_amod_file(
        opts, rev_transforms, swift_transforms)
    
            
if __name__ == '__main__':
    main()            