# simple script to copy image stack gap files  
import sys
import os
import shutil
from os.path import isfile, join

first = int(sys.argv[1])
last = int(sys.argv[2])

print("Filling gap " + str(first) + " - " + str(last))

blank_file = 'blank.tif'
prefix = 'final_vol.'
suffix = '.tif'

for i in range(first, last+1):
    dst = prefix + str(i).zfill(4) + suffix
    if os.path.exists(dst):
        print("Error: file " + dst + " already exists")
        sys.exit(0)
    shutil.copyfile(blank_file, dst)
        
