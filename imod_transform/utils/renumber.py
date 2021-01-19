# simple script to renumber image files  
import sys
import os
from os.path import isfile, join

first = int(sys.argv[1])
last = int(sys.argv[2])
inc = int(sys.argv[3])

print("Renumbering " + str(first) + " - " + str(last) + "(included), addend " + str(inc))


files = [f for f in os.listdir('.') if isfile(f)]

prefix = 'final_vol.'
suffix = '.tif'

for f in sorted(files, reverse=True):
    if f.startswith(prefix) and f.endswith(suffix):
        index_str = f[len(prefix):-len(suffix)]
        index = int(index_str)
        if index >= first and index <= last:
            new_index_str = str(index + inc).zfill(4)
            new_name = prefix + new_index_str + suffix
            print("Renaming " + f + " to " + new_name)
            os.rename(f, new_name)
        
