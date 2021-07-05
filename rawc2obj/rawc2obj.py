#!/usr/bin/python
import sys

if __name__ == '__main__':

  if (len(sys.argv)<2):
    print('\nUsage: %s rawc_filename\n'%(sys.argv[0]))
    exit()

  ifn = sys.argv[1]

  rawc_data = [s for s in open(ifn,'r').read().split('\n') if s != '']
  num_verts = int(rawc_data[0].split()[0])
  num_faces = int(rawc_data[0].split()[1])
  
  for i in range(1,num_verts+1):
    vstr = rawc_data[i].split()[:3]
    print('v %s %s %s' %(vstr[0],vstr[1],vstr[2]))

  for i in range(num_verts+1,num_verts+num_faces+1):
    fstr = rawc_data[i].split()[:3]
    print('f %d %d %d' %(1+int(fstr[0]),1+int(fstr[1]),1+int(fstr[2])))


