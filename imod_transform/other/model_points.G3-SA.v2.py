#!/usr/bin/env python3.5

import os
import glob
import swiftir
import numpy as np
import json


# NOTES:
# expecting that the files are named xxx.id.jpeg

class NumpyArrayEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)
    
def append_transform(transforms, img_filename, afm):
    index = img_filename.split('.')[-2]
    t = {
        'id': int(index),
        'affine_matrix': afm
        }
    transforms.append(t)
    


transforms = []

image_dir = 'images_src_s1/'
align_dir = 'images_modeled/'

img_fn_list = sorted(glob.glob(image_dir + '*'))
#print(img_fn_list)
#img_sta = swiftir.loadImage(im_sta_fn)
#img_sta = swiftir.affineImage(ident, im_sta, rect=rect, grayBorder=True)
#swiftir.saveImage(img_sta,ofn)

w = 1464
h = 1493

xmin = 0
xmax = 1080
ymin = 0
ymax = 1080

p_mdl_4 = np.array([[xmin,ymax],
[xmin,ymin],
[xmax,ymin],
[xmax,ymax]])

p_mdl_3 = np.array([[xmin,ymax],
[xmin,ymin],
#[xmax,ymin],
[xmax,ymax]])

xy_offset = [200,300]

p_mdl_4 += xy_offset
p_mdl_3 += xy_offset

p_mdl_4 = ([0,h] - p_mdl_4)*[-1,1]
p_mdl_3 = ([0,h] - p_mdl_3)*[-1,1]

p_mdl_4 = p_mdl_4.transpose()
p_mdl_3 = p_mdl_3.transpose()

#Sections 1-138: [ [left, lower], [left, upper], [right, upper], [right lower] ]
idx1 = range(0,138)
pa1 = np.array([[238, 1392],
              [171, 302],
              [1231, 308],
              [1309, 1411]])
pa1 = ([0,h] - pa1)*[-1,1]


#Section 139-275:
idx2 = range(138,275)
pa2 = np.array([[252, 1368],
[158, 278],
[1221, 309],
[1328, 1373]])
pa2 = ([0,h] - pa2)*[-1,1]


#Section 276-414:
idx3 = range(275,414)
pa3 = np.array([[223, 1368],
[161, 288],
[1206, 275],
[1271, 1365]])
pa3 = ([0,h] - pa3)*[-1,1]


#Section 415-561:
idx4 = range(414,561)
pa4 = np.array([[235, 1419],
[231, 335],
[1256, 305],
[1276, 1399]])
pa4 = ([0,h] - pa4)*[-1,1]


#Section 562-706:
idx5 = range(561,706)
pa5 = np.array([[223, 1315],
[218, 261],
[1246, 229],
[1266, 1291]])
pa5 = ([0,h] - pa5)*[-1,1]


pal = [pa1,pa2,pa3,pa4,pa5]
idxl = [idx1,idx2,idx3,idx4,idx5]

for i in range(len(pal)):
 
  pa = pal[i]
  idx = idxl[i]

  h1 = np.sqrt(( (pa[1][0]-pa[0][0])**2 +(pa[1][1]-pa[0][1])**2 ))
  w1 = np.sqrt(( (pa[2][0]-pa[1][0])**2 +(pa[2][1]-pa[1][1])**2 ))
  h2 = np.sqrt(( (pa[3][0]-pa[2][0])**2 +(pa[3][1]-pa[2][1])**2 ))
  w2 = np.sqrt(( (pa[0][0]-pa[3][0])**2 +(pa[0][1]-pa[3][1])**2 ))
  havg = np.mean([h1,h2])
  wavg = np.mean([w1,w2])
  lavg = np.mean([havg,wavg])
  print('')
  print(str(havg), str(wavg), str(lavg))

  pa = pa.transpose()

  (afm, err, n) = swiftir.mirIterate(p_mdl_4, pa)

  for j in idx: 
    img_fn = img_fn_list[j]
    img = swiftir.loadImage(img_fn)
    append_transform(transforms, img_fn, afm)
    
#    img = swiftir.affineImage(afm, img, rect=rect, grayBorder=True)
    #print("Image " + img_fn + ", applying " + str(afm))
    img = swiftir.affineImage(afm, img)
    ofn = align_dir + os.path.basename(img_fn)
    swiftir.saveImage(img,ofn)


#Section 707-860:
idx6 = range(706,860)
pa6 = np.array([[135, 1116],
[175, 31],
#[--],
[1159, 1041]])
pa6 = ([0,h] - pa6)*[-1,1]


#Section 861-1001:
idx7 = range(860,1001)
pa7 = np.array([[141, 1128],
[166, 75],
#[--],
[1202, 1013]])
pa7 = ([0,h] - pa7)*[-1,1]


pal = [pa6,pa7]
idxl = [idx6,idx7]

for i in range(len(pal)):

  pa = pal[i]
  idx = idxl[i]

  h1 = np.sqrt(( (pa[1][0]-pa[0][0])**2 +(pa[1][1]-pa[0][1])**2 ))
#  w1 = np.sqrt(( (pa[2][0]-pa[1][0])**2 +(pa[2][1]-pa[1][1])**2 ))
#  h2 = np.sqrt(( (pa[3][0]-pa[2][0])**2 +(pa[3][1]-pa[2][1])**2 ))
  w2 = np.sqrt(( (pa[0][0]-pa[2][0])**2 +(pa[0][1]-pa[2][1])**2 ))
#  havg = np.mean([h1,h2])
#  wavg = np.mean([w1,w2])
  lavg = np.mean([h1,w2])
  print('')
  print(str(h1), str(w2), str(lavg))

  pa = pa.transpose()

  (afm, err, n) = swiftir.mirIterate(p_mdl_3, pa)
  print(str(afm))

  for j in idx: 
    img_fn = img_fn_list[j]
    img = swiftir.loadImage(img_fn)
    append_transform(transforms, img_fn, afm)
#    img = swiftir.affineImage(afm, img, rect=rect, grayBorder=True)
    img = swiftir.affineImage(afm, img)
    ofn = align_dir + os.path.basename(img_fn)
    swiftir.saveImage(img,ofn)




'''
Section 1:

238 1392
171 302
1231 308
1309 1411


Section 139:

252 1368
158 278
1221 309
1328 1373


Section 276:

223 1368
161 288
1206 275
1271 1365


Section 415:

235 1419
231 335
1256 305
1276 1399


Section 562:

223 1315
218 261
1246 229
1266 1291


Section 707:

135 1116
175 31
--
1159 1041


Section 861:

141 1128
166 75
--
1202 1013
'''

with open('rev_transforms.json', 'w') as outfile:
    json.dump(transforms, outfile, cls=NumpyArrayEncoder, indent=4, separators=(',', ': '))


