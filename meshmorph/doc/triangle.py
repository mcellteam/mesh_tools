#!/usr/bin/python

import random
import math


# triangle points
a = [random.random(),random.random()]
b = [random.random(),random.random()]
c = [random.random(),random.random()]

def diff (a,b):
  m = a[0]-b[0]
  n = a[1]-b[1]
  return [m,n]

def getLength (a,b):
  m = a[0]-b[0]
  n = a[1]-b[1]
  return math.sqrt(m*m+n*n)

def getCentroid (a,b,c):
  # from http://en.wikipedia.org/wiki/Centroid
  p1 = (a[0]+b[0]+c[0])/3.0
  p2 = (a[1]+b[1]+c[1])/3.0
  return [p1,p2]

def cross (a,b):
  return a[0]*b[1]-b[0]*a[1]
  #return (p[1]*v.p[2]-p[2]*v.p[1],
  #        p[2]*v.p[0]-p[0]*v.p[2],
  #        p[0]*v.p[1]-p[1]*v.p[0]);

def getArea (a,b,c):
  # from http://en.wikipedia.org/wiki/Triangle
  # area = 0.5*|ABxAC|
  AB = diff(a,b)
  AC = diff(a,c)
  mycross = cross(AB,AC)
  return math.fabs(0.5*mycross)

def getArea2 (a,b,c):
  m = a[0]-c[0]
  n = b[1]-a[1]
  o = a[0]-b[0]
  p = c[1]-a[1]
  return math.fabs(0.5*(m*n-o*p))

def getMidpoint (a,b):
  x = 0.5*(a[0]+b[0])
  y = 0.5*(a[1]+b[1])
  return [x,y]

print "a = ",a
print "b = ",b
print "c = ",c

#incenter = getIncenter(a,b,c)
#print "incenter = ",incenter

centroid = getCentroid(a,b,c)
print "centroid = ",centroid

total_area = getArea(a,b,c)
print "total area = ",total_area
print "total area2 = ",getArea2(a,b,c)

mAB = getMidpoint(a,b)
mAC = getMidpoint(a,c)
mBC = getMidpoint(b,c)

area_A_mAB_centroid = getArea(a,mAB,centroid)
area_A_mAC_centroid = getArea(a,mAC,centroid)
print "area_A_mAB_centroid = ",area_A_mAB_centroid
print "area_A_mAC_centroid = ",area_A_mAC_centroid
area_B_mAB_centroid = getArea(b,mAB,centroid)
area_B_mBC_centroid = getArea(b,mBC,centroid)
print "area_B_mAB_centroid = ",area_B_mAB_centroid
print "area_B_mBC_centroid = ",area_B_mBC_centroid
area_C_mBC_centroid = getArea(c,mBC,centroid)
area_C_mAC_centroid = getArea(c,mAC,centroid)
print "area_C_mBC_centroid = ",area_C_mBC_centroid
print "area_C_mAC_centroid = ",area_C_mAC_centroid


cumulative_area = area_A_mAB_centroid + area_A_mAC_centroid + area_B_mAB_centroid + area_B_mBC_centroid + area_C_mBC_centroid + area_C_mAC_centroid

print "diff = ",cumulative_area-total_area

# June 2010
# the conclusion was that the triangle is divided
# into six equal pieces, exactly one-third
# of the triangle area to each vertex

# http://en.wikipedia.org/wiki/Triangle
# A median of a triangle is a straight line through a vertex and the
# midpoint of the opposite side, and divides the triangle into two equal
# areas. The three medians intersect in a single point, the triangle's
# centroid or geometric barycenter. The centroid of a rigid triangular
# object (cut out of a thin sheet of uniform density) is also its center
# of mass: the object can be balanced on its centroid in a uniform
# gravitational field. The centroid cuts every median in the ratio 2:1,
# i.e. the distance between a vertex and the centroid is twice the distance
# between the centroid and the midpoint of the opposite side.
