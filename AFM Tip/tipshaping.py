#!/usr/bin/env python3

import math
import MDAnalysis as mda
#import nglview as nv
import numpy as np
from IPython.core.display import Image
import MDAnalysis.analysis.rdf
import matplotlib.pyplot as plt
from MDAnalysis.analysis import distances

u = mda.Universe('Tip20nm_RUN.gro')

SiO2 = u.select_atoms('name Si')
minx, miny, minz = np.min(SiO2.positions,0) #find the minimum coordinates of the box
maxx, maxy, maxz = np.max(SiO2.positions,0) #find the maximum coordinates of the box

#A = np.array([(minx+maxx)/2, (miny+maxy)/2,minz])
#B = np.array([minx, miny, maxz])
#C = np.array([minx, maxy, maxz])
#D = np.array([maxx, miny, maxz])
#E = np.array([maxx, maxy, maxz])

#setting the co-ordinates of all the vertices
A = np.array([(minx+maxx)/2, (miny+maxy)/2,minz])
B = np.array([minx, miny, maxz])
C = np.array([minx, maxy, maxz])
E = np.array([maxx, miny, maxz])
D = np.array([maxx, maxy, maxz])

#crossproduct to find the normal onto the plane
nabc = -1*np.cross(C - A, B - A)
nacd = -1*np.cross(D - A, C - A)
nade = -1*np.cross(E - A, D - A)
naeb = -1*np.cross(B - A, E - A)
nbcd = -1*np.cross(C - B, D - B)

inside_or_not = np.zeros([len(SiO2.positions),1]) #create a zero array to store beads that are inside the pyramid

print('A:' + str(A)) #print all coordinates
print('B:' + str(B))
print('C:' + str(C))
print('D:' + str(D))
print('E:' + str(E))

#dotproduct to find if the bead is inside or outside of the pyramid
for ii in range(len(SiO2.positions)):
    point = SiO2.positions[ii]
    r1 = np.dot(point - A, nabc)
    r2 = np.dot(point - A, nacd)
    r3 = np.dot(point - A, nade)
    r4 = np.dot(point - A, naeb)
    r5 = np.dot(point - B, nbcd)
    
    all_rs = np.array([r1,r2,r3,r4,r5])
    if np.max(all_rs)<0: #if all the dotproducts are negative then the bead is inside the pyramid or else outside the pyramid
        inside_or_not[ii]=1
    else:
        inside_or_not[ii]=0

print(np.sum(inside_or_not))

# this is an array of zeros, that will keep track of all the positions that are inside the array
centers_in_pyramid = np.zeros([int(np.sum(inside_or_not)),3])
iiii=0 # a new counter
for ii in range(len(SiO2.positions)): #go through all original SiO2 positions
    if inside_or_not[ii]==1: #if that SiO2 position is inside the pyramid
        centers_in_pyramid[iiii,:] = SiO2.positions[ii] #store that SiO2 position in the new rray
        iiii +=1 #increase the counter for the new array

Si_selection = u.select_atoms('name Si') #list of silicon atoms with position
Si_selection = Si_selection[0:int(np.sum(inside_or_not))] #select only the first N SiO2 atoms, because those are the ones
#whose positions we will replace. Here N is number of atoms inside the pyramid
Si_selection.positions = centers_in_pyramid #replace the positions with new positions
print(Si_selection.positions)

Si_selection.write('Sqaure_pyramid_tip.gro')
