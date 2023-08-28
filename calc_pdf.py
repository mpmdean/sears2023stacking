#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  2 17:15:09 2023

@author: jsears
"""

# Compute the PDF for RuCl3, original and distorted structures

import matplotlib.pyplot as plt
import numpy as np
from ase.io import read

# Load the files 

ciffilename1 = '../../c2m.cif'
ciffilename2 = '../../distort_p1.cif'
struc1 = read(ciffilename1, format="cif")
struc2 = read(ciffilename2, format="cif")



# Get the structural information
pos1 = struc1.get_scaled_positions()
pos2 = struc2.get_scaled_positions()

symbol1 = struc1.get_chemical_symbols()
symbol2 = struc2.get_chemical_symbols()

weights1 = struc1.get_atomic_numbers()
weights2 = struc2.get_atomic_numbers()



# Make list of all difference vectors
veclist1 = []
atomlist1 = []
weightlist1 = []
for ii in np.arange(0,len(pos1)):
    for jj in np.arange(0,len(pos1)):
        veclist1.append( (pos1[ii]-pos1[jj])%1 )        
        atomlist1.append([symbol1[ii], symbol1[jj] ])
        weightlist1.append(weights1[ii]*weights1[jj])
        
veclist2 = []
atomlist2 = []
weightlist2 = []
for ii in np.arange(0,len(pos2)):
    for jj in np.arange(0,len(pos2)):
        veclist2.append( (pos2[ii]-pos2[jj])%1 )        
        atomlist2.append([symbol2[ii], symbol2[jj] ])
        weightlist2.append(weights2[ii]*weights2[jj])
        
# Choose the z and dz values for the cut
z = .0
dz = .05



# Change basis to hexagonal unit cell
a = np.array([5.98, 0, 0])
b = np.array([0, 10.357, 0])
c = np.array([np.cos(108.8*np.pi/180), 0, np.sin(108.8*np.pi/180)])*6.014
c2m_orth = np.array([a,b,c]).transpose()

a1 = np.array([5.98, 0, 0])
b1 = 5.98*np.array([np.cos(np.pi*120/180), np.sin(np.pi*120/180), 0])
c1 = np.array([0, 0, 17.05/3])
vol1 = a1.dot(np.cross(b1,c1))
ast = np.cross(b1,c1)/vol1
bst = np.cross(c1,a1)/vol1
cst = np.cross(a1,b1)/vol1
tt = np.array([ast,bst,cst])
metric = tt.dot(tt.transpose())

orth_hex = np.linalg.inv( np.array([a1,b1,c1]).transpose() )

c2m_hex = orth_hex.dot(c2m_orth)
hex_c2m = np.linalg.inv(c2m_hex)
# Make arrays in hexagonal coordinates
pdf_arr = np.zeros([25,25])
x = np.linspace(0,1,25)
y = np.linspace(0,1,25)
xv, yv = np.meshgrid(x,y)
zv = np.zeros([25,25])
# Convert the coordinate arrays to monoclinic coordinates

for t1 in np.arange(0,len(x)):
    for t2 in np.arange(0,len(y)):
        vec = np.array([x[t1], y[t2], z])
        vec = hex_c2m.dot(vec)
        xv[t1,t2] = vec[0]
        yv[t1,t2] = vec[1]
        zv[t1,t2] = vec[2]
       
# Function to place a 3D Gaussian of the appropriate brightness at a given position vec
def place_peak(arr, xarr, yarr, z, vec, weight):
    dist = (xarr-vec[0])**2 + 1.7319**2*(yarr-vec[1])**2 + (z-vec[2])**2
    sigma = 0.03
    arr += weight*np.exp( -dist/2/sigma**2 )
    



# Make the PDF array
for jj in np.arange(0,len(veclist1)):
    vec = veclist1[jj]
    weight = weightlist1[jj]
    zdist = (vec[-1]-z)
    if abs(zdist)<dz:
        place_peak(pdf_arr, xv,yv,z,vec,weight**2)
        #dz=-1e5
        place_peak(pdf_arr, xv+1,yv,z,vec,weight**2)
        place_peak(pdf_arr, xv,yv+1,z,vec,weight**2)
        place_peak(pdf_arr, xv+1,yv+1,z,vec,weight**2)
        place_peak(pdf_arr, xv-1,yv,z,vec,weight**2)
        place_peak(pdf_arr, xv-1,yv-1,z,vec,weight**2)
        place_peak(pdf_arr, xv,yv-1,z,vec,weight**2)
        place_peak(pdf_arr, xv+1,yv-1,z,vec,weight**2)
        place_peak(pdf_arr, xv-1,yv+1,z,vec,weight**2)
       
pdf_arr2 = np.zeros([25,25])

# Make the PDF array
for jj in np.arange(0,len(veclist2)):
    vec = veclist2[jj]
    weight = weightlist2[jj]
    zdist = (vec[-1]-z)
    if abs(zdist)<dz:
        place_peak(pdf_arr2, xv,yv,z,vec,weight**2)
        place_peak(pdf_arr2, xv+1,yv,z,vec,weight**2)
        place_peak(pdf_arr2, xv,yv+1,z,vec,weight**2)
        place_peak(pdf_arr2, xv+1,yv+1,z,vec,weight**2)
        place_peak(pdf_arr2, xv-1,yv,z,vec,weight**2)
        place_peak(pdf_arr2, xv-1,yv-1,z,vec,weight**2)
        place_peak(pdf_arr2, xv,yv-1,z,vec,weight**2)
        place_peak(pdf_arr2, xv+1,yv-1,z,vec,weight**2)
        place_peak(pdf_arr2, xv-1,yv+1,z,vec,weight**2)
        
# Make the plots using dots to represent intensity
fig = plt.figure(figsize=(5,10))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
ax1.set_title('Structure 1')
ax2.set_title('Structure 2')
tt = 3e-3
im2 = ax1.pcolormesh(x,y,pdf_arr)
#import matplotlib.transforms as mtransforms
#transf = mtransforms.Affine2D().skew_deg(-30,0)
#trans_data = transf + ax1.transData
#im2.set_transform(trans_data)

for jj in np.arange(0,len(veclist1)):
    vec = c2m_hex.dot(veclist1[jj])
    #vec = veclist1[jj]
    intensity = weightlist1[jj]
    zdist = (vec[-1]-z)
    if abs(zdist)<dz:
        ax1.plot(vec[0], vec[1], marker='o', ms=tt*intensity, color='r')
        ax1.plot(vec[0]+1, vec[1], marker='o', ms=tt*intensity, color='r')
        ax1.plot(vec[0], vec[1]+1, marker='o', ms=tt*intensity, color='r')
        ax1.plot(vec[0]+1, vec[1]+1, marker='o', ms=tt*intensity, color='r')
        ax1.plot(vec[0]-1, vec[1], marker='o', ms=tt*intensity, color='r')
        ax1.plot(vec[0], vec[1]-1, marker='o', ms=tt*intensity, color='r')
        ax1.plot(vec[0]-1, vec[1]-1, marker='o', ms=tt*intensity, color='r')
        ax1.plot(vec[0]-1, vec[1]+1, marker='o', ms=tt*intensity, color='r')
        ax1.plot(vec[0]+1, vec[1]-1, marker='o', ms=tt*intensity, color='r')
    vec = -c2m_hex.dot(veclist1[jj])
    #vec = veclist1[jj]
    intensity = weightlist1[jj]
    zdist = (vec[-1]-z)
    if abs(zdist)<dz:
        ax1.plot(vec[0], vec[1], marker='o', ms=tt*intensity, color='r')
        ax1.plot(vec[0]+1, vec[1], marker='o', ms=tt*intensity, color='r')
        ax1.plot(vec[0], vec[1]+1, marker='o', ms=tt*intensity, color='r')
        ax1.plot(vec[0]+1, vec[1]+1, marker='o', ms=tt*intensity, color='r')
        ax1.plot(vec[0]-1, vec[1], marker='o', ms=tt*intensity, color='r')
        ax1.plot(vec[0], vec[1]-1, marker='o', ms=tt*intensity, color='r')
        ax1.plot(vec[0]-1, vec[1]-1, marker='o', ms=tt*intensity, color='r')
        ax1.plot(vec[0]-1, vec[1]+1, marker='o', ms=tt*intensity, color='r')
        ax1.plot(vec[0]+1, vec[1]-1, marker='o', ms=tt*intensity, color='r')


ax2.pcolormesh(x,y,pdf_arr2)

for jj in np.arange(0,len(veclist2)):
    vec = c2m_hex.dot(veclist2[jj])
    intensity = weightlist2[jj]
    zdist = (vec[-1]-z)
    if abs(zdist)<dz:
        ax2.plot(vec[0], vec[1], marker='o', ms=tt*intensity, color='r')
        ax2.plot(vec[0]+1, vec[1], marker='o', ms=tt*intensity, color='r')
        ax2.plot(vec[0], vec[1]+1, marker='o', ms=tt*intensity, color='r')
        ax2.plot(vec[0]+1, vec[1]+1, marker='o', ms=tt*intensity, color='r')
        ax2.plot(vec[0]-1, vec[1], marker='o', ms=tt*intensity, color='r')
        ax2.plot(vec[0], vec[1]-1, marker='o', ms=tt*intensity, color='r')
        ax2.plot(vec[0]-1, vec[1]-1, marker='o', ms=tt*intensity, color='r')
        ax2.plot(vec[0]-1, vec[1]+1, marker='o', ms=tt*intensity, color='r')
        ax2.plot(vec[0]+1, vec[1]-1, marker='o', ms=tt*intensity, color='r')
    vec = -c2m_hex.dot(veclist2[jj])
    intensity = weightlist2[jj]
    zdist = (vec[-1]-z)
    if abs(zdist)<dz:
        ax2.plot(vec[0], vec[1], marker='o', ms=tt*intensity, color='r')
        ax2.plot(vec[0]+1, vec[1], marker='o', ms=tt*intensity, color='r')
        ax2.plot(vec[0], vec[1]+1, marker='o', ms=tt*intensity, color='r')
        ax2.plot(vec[0]+1, vec[1]+1, marker='o', ms=tt*intensity, color='r')
        ax2.plot(vec[0]-1, vec[1], marker='o', ms=tt*intensity, color='r')
        ax2.plot(vec[0], vec[1]-1, marker='o', ms=tt*intensity, color='r')
        ax2.plot(vec[0]-1, vec[1]-1, marker='o', ms=tt*intensity, color='r')
        ax2.plot(vec[0]-1, vec[1]+1, marker='o', ms=tt*intensity, color='r')
        ax2.plot(vec[0]+1, vec[1]-1, marker='o', ms=tt*intensity, color='r')
       
ax1.set_xlim(-0.1,1.1)
ax2.set_xlim(-0.1,1.1)
ax1.set_ylim(-0.1,1.1)
ax2.set_ylim(-0.1,1.1)

plt.show()
