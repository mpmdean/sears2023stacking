#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 31 10:43:34 2023

@author: jsears
"""
import matplotlib.pyplot as plt
import numpy as np
from nexusformat.nexus import nxload
import matplotlib.colors as colors
import matplotlib
from matplotlib.colors import LogNorm
import matplotlib.transforms as mtransforms

fig = plt.figure(figsize=(12,6))
ax1 = fig.add_axes([0.1,0.1,.8,.8])

pdf_arr = np.zeros([25,25])
x = np.linspace(0,1,25)
y = np.linspace(0,1,25)
xv, yv = np.meshgrid(x,y)
zv = np.zeros([25,25])

def place_peak(arr, xarr, yarr, z, vec, weight):
    dist = (xarr-vec[0])**2 + (yarr-vec[1])**2 + (z-vec[2])**2
    sigma = 0.03
    arr += weight*np.exp( -dist/2/sigma**2 )
    
vec = np.array([0.5,0.5,0])
weight = 1
z = 0
zdist = (vec[-1]-z)
dz = 0.05
if abs(zdist)<dz:
    place_peak(pdf_arr, xv,yv,z,vec,weight**2)
    

im1 = ax1.pcolormesh(x, y,  pdf_arr, shading='nearest', cmap='seismic')

transf = mtransforms.Affine2D().skew_deg(-26.5651,0).translate(0, 0)
trans_data = transf + ax1.transData
im1.set_transform(trans_data)
ax1.set_xlim(-1,2)
rr = 1./3**.5
ax1.set_ylim(-rr*1,rr*2)
plt.show()