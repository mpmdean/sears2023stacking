#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 16:43:51 2023

@author: jsears
"""
import matplotlib.pyplot as plt
import numpy as np
from nexusformat.nexus import nxload
import matplotlib.colors as colors
import matplotlib
from matplotlib.colors import LogNorm
import matplotlib.transforms as mtransforms
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.ticker as ticker

from ase.io import read

# Load the files 

ciffilename1 = '../../../c2m.cif'
ciffilename2 = '../../../distort_p1.cif'
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
z = 0.24
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
        zv[t1,t1] = vec[2]
        
# Function to place a 3D Gaussian of the appropriate brightness at a given position vec
def place_peak(arr, xarr, yarr, z, vec, weight):
    dist = (xarr-vec[0])**2 + (yarr-vec[1])**2 + (z-vec[2])**2
    sigma = 0.03
    arr += weight*np.exp( -dist/2/sigma**2 )
    



# Make the PDF array
for jj in np.arange(0,len(veclist1)):
    vec = veclist1[jj]
    weight = weightlist1[jj]
    zdist = (vec[-1]-z)
    if abs(zdist)<dz:
        place_peak(pdf_arr, xv,yv,z,vec,weight**2)
        place_peak(pdf_arr, xv+1,yv,z,vec,weight**2)
        place_peak(pdf_arr, xv,yv+1,z,vec,weight**2)
        place_peak(pdf_arr, xv+1,yv+1,z,vec,weight**2)
        place_peak(pdf_arr, xv-1,yv,z,vec,weight**2)
        place_peak(pdf_arr, xv-1,yv-1,z,vec,weight**2)
        place_peak(pdf_arr, xv,yv-1,z,vec,weight**2)
        place_peak(pdf_arr, xv+1,yv-1,z,vec,weight**2)
        place_peak(pdf_arr, xv-1,yv+1,z,vec,weight**2)
    vec = -veclist1[jj]
    zdist = (vec[-1]-z)
    if abs(zdist)<dz:
        place_peak(pdf_arr, xv,yv,z,vec,weight**2)
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
    vec = -vec
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


# Compute the deltaPDF for a long range ordered structure with respect to the structure averaged onto a 1/3,1/3 unit cell
# Compare the single layer structure for 30K vs 300K; LRO vs SRO
# Plot in the shape of the unit cell


#font = {'size' : 10}
#matplotlib.rc('font', **font)

fig = plt.figure(figsize=(1.08*2*(3+3/8)+0.25, 3+3/8+0.125))
ax1 = fig.add_axes([0,0,1,1])
ax1.set_axis_off()

# Choose which cell
x = 6.0
y = 3.0
z=0.333




# 30K, SRO
a = nxload('../../yan1/rucl3_yan1_30K.nxs')
dat = a['processed/fft_interp_p16'][z-0.02:z+0.02, y:y+1, x:x+1].sum(axis=0)
sro_dat = dat.data.nxvalue
im1 = ax1.pcolormesh(dat['x'], dat['y'],  sro_dat/np.amax(sro_dat), shading='nearest', cmap='seismic')
transf = mtransforms.Affine2D().skew_deg(-26.5651,0).translate(0, 0)
trans_data = transf + ax1.transData
im1.set_transform(trans_data)


z1 = 0
dat = a['processed/fft_interp_p16'][z1-0.02:z1+0.02, y:y+1, x:x+1].sum(axis=0)
sro_z0 = dat.data.nxvalue
combo = np.roll(sro_z0,[16,8], axis=(0,1)) + np.roll(sro_z0,[8,16], axis=(0,1))
combo += 0.25*(np.roll(sro_z0,[16,0], axis=(0,1)) + np.roll(sro_z0,[8,8], axis=(0,1))+ np.roll(sro_z0,[0,16], axis=(0,1)))
######### difference
#combo = sro_dat/np.amax(sro_dat) - combo/np.amax(combo)
#im2 = ax1.pcolormesh(dat['x'], dat['y'],  combo, shading='nearest', cmap='seismic')

im2 = ax1.pcolormesh(dat['x'], dat['y'],  combo/np.amax(combo), shading='nearest', cmap='seismic')
transf = mtransforms.Affine2D().skew_deg(-26.5651,0).translate(1.7, 0)
trans_data = transf + ax1.transData
im2.set_transform(trans_data)



# 300K, SRO
a = nxload('../../yan1/rucl3_yan1_300K.nxs')
dat = a['processed/fft_interp_p16'][z-0.02:z+0.02, y:y+1, x:x+1].sum(axis=0)
sro_dat = dat.data.nxvalue
im3 = ax1.pcolormesh(dat['x'], dat['y'],  sro_dat/np.amax(sro_dat), shading='nearest', cmap='seismic')
transf = mtransforms.Affine2D().skew_deg(-26.5651,0).translate(0, 1.6)
trans_data = transf + ax1.transData
im3.set_transform(trans_data)

z1 = 0
dat = a['processed/fft_interp_p16'][z1-0.02:z1+0.02, y:y+1, x:x+1].sum(axis=0)
sro_z0 = dat.data.nxvalue
combo = np.roll(sro_z0,[0,16], axis=(0,1)) #main monoclinic
combo += 0.22*(np.roll(sro_z0,[16,0], axis=(0,1))+ np.roll(sro_z0,[8,8], axis=(0,1))) #secondary monoclinic
combo += 0.17*(np.roll(sro_z0,[8,16], axis=(0,1))+ np.roll(sro_z0,[16,8], axis=(0,1))) #rhombohedral
combo += 0.1*np.roll(sro_z0,[0,0], axis=(0,1)) 
#######difference
#combo = sro_dat/np.amax(sro_dat) - combo/np.amax(combo)
#im4 = ax1.pcolormesh(dat['x'], dat['y'],  combo, shading='nearest', cmap='seismic')


im4 = ax1.pcolormesh(dat['x'], dat['y'],  combo/np.amax(combo), shading='nearest', cmap='seismic')
transf = mtransforms.Affine2D().skew_deg(-26.5651,0).translate(1.7, 1.6)
trans_data = transf + ax1.transData
im4.set_transform(trans_data)




mval = 1.
im1.set_clim(-mval,mval)
im3.set_clim(-mval,mval)


mval = 1.
im2.set_clim(-mval,mval)
im4.set_clim(-mval,mval)



x0 = 1.4896+3
y0 = 4.57916
d=6
ax1.set_xlim(3.,3.+d*1.08)
ax1.set_ylim(2.6,2.6+d*3**-.5)

def drawcell(x,y,ax):
    dx = 1./24
    off = 0.5 +0.02061
    ax.plot([x,x+1+dx,x+1+dx-off,x-off,x], [y, y, y+1+dx, y+1+dx,y], color='k', lw=2)
    ax.plot([x,x], [y, y-dx], color='k', lw=1.5)
    ax.plot([x+1+dx,x+1+dx], [y, y-dx], color='k', lw=1.5)
    ax.plot([x-dx,x], [y, y], color='k', lw=1.5)
    ax.plot([x-dx-off,x-off], [y+1+dx, y+1+dx], color='k', lw=1.5)
    ax.text(x-0.033,y-4*dx,'0')
    ax.text(x-3.5*dx,y-dx,'0')
    ax.text(x+1+dx-0.027,y-4*dx,'1')
    ax.text(x-3.5*dx-off,y+1-0.01,'1')
    return

drawcell(x0,y0,ax1)
drawcell(x0+1.7,y0,ax1)
drawcell(x0+1.7,y0-1.6,ax1)
drawcell(x0,y0-1.6,ax1)

def drawcell1(x,y,ax):
    dx = 1./24
    off = 0.5 +0.02061
    ax.plot([x,x+1+dx,x+1+dx-off,x-off,x], [y, y, y+1+dx, y+1+dx,y], color='k', lw=2)
    ax.plot([x-off*0.3333,x+1+dx-off*0.3333], [y+(1+dx)/3,y+(1+dx)/3], ls='dashed', color='k')
    ax.plot([x-off*0.666,x+1+dx-off*0.666], [y+2*(1+dx)/3,y+2*(1+dx)/3], ls='dashed', color='k')
    ax.plot([x+0.3333+dx/3, x+0.3333-off+dx/2], [y,y+1+dx], ls='dashed', color='k')
    ax.plot([x+0.6666+dx*2./3, x+0.6666-off+dx*2/3], [y,y+1+dx], ls='dashed', color='k')
    ax.plot([x,x], [y, y-dx], color='k', lw=1.5)
    ax.plot([x+1+dx,x+1+dx], [y, y-dx], color='k', lw=1.5)
    ax.plot([x-dx,x], [y, y], color='k', lw=1.5)
    ax.plot([x-dx-off,x-off], [y+1+dx, y+1+dx], color='k', lw=1.5)
    ax.text(x-0.033,y-4*dx,'0')
    ax.text(x-3.5*dx,y-dx,'0')
    ax.text(x+1+dx-0.027,y-4*dx,'1')
    ax.text(x-3.5*dx-off,y+1-0.01,'1')

    return
drawcell1(x0+1.7*2,y0-1.6,ax1)
drawcell1(x0+1.7*2,y0,ax1)

x1 = 1.18+3
ax1.text(x1, 5.7, '300 K, z=1/3')
ax1.text(x1, 5.5-1.4, '30 K , z=1/3')
ax1.text(x1+1.39, 5.7, '300 K, weighted sum z=0')
ax1.text(x1+1.39, 5.5-1.4, '30 K, weighted sum z=0')
ax1.text(x1+1.71*2-.1, 5.7, '300 K stacking')
ax1.text(x1+1.76*2-.15, 5.5-1.4, '30 K stacking')

ax1.text(x0+.47,4.46, 'a')
ax1.text(x0+.47+1.7,4.46, 'a')
ax1.text(x0+.47+1.7,4.26-1.4, 'a')
ax1.text(x0+.47,4.26-1.4, 'a')
ax1.text(x0+.47+1.7*2,4.46, 'a')
ax1.text(x0+.47+1.7*2,4.26-1.4, 'a')

ax1.text(x0-0.4,y0+0.45, 'b')
ax1.text(x0+1.3,y0+0.45, 'b')
ax1.text(x0+1.3,y0-1.15, 'b')
ax1.text(x0-0.4,y0-1.15, 'b')
ax1.text(x0+1.3+1.7,y0+0.45, 'b')
ax1.text(x0+1.3+1.7,y0-1.15, 'b')

da = 0.333*np.array([1,0])
db = 0.333*np.array([-0.5, 1])
vec = da*2
vec *= 0.85
ax1.arrow(x0+1.7*2, y0, vec[0], vec[1], color='g', head_width=0.08, zorder=20)
vec = db*2
vec *= 0.85
ax1.arrow(x0+1.7*2, y0, vec[0], vec[1], color='g', head_width=0.08, zorder=20, linestyle='dotted')
vec = db+da
vec *= 0.75
ax1.arrow(x0+1.7*2, y0, vec[0], vec[1], color='g', head_width=0.08, zorder=20, linestyle='dotted')


vec = da*2+db
vec *= 0.9
ax1.arrow(x0+1.7*2, y0, vec[0], vec[1], color='k', head_width=0.05)
vec = da+2*db
vec *= 0.95
ax1.arrow(x0+1.7*2, y0, vec[0], vec[1], color='k', head_width=0.05)
ax1.plot(x0+1.7*2, y0, marker='o', color='r', ms=3)
ax1.plot(x0+1.7*2, y0, marker='o', mec='r', mfc='None',ms=8)

# Legend arrows
dxl = 0.1
p1 = 8.65
p2 = 5.65
ax1.arrow(p1, p2, dxl, 0,  color='g', head_width=0.03)
ax1.arrow(p1, p2-.12, dxl, 0,  color='g', head_width=0.03, ls='dotted')
ax1.arrow(p1, p2-.24, dxl, 0,  color='k', head_width=0.03)
ax1.arrow(p1, p2-.36, dxl, 0,  color='r', head_width=0.03)
ax1.text(p1+0.18, p2-.03, '53% C2/m')
ax1.text(p1+0.18, p2-0.12-.03, '23% C2/m')
ax1.text(p1+0.18, p2-0.24-.03, '18%, R-3')
ax1.text(p1+0.18, p2-0.36-.03, '5%, no shift')



ax1.arrow(p1, p2-1.6, dxl, 0,  color='g', head_width=0.03)
ax1.arrow(p1, p2-.12-1.6, dxl, 0,  color='k', head_width=0.03)

ax1.text(p1+0.18, p2-.03-1.6, '73%, R-3')
ax1.text(p1+0.18, p2-0.12-.03-1.6, '27%, C2/m')

vec = da*2 +db
vec *= 0.85
ax1.arrow(x0+1.7*2, y0-1.6, vec[0], vec[1], color='g', head_width=0.08)
vec = db*2 +da
vec *= 0.85
ax1.arrow(x0+1.7*2, y0-1.6, vec[0], vec[1], color='g', head_width=0.08)


vec = 2*da
vec *= 0.9
ax1.arrow(x0+1.7*2, y0-1.6, vec[0], vec[1], color='k', head_width=0.05)
vec = 2*db
vec *= 0.95
ax1.arrow(x0+1.7*2, y0-1.6, vec[0], vec[1], color='k', head_width=0.05)

vec = db+da
vec *= 0.9
ax1.arrow(x0+1.7*2, y0-1.6, vec[0], vec[1], color='k', head_width=0.05)

cbaxes = inset_axes(ax1, width="1.1%", height="45%", loc='center left')
cb = fig.colorbar(im1,cax=cbaxes, ticks=[-1.,0,1], orientation='vertical')
cb.set_label('Electron-electron density correl. (Arb. Units)')
cb.set_ticks(ticker.MultipleLocator(0.5))


fig.text( 0.1, 0.9, '(a)')
fig.text(0.36, 0.9, '(c)')
fig.text(0.625, 0.9, '(e)')

fig.text( 0.1, 0.44, '(b)')
fig.text(0.36, 0.44, '(d)')
fig.text(0.625, 0.44, '(f)')

plt.show()
