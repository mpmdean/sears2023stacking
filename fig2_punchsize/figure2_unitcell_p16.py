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
plt.rcParams['text.usetex'] = False

# Diffuse scattering deltaPDF, 30K and 300K
# Plot the z=1/3 cut encoding the interlayer shifts
# Compare with weighted sum of z=0 cut, shifted by different amounts

#font = {'size' : 10}
#matplotlib.rc('font', **font)

fig = plt.figure(figsize=(1.08*2*(3+3/8)+0.25, 3+3/8+0.125))
ax1 = fig.add_axes([0,0,1,1])
ax1.set_axis_off()

# Choose which cell
x = 3.0
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



x0 = 1.4896
y0 = 4.57916
d=6
ax1.set_xlim(0.,0.+d*1.08)
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

x1 = 1.18
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
p1 = 5.65
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
