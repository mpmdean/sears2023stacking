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


# Compute the deltaPDF for a long range ordered structure with respect to the structure averaged onto a 1/3,1/3 unit cell
# Compare the single layer structure for 30K vs 300K; LRO vs SRO


font = {'size' : 10}
matplotlib.rc('font', **font)

fig = plt.figure(figsize=(8,8))
dd = .43
ax1 = fig.add_axes([0.08, 0.08, dd, dd])
ax3 = fig.add_axes([0.08, 0.55, dd, dd])
ax2 = fig.add_axes([0.56, 0.08, dd, dd])
ax4 = fig.add_axes([0.56, 0.55, dd, dd])

# Choose which cell
x = 6.0
y = 3.0
z=0.0

x1 = 6.03
y1 = 3.85
ll = 0.2
ul = 0.46

# 30K, LRO
a = nxload('../yan1_recycle/rucl3_yan1_30K_atten.nxs')
dat = a['processed/fft_interp'][z-0.1:z+0.1, y:y+1, x:x+1]

# Make averaged piece
av_arr = np.zeros([len(dat.z), 9,9])
for xoff in [0,8,16]:
    for yoff in [0,8,16]:
        av_arr += dat.data.nxvalue[:,xoff:xoff+9,yoff:yoff+9]
# Tile averaged piece
tile_av = np.zeros([len(dat.z),25,25])
for xoff in [0,8,16]:
    for yoff in [0,8,16]:
        tile_av[:,xoff:xoff+9,yoff:yoff+9] += av_arr

tile_av[:,:,8] = tile_av[:,:,8]/2
tile_av[:,:,16] = tile_av[:,:,16]/2
tile_av[:,8,:] = tile_av[:,8,:]/2
tile_av[:,16,:] = tile_av[:,16,:]/2
dat.data -= tile_av/9
delpdf = dat[:,y+ll:y+ul,:].sum(axis=1).data.nxdata
im1 = ax1.pcolormesh(dat['x'], dat['z'],  delpdf/np.amax(delpdf), shading='nearest', cmap='seismic')
ax1.text(x1,y1,'30K, LRO', fontsize=12)


# 30K, SRO
a = nxload('../yan1/rucl3_yan1_30K.nxs')
dat = a['processed/fft_interp_p12'][z-0.1:z+0.1, y+ll:y+ul, x:x+1].sum(axis=1)
sro_dat = dat.data.nxvalue
im2 = ax2.pcolormesh(dat['x'], dat['z'],  sro_dat/np.amax(sro_dat), shading='nearest', cmap='seismic')
ax2.text(x1,y1,'30K, SRO', fontsize=12)


# 300K, LRO
a = nxload('../yan1_recycle/rucl3_yan1_301K_atten.nxs')
dat = a['processed/fft_interp'][z-0.1:z+0.1, y:y+1, x:x+1]

# Make averaged piece
av_arr = np.zeros([len(dat.z), 9,9])
for xoff in [0,8,16]:
    for yoff in [0,8,16]:
        av_arr += dat.data.nxvalue[:,xoff:xoff+9,yoff:yoff+9]
# Tile averaged piece
tile_av = np.zeros([len(dat.z),25,25])
for xoff in [0,8,16]:
    for yoff in [0,8,16]:
        tile_av[:,xoff:xoff+9,yoff:yoff+9] += av_arr

tile_av[:,:,8] = tile_av[:,:,8]/2
tile_av[:,:,16] = tile_av[:,:,16]/2
tile_av[:,8,:] = tile_av[:,8,:]/2
tile_av[:,16,:] = tile_av[:,16,:]/2
dat.data -= tile_av/9
delpdf = dat[:,y+ll:y+ul,:].sum(axis=1).data.nxdata
im3 = ax3.pcolormesh(dat['x'], dat['z'],  delpdf/np.amax(delpdf), shading='nearest', cmap='seismic')
ax3.text(x1,y1,'300K, LRO', fontsize=12)


# 300K, SRO
a = nxload('../yan1/rucl3_yan1_300K.nxs')
dat = a['processed/fft_interp_p12'][z-0.1:z+0.1, y+ll:y+ul, x:x+1].sum(axis=1)
sro_dat = dat.data.nxvalue
im4 = ax4.pcolormesh(dat['x'], dat['z'],  sro_dat/np.amax(sro_dat), shading='nearest', cmap='seismic')
ax4.text(x1,y1,'300K, SRO', fontsize=12)




mval = 1.3
im1.set_clim(-mval,mval)
im3.set_clim(-mval,mval)


mval = 1.3
im2.set_clim(-mval,mval)
im4.set_clim(-mval,mval)


ax1.text(6.1,6.9,'z=0')

ax1.set_ylabel('b (real space)')
#ax2.set_ylabel('b (real space)')
ax3.set_ylabel('b (real space)')
#ax4.set_ylabel('b (real space)')
ax1.set_xlabel('a (real space)')
ax2.set_xlabel('a (real space)')
#ax3.set_xlabel('a (real space)')
#ax4.set_xlabel('a (real space)')

"""
ax1.axvline(x=x+.3333,ls='--', color='k', lw=1)
ax2.axvline(x=x+.3333,ls='--', color='k', lw=1)
ax3.axvline(x=x+.3333,ls='--', color='k', lw=1)
ax4.axvline(x=x+.3333,ls='--', color='k', lw=1)

ax1.axvline(x=x+.6666,ls='--', color='k', lw=1)
ax2.axvline(x=x+.6666,ls='--', color='k', lw=1)
ax3.axvline(x=x+.6666,ls='--', color='k', lw=1)
ax4.axvline(x=x+.6666,ls='--', color='k', lw=1)


ax1.axhline(y=y+.3333,ls='--', color='k', lw=1)
ax2.axhline(y=y+.3333,ls='--', color='k', lw=1)
ax3.axhline(y=y+.3333,ls='--', color='k', lw=1)
ax4.axhline(y=y+.3333,ls='--', color='k', lw=1)

ax1.axhline(y=y+.6666,ls='--', color='k', lw=1)
ax2.axhline(y=y+.6666,ls='--', color='k', lw=1)
ax3.axhline(y=y+.6666,ls='--', color='k', lw=1)
ax4.axhline(y=y+.6666,ls='--', color='k', lw=1)
"""
plt.show()
