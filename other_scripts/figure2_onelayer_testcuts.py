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


#fig, axs = plt.subplots(nrows=2, ncols=4)


# Choose which cell
x = 6.0
y = 3.0
z=0.08

x1 = 6.03
y1 = 3.85

i1 = 13
i2 = 19

t1 = 1.5

# 30K, LRO
a = nxload('../yan1_recycle/rucl3_yan1_30K_atten.nxs')
dat = a['processed/fft_interp'][z-0.02:z+0.02, y:y+1, x:x+1].sum(axis=0)

# Make averaged piece
av_arr = np.zeros([9,9])
for xoff in [0,8,16]:
    for yoff in [0,8,16]:
        av_arr += dat.data.nxvalue[xoff:xoff+9,yoff:yoff+9]
# Tile averaged piece
tile_av = np.zeros([25,25])
for xoff in [0,8,16]:
    for yoff in [0,8,16]:
        tile_av[xoff:xoff+9,yoff:yoff+9] += av_arr

tile_av[:,8] = tile_av[:,8]/2
tile_av[:,16] = tile_av[:,16]/2
tile_av[8,:] = tile_av[8,:]/2
tile_av[16,:] = tile_av[16,:]/2
delpdf = dat.data.nxvalue - tile_av/9/t1
delpdf = delpdf/np.amax(delpdf)

delpdf = delpdf[i1:i2,:].sum(axis=0)
ax1.plot(dat['x'].nxvalue, delpdf)



# 30K, SRO
a = nxload('../yan1/rucl3_yan1_30K.nxs')
dat = a['processed/fft_interp_p12'][z-0.02:z+0.02, y:y+1, x:x+1].sum(axis=0)
sro_dat = dat.data.nxvalue
sro_dat = sro_dat/np.amax(sro_dat)
sro_dat = sro_dat[i1:i2,:].sum(axis=0)
ax1.plot(dat['x'].nxvalue, sro_dat)
#im2 = ax2.pcolormesh(dat['x'], dat['y'],  sro_dat/np.amax(sro_dat), shading='nearest', cmap='seismic')
#ax1.text(x1,y1,'30K, SRO', fontsize=12)


# 300K, LRO
a = nxload('../yan1_recycle/rucl3_yan1_301K_atten.nxs')
dat = a['processed/fft_interp'][z-0.02:z+0.02, y:y+1, x:x+1].sum(axis=0)

# Make averaged piece
av_arr = np.zeros([9,9])
for xoff in [0,8,16]:
    for yoff in [0,8,16]:
        av_arr += dat.data.nxvalue[xoff:xoff+9,yoff:yoff+9]
# Tile averaged piece
tile_av = np.zeros([25,25])
for xoff in [0,8,16]:
    for yoff in [0,8,16]:
        tile_av[xoff:xoff+9,yoff:yoff+9] += av_arr

tile_av[:,8] = tile_av[:,8]/2
tile_av[:,16] = tile_av[:,16]/2
tile_av[8,:] = tile_av[8,:]/2
tile_av[16,:] = tile_av[16,:]/2
delpdf = dat.data.nxvalue - tile_av/9/t1
delpdf = delpdf/np.amax(delpdf)
delpdf = delpdf[i1:i2,:].sum(axis=0)
ax3.plot(dat['x'].nxvalue, delpdf)
#im3 = ax3.pcolormesh(dat['x'], dat['y'],  delpdf/np.amax(delpdf), shading='nearest', cmap='seismic')
#ax3.text(x1,y1,'300K, LRO', fontsize=12)


# 300K, SRO
a = nxload('../yan1/rucl3_yan1_300K.nxs')
dat = a['processed/fft_interp_p12'][z-0.02:z+0.02, y:y+1, x:x+1].sum(axis=0)
sro_dat = dat.data.nxvalue
sro_dat = sro_dat/np.amax(sro_dat)
sro_dat = sro_dat[i1:i2,:].sum(axis=0)
ax3.plot(dat['x'].nxvalue, sro_dat)
#im4 = ax4.pcolormesh(dat['x'], dat['y'],  sro_dat/np.amax(sro_dat), shading='nearest', cmap='seismic')
#ax4.text(x1,y1,'300K, SRO', fontsize=12)







plt.show()
