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
dd = .8
ax1 = fig.add_axes([0.08, 0.08, dd, dd])



#fig, axs = plt.subplots(nrows=2, ncols=4)


# Choose which cell
x = 6.
y = 0.0
z=0.0

llim = 2.
ulim = 12.

yscan = 1
# 30K, LRO
a = nxload('../yan1_recycle/rucl3_yan1_30K_atten.nxs')
if yscan:
    dat = a['processed/fft_interp'][z-0.02:z+0.02, llim:ulim, x-0.2:x+0.2].sum(axis=[0,2])
    ax1.plot(dat['y'].nxvalue, dat.data.nxvalue/np.amax(dat.data.nxvalue))
else:
    dat = a['processed/fft_interp'][z-0.02:z+0.02, y-0.2:y+0.2, llim:ulim].sum(axis=[0,1])
    ax1.plot(dat['x'].nxvalue, dat.data.nxvalue/np.amax(dat.data.nxvalue))



# 30K, SRO
a = nxload('../yan1/rucl3_yan1_30K.nxs')
if yscan:
    dat = a['processed/fft_interp_p12'][z-0.02:z+0.02, llim:ulim, x-0.2:x+0.2].sum(axis=[0,2])
    ax1.plot(dat['y'].nxvalue, dat.data.nxvalue/np.amax(dat.data.nxvalue))

else:
    dat = a['processed/fft_interp_p12'][z-0.02:z+0.02, y-0.2:y+0.2, llim:ulim].sum(axis=[0,1])
    ax1.plot(dat['x'].nxvalue, dat.data.nxvalue/np.amax(dat.data.nxvalue))



# 300K, LRO
a = nxload('../yan1_recycle/rucl3_yan1_301K_atten.nxs')
if yscan:
    dat = a['processed/fft_interp'][z-0.02:z+0.02, llim:ulim,x-0.2:x+0.2].sum(axis=[0,2])
    ax1.plot(dat['y'].nxvalue, dat.data.nxvalue/np.amax(dat.data.nxvalue))

else:
    dat = a['processed/fft_interp'][z-0.02:z+0.02, y-0.2:y+0.2, llim:ulim].sum(axis=[0,1])
    ax1.plot(dat['x'].nxvalue, dat.data.nxvalue/np.amax(dat.data.nxvalue))



# 300K, SRO
a = nxload('../yan1/rucl3_yan1_300K.nxs')
if yscan:
    dat = a['processed/fft_interp_p12'][z-0.02:z+0.02, llim:ulim, x-0.2:x+0.2].sum(axis=[0,2])
    ax1.plot(dat['y'].nxvalue, dat.data.nxvalue/np.amax(dat.data.nxvalue))
    
else:
    dat = a['processed/fft_interp_p12'][z-0.02:z+0.02, y-0.2:y+0.2, llim:ulim].sum(axis=[0,1])
    ax1.plot(dat['x'].nxvalue, dat.data.nxvalue/np.amax(dat.data.nxvalue))







plt.show()
