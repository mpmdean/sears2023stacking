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
from scipy.optimize import leastsq


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
z = 0.333


# 30K, SRO
# Load and plot the z=1/3 data
a = nxload('../../yan1/rucl3_yan1_30K.nxs')
dat = a['processed/fft_interp_p16'][z-0.02:z+0.02, y:y+1, x:x+1].sum(axis=0)
sro_dat = dat.data.nxvalue
im1 = ax1.pcolormesh(dat['x'], dat['y'],  sro_dat/np.amax(sro_dat), shading='nearest', cmap='seismic')
transf = mtransforms.Affine2D().skew_deg(-26.5651,0).translate(0, 0)
trans_data = transf + ax1.transData
im1.set_transform(trans_data)

# Load and plot the z=0 data
z1 = 0
dat = a['processed/fft_interp_p16'][z1-0.02:z1+0.02, y:y+1, x:x+1].sum(axis=0)
sro_z0 = dat.data.nxvalue
im1a = ax1.pcolormesh(dat['x'], dat['y'],  sro_z0/np.amax(sro_z0), shading='nearest', cmap='seismic')
transf = mtransforms.Affine2D().skew_deg(-26.5651, 0).translate(-1.7, 0)
trans_data = transf + ax1.transData
im1a.set_transform(trans_data)

def reduce_dat(dat):
    def box(x,y):
        if x==0:
            if y==0:
                boxval = dat[x-2:, y-2:].ravel().sum() 
                boxval += dat[x-2:, :y+3].ravel().sum()
                boxval += dat[:x+3, y-2:].ravel().sum()
                boxval += dat[:x+3, :y+3].ravel().sum()
            else:
                boxval = dat[x-2:, y-2:y+3].ravel().sum() 
                boxval = dat[:x+3, y-2:y+3].ravel().sum()
        else:
            if y==0:
                boxval = dat[x-2:x+3, y-2:].ravel().sum() 
                boxval = dat[x-2:x+3, :y+3].ravel().sum()   
            else:
                boxval = dat[x-2:x+3, y-2:y+3].ravel().sum() 
        return boxval
    
    datvec = np.array([box(0, 0), box(0, 8), box(0, 16), box(8,0), box(8, 8), box(8, 16), box(16, 0), box(16, 8), box(16, 16) ])
    return datvec
    
# Fit the z=1/3 data set
def costfunc(pars, zerodat, thirddat):
    #sro_z0a = sro_z0[:-1,:-1]
    sro_z0a = zerodat[:-1,:-1]
    combo = np.roll(sro_z0a, [16, 8], axis=(0, 1)) + np.roll(sro_z0a,[8, 16], axis=(0,1))
    combo += pars[0]*(np.roll(sro_z0a,[16,0], axis=(0,1)) + np.roll(sro_z0a,[8, 8], axis=(0,1))+ np.roll(sro_z0a,[0,16], axis=(0,1)))
    combo1 = np.zeros((25,25))
    combo1[:-1,:-1] += combo
    combo1[:-1,-1] += combo[:,0]
    combo1[-1,:-1] += combo[0,:]
    combo1[-1,-1] += combo[0,0]
    combo1 = combo1/np.amax(combo1)
    thirddat = thirddat/np.amax(thirddat)
    cost = reduce_dat(thirddat[:-1,:-1]).ravel() - reduce_dat(combo1[:-1,:-1]).ravel()
    cost = cost**2
    return cost

fit1 = leastsq(costfunc, x0=[0.25], args=(sro_z0, sro_dat), full_output=1)
residuals = costfunc(fit1[0], sro_z0, sro_dat)
hess = fit1[1]
hess_diag = np.asarray([ hess[n][n] for n in np.arange(0,len(hess)) ])
err = hess_diag*np.var(residuals)

# Make the sum of the z=0 data set and plot
sro_z0a = sro_z0[:-1,:-1]
combo = np.roll(sro_z0a, [16, 8], axis=(0, 1)) + np.roll(sro_z0a,[8, 16], axis=(0,1))
combo += fit1[0][0]*(np.roll(sro_z0a,[16,0], axis=(0,1)) + np.roll(sro_z0a,[8, 8], axis=(0,1))+ np.roll(sro_z0a,[0,16], axis=(0,1)))
combo1 = np.zeros((25,25))
combo1[:-1,:-1] += combo
combo1[:-1,-1] += combo[:,0]
combo1[-1,:-1] += combo[0,:]
combo1[-1,-1] += combo[0,0]
im2 = ax1.pcolormesh(dat['x'], dat['y'],  combo1/np.amax(combo1), shading='nearest', cmap='seismic')
transf = mtransforms.Affine2D().skew_deg(-26.5651,0).translate(1.7, 0)
trans_data = transf + ax1.transData
im2.set_transform(trans_data)


# 300K, SRO
# Load and plot the z=1/3 data
a = nxload('../../yan1/rucl3_yan1_300K.nxs')
dat = a['processed/fft_interp_p16'][z-0.02:z+0.02, y:y+1, x:x+1].sum(axis=0)
sro_dat = dat.data.nxvalue
im3 = ax1.pcolormesh(dat['x'], dat['y'],  sro_dat/np.amax(sro_dat), shading='nearest', cmap='seismic')
transf = mtransforms.Affine2D().skew_deg(-26.5651,0).translate(0, 1.6)
trans_data = transf + ax1.transData
im3.set_transform(trans_data)

# Load and plot the z=1/3 data
z1 = 0
dat = a['processed/fft_interp_p16'][z1-0.02:z1+0.02, y:y+1, x:x+1].sum(axis=0)
sro_z0 = dat.data.nxvalue
im3a = ax1.pcolormesh(dat['x'], dat['y'],  sro_z0/np.amax(sro_z0), shading='nearest', cmap='seismic')
transf = mtransforms.Affine2D().skew_deg(-26.5651,0).translate(-1.7, 1.6)
trans_data = transf + ax1.transData
im3a.set_transform(trans_data)

# Fit the z=1/3 data set
def costfunc(pars, zerodat, thirddat):
    sro_z0a = zerodat[:-1,:-1]
    combo = np.roll(sro_z0a,[0,16], axis=(0,1)) #main monoclinic
    combo += pars[0]*(np.roll(sro_z0a,[16,0], axis=(0,1))+ np.roll(sro_z0a,[8,8], axis=(0,1))) #secondary monoclinic
    combo += pars[1]*(np.roll(sro_z0a,[8,16], axis=(0,1))+ np.roll(sro_z0a,[16,8], axis=(0,1))) #rhombohedral
    combo += pars[2]*np.roll(sro_z0a,[0,0], axis=(0,1)) 
    combo1 = np.zeros((25,25))
    combo1[:-1,:-1] += combo
    combo1[:-1,-1] += combo[:,0]
    combo1[-1,:-1] += combo[0,:]
    combo1[-1,-1] += combo[0,0]
    combo1 = combo1/np.amax(combo1)
    thirddat = thirddat/np.amax(thirddat)
    cost = reduce_dat(thirddat[:-1,:-1]).ravel() - reduce_dat(combo1[:-1,:-1]).ravel()
    cost = cost**2
    return cost

fit2 = leastsq(costfunc, x0=[0.22, 0.17, 0.1], args=(sro_z0, sro_dat), full_output=1)
residuals = costfunc(fit2[0], sro_z0, sro_dat)
hess = fit2[1]
hess_diag = np.asarray([ hess[n][n] for n in np.arange(0,len(hess)) ])
err2 = hess_diag*np.var(residuals)

sro_z0a = sro_z0[:-1,:-1]
combo = np.roll(sro_z0a,[0,16], axis=(0,1)) #main monoclinic
combo += fit2[0][0]*(np.roll(sro_z0a,[16,0], axis=(0,1))+ np.roll(sro_z0a,[8,8], axis=(0,1))) #secondary monoclinic
combo += fit2[0][1]*(np.roll(sro_z0a,[8,16], axis=(0,1))+ np.roll(sro_z0a,[16,8], axis=(0,1))) #rhombohedral
combo += fit2[0][2]*np.roll(sro_z0a,[0,0], axis=(0,1)) 
combo1 = np.zeros((25,25))
combo1[:-1,:-1] += combo
combo1[:-1,-1] += combo[:,0]
combo1[-1,:-1] += combo[0,:]
combo1[-1,-1] += combo[0,0]

im4 = ax1.pcolormesh(dat['x'], dat['y'],  combo1/np.amax(combo1), shading='nearest', cmap='seismic')
transf = mtransforms.Affine2D().skew_deg(-26.5651,0).translate(1.7, 1.6)
trans_data = transf + ax1.transData
im4.set_transform(trans_data)




mval = 1.
im1.set_clim(-mval,mval)
im3.set_clim(-mval,mval)
im1a.set_clim(-mval,mval)
im3a.set_clim(-mval,mval)

mval = 1.
im2.set_clim(-mval,mval)
im4.set_clim(-mval,mval)



x0 = 1.4896
y0 = 4.57916
d = 7
ax1.set_xlim(-1,-1+d*1.08)
ax1.set_ylim(2,2+d*3**-.5)

def drawcell(x,y,ax):
    dx = 1./24
    off = 0.5 +0.02061
    ax.plot([x,x+1+dx,x+1+dx-off,x-off,x], [y, y, y+1+dx, y+1+dx,y], color='k', lw=2)
    ax.plot([x,x], [y, y-dx], color='k', lw=1.5)
    ax.plot([x+1+dx,x+1+dx], [y, y-dx], color='k', lw=1.5)
    ax.plot([x-dx,x], [y, y], color='k', lw=1.5)
    ax.plot([x-dx-off,x-off], [y+1+dx, y+1+dx], color='k', lw=1.5)
    ax.text(x-0.037,y-4*dx,'0')
    ax.text(x-3.5*dx,y-dx,'0')
    ax.text(x+1+dx-0.039,y-4*dx,'1')
    ax.text(x-3.5*dx-off,y+1-0.01,'1')
    return

drawcell(x0,y0,ax1)
drawcell(x0+1.7,y0,ax1)
drawcell(x0+1.7,y0-1.6,ax1)
drawcell(x0,y0-1.6,ax1)
drawcell(x0-1.7,y0,ax1)
drawcell(x0-1.7,y0-1.6,ax1)

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
    ax.text(x-0.037,y-4*dx,'0')
    ax.text(x-3.5*dx,y-dx,'0')
    ax.text(x+1+dx-0.037,y-4*dx,'1')
    ax.text(x-3.5*dx-off,y+1-0.01,'1')

    return
drawcell1(x0+3.4,y0-1.6,ax1)
drawcell1(x0+3.4,y0,ax1)

x1 = 1.18
ax1.text(x1-1.7, 5.7, '300 K, z=0')
ax1.text(x1-1.7, 5.5-1.4, '30 K , z=0')
ax1.text(x1-.06, 5.7, '300 K, z=1/3')
ax1.text(x1-.06, 5.5-1.4, '30 K , z=1/3')
ax1.text(x1+1.27, 5.7, '300 K, weighted sum z=0')
ax1.text(x1+1.31, 5.5-1.4, '30 K, weighted sum z=0')
ax1.text(x1+1.71*2-.1, 5.7, '300 K stacking')
ax1.text(x1+1.76*2-.15, 5.5-1.4, '30 K stacking')

ax1.text(x0+.47,4.46, 'a')
ax1.text(x0+.47+1.7,4.46, 'a')
ax1.text(x0+.47+1.7,4.26-1.4, 'a')
ax1.text(x0+.47,4.26-1.4, 'a')
ax1.text(x0+.47+1.7*2,4.46, 'a')
ax1.text(x0+.47+1.7*2,4.26-1.4, 'a')
ax1.text(x0+.47-1.7,4.46, 'a')
ax1.text(x0+.47-1.7,4.26-1.4, 'a')


ax1.text(x0-0.4,y0+0.45, 'b')
ax1.text(x0+1.3,y0+0.45, 'b')
ax1.text(x0+1.3,y0-1.15, 'b')
ax1.text(x0-0.4,y0-1.15, 'b')
ax1.text(x0+1.3+1.7,y0+0.45, 'b')
ax1.text(x0+1.3+1.7,y0-1.15, 'b')
ax1.text(x0-2.1,y0+0.45, 'b')
ax1.text(x0-2.1,y0-1.15, 'b')

da = 0.333*np.array([1,0])
db = 0.333*np.array([-0.5, 1])
vec = da*2
vec *= 0.85
ax1.arrow(x0+2*1.7, y0, vec[0], vec[1], color='g', head_width=0.08, zorder=20)
vec = db*2
vec *= 0.85
ax1.arrow(x0+2*1.7, y0, vec[0], vec[1], color='b', head_width=0.08, zorder=20)
vec = db+da
vec *= 0.75
ax1.arrow(x0+2*1.7, y0, vec[0], vec[1], color='b', head_width=0.08, zorder=20)


vec = da*2+db
vec *= 0.9
ax1.arrow(x0+2*1.7, y0, vec[0], vec[1], color='k', head_width=0.05)
vec = da+2*db
vec *= 0.95
ax1.arrow(x0+2*1.7, y0, vec[0], vec[1], color='k', head_width=0.05)
ax1.plot(x0+2*1.7, y0, marker='o', color='r', ms=3)
ax1.plot(x0+2*1.7, y0, marker='o', mec='r', mfc='None',ms=8)

# Legend arrows
dxl = 0.1
p1 = 5.65
p2 = 5.65
ax1.arrow(p1, p2, dxl, 0,  color='g', head_width=0.03)
ax1.arrow(p1, p2-.12, dxl, 0,  color='b', head_width=0.03)
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
ax1.arrow(x0+2*1.7, y0-1.6, vec[0], vec[1], color='g', head_width=0.08)
vec = db*2 +da
vec *= 0.85
ax1.arrow(x0+2*1.7, y0-1.6, vec[0], vec[1], color='g', head_width=0.08)


vec = 2*da
vec *= 0.9
ax1.arrow(x0+2*1.7, y0-1.6, vec[0], vec[1], color='k', head_width=0.05)
vec = 2*db
vec *= 0.95
ax1.arrow(x0+2*1.7, y0-1.6, vec[0], vec[1], color='k', head_width=0.05)

vec = db+da
vec *= 0.9
ax1.arrow(x0+2*1.7, y0-1.6, vec[0], vec[1], color='k', head_width=0.05)

cbaxes = inset_axes(ax1, width="20%", height="2%", loc='lower center', borderpad=4)
cb = fig.colorbar(im1,cax=cbaxes, ticks=[-1.,0,1], orientation='horizontal')
cb.set_label('Electron-electron density correl. (Arb. Units)')
cb.set_ticks(ticker.MultipleLocator(0.5))


fig.text( 0.005, 0.96, '(a)')
fig.text( 0.23, 0.96, '(c)')
fig.text(0.45, 0.96, '(e)')
fig.text(0.676, 0.96, '(g)')


fig.text( 0.005, 0.55, '(b)')
fig.text( 0.23, 0.55, '(d)')
fig.text(0.45, 0.55, '(f)')
fig.text(0.676, 0.55, '(h)')

plt.show()
