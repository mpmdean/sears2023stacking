#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  5 12:13:01 2023

@author: jsears
"""
import numpy as np

# Atom basis for C2/m
zoff = 0.0
a1 = 0.24/5.98
c1 = 0.15/6.014
ru1 = np.array([0, 0.16653, 0.5]) 
cl1 = np.array([0.2273-a1, 0, 0.7359-zoff])
cl2 = np.array([0.2504-a1*2**-.5, 0.17388+a1*2**-.5, 0.2662+zoff])

# Symmetry operators for C2/m
op0 = np.array([[1,0,0],[0,1,0],[0,0,1]])
op1 = np.array([[-1,0,0],[0,1,0],[0,0,-1]])
op2 = np.array([[-1,0,0],[0,-1,0],[0,0,-1]])
op3 = np.array([[1,0,0],[0,-1,0],[0,0,1]])
oplist = [op0,op1,op2,op3]

# Make the full basis for C2/m, removing duplicates
rulist = [ru1]
cl1list = [cl1]
cl2list = [cl2]

def checknew(vec,veclist):
    found = False
    for v0 in veclist:
        if np.linalg.norm(v0%1-vec%1)<0.001:
            found = True
    return not found

for op in oplist:
    ru_new = op.dot(ru1)
    if checknew(ru_new, rulist):
        rulist.append(ru_new%1)
    cl1_new = op.dot(cl1)
    if checknew(cl1_new, cl1list):
        cl1list.append(cl1_new%1)
    cl2_new = op.dot(cl2)
    if checknew(cl2_new, cl2list):
        cl2list.append(cl2_new%1)
        print(cl2_new)
    else:
        print(cl2_new)

for op in oplist:
    ru_new = op.dot(ru1) + np.array([0.5, 0.5,0])
    if checknew(ru_new, rulist):
        rulist.append(ru_new%1)
    cl1_new = op.dot(cl1)  + np.array([0.5, 0.5,0])
    if checknew(cl1_new, cl1list):
        cl1list.append(cl1_new%1)
    cl2_new = op.dot(cl2)  + np.array([0.5, 0.5,0])
    if checknew(cl2_new, cl2list):
        cl2list.append(cl2_new%1)
        #print(cl2_new)

    else:
        print(cl2_new)

"""        
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

rulist = [c2m_hex.dot(tt) for tt in rulist]
cl1list = [c2m_hex.dot(tt) for tt in cl1list]
cl2list = [c2m_hex.dot(tt) for tt in cl2list]
cl3list = [c2m_hex.dot(tt) for tt in cl3list]
"""
rulist[0] +=  c1*np.array([-1/np.tan(108.8*np.pi/180), 0, 1/np.sin(108.8*np.pi/180)])
rulist[1] -=  c1*np.array([-1/np.tan(108.8*np.pi/180), 0, 1/np.sin(108.8*np.pi/180)])
rulist[2] +=  c1*np.array([-1/np.tan(108.8*np.pi/180), 0, 1/np.sin(108.8*np.pi/180)])
rulist[3] -=  c1*np.array([-1/np.tan(108.8*np.pi/180), 0, 1/np.sin(108.8*np.pi/180)])

n = 1
for ru in rulist:
    print('Ru'+str(n)+' Ru3+ 2 i '+str(ru[0])+' '+str(ru[1])+' '+str(ru[2])+' 0.005 1')
    n+=1

n = 1
for ru in cl1list:
    print('Cl'+str(n)+' Cl1- 2 i '+str(ru[0])+' '+str(ru[1])+' '+str(ru[2])+' 0.005 1')
    n+=1
    
for ru in cl2list:
    print('Cl'+str(n)+' Cl1- 2 i '+str(ru[0])+' '+str(ru[1])+' '+str(ru[2])+' 0.005 1')
    n+=1
