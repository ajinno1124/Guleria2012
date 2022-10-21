#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pylab as pl
from pylab import *
import cmath


#plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['text.usetex'] = True
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['figure.subplot.bottom'] = 0.15

d1_dense=pd.read_csv('../data/Z8N7L1_SLy4GKW3_medium(rho1.5)+Kohno3(k1.0)_a3L30/density.csv',comment='#')
d2_dense=pd.read_csv('../data/Z20N19L1_SLy4GKW3_medium(rho1.5)+Kohno3(k1.0)_a3L30/density.csv',comment='#')
d3_dense=pd.read_csv('../data/Z39N49L1_SLy4GKW3_medium(rho1.5)+Kohno3(k1.0)_a3L30/density.csv',comment='#')
d4_dense=pd.read_csv('../data/Z57N81L1_SLy4GKW3_medium(rho1.5)+Kohno3(k1.0)_a3L30/density.csv',comment='#')
d5_dense=pd.read_csv('../data/Z82N125L1_SLy4GKW3_medium(rho1.5)+Kohno3(k1.0)_a3L30/density.csv',comment='#')

#print(d1_dense.head)

fig=plt.figure()
subplots_adjust(hspace=0.0,wspace=0.0,top=0.9,left=0.2,right=0.85)
ax = subplot(1,1,1)

ax.plot(d1_dense["r(fm)"],(d1_dense["taul"]/d1_dense["Rhol"])**0.5,label=r"$^{16}_\Lambda$O",linewidth=2,color="k")
ax.plot(d2_dense["r(fm)"],(d2_dense["taul"]/d2_dense["Rhol"])**0.5,label=r"$^{40}_\Lambda$Ca",linewidth=2,color="red")
ax.plot(d3_dense["r(fm)"],(d3_dense["taul"]/d3_dense["Rhol"])**0.5,label=r"$^{89}_\Lambda$Y",linewidth=2,color='b')
ax.plot(d4_dense["r(fm)"],(d4_dense["taul"]/d4_dense["Rhol"])**0.5,label=r"$^{139}_\Lambda$La",linewidth=2,color='m')
ax.plot(d5_dense["r(fm)"],(d5_dense["taul"]/d5_dense["Rhol"])**0.5,label=r"$^{208}_\Lambda$Pb",linewidth=2,color='darkgreen')

#ax.plot(d1_dense["r(fm)"],d1_dense["taul"],label=r"$^{16}_\Lambda$O",linewidth=2,color="k")
#ax.plot(d2_dense["r(fm)"],d2_dense["taul"],label=r"$^{40}_\Lambda$Ca",linewidth=2,color="red")
#ax.plot(d3_dense["r(fm)"],d3_dense["taul"],label=r"$^{89}_\Lambda$Y",linewidth=2,color='b')
#ax.plot(d4_dense["r(fm)"],d4_dense["taul"],label=r"$^{139}_\Lambda$La",linewidth=2,color='m')
#ax.plot(d5_dense["r(fm)"],d5_dense["taul"],label=r"$^{208}_\Lambda$Pb",linewidth=2,color='darkgreen')


#ax.text(7,0.0050,'SLy4',{'color':'k','fontsize':14})
#ax.text(7,0.0025,r'GKW3 (u<1.5) + Kohno3 (k<1.0 fm$^{-1}$) $a^\Lambda_3=30$ MeV $\cdot$ fm$^{-5}$)',{'color':'k','fontsize':14})
#ax.text(6,0.0025,'$^{208}_\Lambda$Pb',{'color':'k','fontsize':16})

ax.set_title(r"GKW3 ($u<1.5$) + Kohno3 ($k<1.0$ fm$^{-1}$) $a^\Lambda_3=30$ MeV$\cdot$fm$^{-5}$")
ax.legend(loc='lower right',frameon=0,numpoints=1,fontsize=14)
ax.set_xlim(0,20)
ax.set_ylim(-0.05,1.4)
#plt.yticks(arange(0.01,0.08,0.02), fontsize=14)
ax.set_ylabel(r"$(\tau_\Lambda/\rho_\Lambda)^{1/2}$ (fm$^{-1}$)",fontsize=16)
ax.set_xlabel(r'$r$ (fm)',fontsize=16)
ax.tick_params(axis='x', which='both', top=True, bottom=True, direction='in',labelsize=14)
ax.tick_params(axis='y', which='both', left=True, right=True, direction='in',labelsize=14)
#plt.tight_layout()

#ax.set_xticks([0,1,2,3,4,5])
#ax.set_yticks([-50,0,50,100,200,300])
#ax.tick_params(labelsize=12)

#plt.axis([0.0,1.0, 0.0,0.2])
#plot.ylabel(r'$v_2$',fontsize=20)
#plt.text(1.25,9,'Au+Au b<3.4 fm',{'color':'b','fontsize':20})
#plt.text(0.17,11.5,'$\mu=0.0$',{'color':'b','fontsize':20})
#plt.text(0.17,10,r'hadrons=$(n,\Delta,\pi,\rho$)',{'color':'b','fontsize':20})

#plt.savefig("v2ch.eps",format='eps',dpi=1000)
#plt.savefig("pxe895qmdmsv.eps",format='eps',bbox__inches='tight')
plt.savefig("MomentumLam_GKWa3.pdf",dpi=300)
plt.savefig("MomentumLam_GKWa3.png",dpi=300)
plt.show()

#指定可能なファイル形式は emf, eps, jpeg, jpg, pdf, png, ps, raw, rgba, svg,
#svgz, tif, tiff です。拡張子を指定すると勝手に判断されます。
