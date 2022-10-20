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
plt.rcParams['text.usetex'] = False
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['figure.subplot.bottom'] = 0.15

d1_pot=pd.read_csv('../data/Z8N7L1_SLy4GKW3_medium(rho1.5)/potential.csv',comment='#')
d2_pot=pd.read_csv('../data/Z20N19L1_SLy4GKW3_medium(rho1.5)/potential.csv',comment='#')
d3_pot=pd.read_csv('../data/Z39N49L1_SLy4GKW3_medium(rho1.5)/potential.csv',comment='#')
d4_pot=pd.read_csv('../data/Z57N81L1_SLy4GKW3_medium(rho1.5)/potential.csv',comment='#')
d5_pot=pd.read_csv('../data/Z82N125L1_SLy4GKW3_medium(rho1.5)/potential.csv',comment='#')
#print(d1_dense.head)

fig=plt.figure()
subplots_adjust(hspace=0.0,wspace=0.0,top=0.9,left=0.2,right=0.85)
ax = subplot(1,1,1)

ax.plot(d1_pot["r(fm)"],d1_pot["VNp(MeV)"],label=r"$^{16}_\Lambda$O",linewidth=2,color="k", ls="-")
ax.plot(d2_pot["r(fm)"],d2_pot["VNp(MeV)"],label=r"$^{40}_\Lambda$Ca",linewidth=2,color="red", ls="-")
ax.plot(d3_pot["r(fm)"],d3_pot["VNp(MeV)"],label=r"$^{89}_\Lambda$Y",linewidth=2,color='b',ls='-')
ax.plot(d4_pot["r(fm)"],d4_pot["VNp(MeV)"],label=r"$^{139}_\Lambda$La",linewidth=2,color='m',ls='-')
ax.plot(d5_pot["r(fm)"],d5_pot["VNp(MeV)"],label=r"$^{208}_\Lambda$Pb",linewidth=2,color='darkgreen',ls='-')

ax.text(7.5,-60,'SLy4',{'color':'k','fontsize':14})
ax.text(7,-70,r'GKW3(u<1.5)',{'color':'k','fontsize':14})
#ax.text(0.5,-15,'$^{208}_\Lambda$Pb',{'color':'k','fontsize':14})

ax.legend(loc='upper left',frameon=0,numpoints=1,fontsize=14)
ax.set_xlim(0,10)
ax.set_ylim(-85,0)
#plt.yticks(arange(0.01,0.08,0.02), fontsize=14)
ax.set_ylabel('Potential (MeV)',fontsize=16)
ax.set_xlabel('$r$ (fm)',fontsize=16)
ax.tick_params(axis='x', which='both', direction='in',labelsize=14)
ax.tick_params(axis='y', which='both', direction='in',labelsize=14)
plt.tight_layout()

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
plt.savefig("PotentialN_GKW3.pdf",dpi=300)
plt.savefig("PotentialN_GKW3.png",dpi=300)
plt.show()

#指定可能なファイル形式は emf, eps, jpeg, jpg, pdf, png, ps, raw, rgba, svg,
#svgz, tif, tiff です。拡張子を指定すると勝手に判断されます。
