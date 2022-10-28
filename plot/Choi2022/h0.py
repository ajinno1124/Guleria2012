#!/usr/bin/env python
# -*- coding: utf-8 -*-

from decimal import ROUND_HALF_DOWN
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

d1=pd.read_csv('../../Lambda Parameters.csv',comment='#')

RhoN=np.arange(0,0.6,0.01)

d_KIDS0=d1["a1"][21]*RhoN
d_KIDSA=d1["a1"][22]*RhoN
d_KIDSB=d1["a1"][23]*RhoN
d_KIDSC=d1["a1"][24]*RhoN
d_KIDSD=d1["a1"][25]*RhoN


#print(d1_dense.head)

fig=plt.figure()
subplots_adjust(hspace=0.0,wspace=0.0,top=0.9,left=0.2,right=0.85)
ax = subplot(1,1,1)

ax.plot(RhoN,d_KIDS0,label="KIDS0",ls="-",c='b')
ax.plot(RhoN,d_KIDSA,label="KIDSA",ls="--",c='k')
ax.plot(RhoN,d_KIDSB,label="KIDSB",ls="-.",c='r')
ax.plot(RhoN,d_KIDSC,label="KIDSC",ls="-",c='darkgreen')
ax.plot(RhoN,d_KIDS0,label="KIDSD",ls="-",c='darkorange')


#ax.text(7,0.0050,'SLy4',{'color':'k','fontsize':14})
#ax.text(7,0.0025,'GKW3 (u<1.5)',{'color':'k','fontsize':14})
#ax.text(6,0.0025,'$^{208}_\Lambda$Pb',{'color':'k','fontsize':16})

ax.legend(loc='upper right',frameon=0,numpoints=1,fontsize=14)
ax.set_xlim(0,0.6)
#ax.set_ylim(0,0.031)
#plt.yticks(arange(0.01,0.08,0.02), fontsize=14)
ax.set_xlabel(r'$\rho_N$ (fm$^{-3})$',fontsize=16)
ax.set_ylabel(r'$h_0$ (fm)',fontsize=16)
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
plt.savefig("h0.pdf",dpi=300)
plt.savefig("h0.png",dpi=300)
plt.show()

#指定可能なファイル形式は emf, eps, jpeg, jpg, pdf, png, ps, raw, rgba, svg,
#svgz, tif, tiff です。拡張子を指定すると勝手に判断されます。
