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

d1=pd.read_csv('../JLK.csv',comment='#')
d2=pd.read_csv('../data/BindingEnergyLam/ChiSquared.csv',comment='#')

J=np.zeros(len(d2))
L=np.zeros(len(d2))
K=np.zeros(len(d2))
ms_m=np.zeros(len(d2))

for i in range(len(d2)):
	#Delete KIDS
	#if i<=20 and i>=26:
	J[i]=d1["J (MeV)"][d2["index"][i]-1]
	L[i]=d1["L (MeV)"][d2["index"][i]-1]
	K[i]=d1["K (MeV)"][d2["index"][i]-1]
	ms_m[i]=d1["m*/m"][d2["index"][i]-1]

fig=plt.figure()
subplots_adjust(hspace=0.0,wspace=0.0,top=0.9,left=0.2,right=0.85)
ax = subplot(1,1,1)

ax.scatter(J,ms_m,c=d2["ChiSquare2"])
#cax=plt.axes([0.85,0.1,0.075,0.8])
#plt.colorbar(cax=cax)

#ax.text(7.5,-20,'SLy4',{'color':'k','fontsize':14})
#ax.text(7,-25,r'GKW3(u<1.5)',{'color':'k','fontsize':14})
#ax.text(0.5,-15,'$^{208}_\Lambda$Pb',{'color':'k','fontsize':14})

ax.legend(loc='upper left',frameon=0,numpoints=1,fontsize=14)
#ax.set_xlim(0,10)
#ax.set_ylim(-32,0)
#plt.yticks(arange(0.01,0.08,0.02), fontsize=14)
ax.set_xlabel(r'$J$ (MeV)',fontsize=16)
ax.set_ylabel(r'$m^*/m$ (MeV)',fontsize=16)
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
plt.savefig("Jm_Chi.pdf",dpi=300)
plt.savefig("Jm_Chi.png",dpi=300)
plt.show()

#指定可能なファイル形式は emf, eps, jpeg, jpg, pdf, png, ps, raw, rgba, svg,
#svgz, tif, tiff です。拡張子を指定すると勝手に判断されます。
