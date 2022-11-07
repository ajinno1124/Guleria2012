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
#plt.rc('text.latex', preamble=r'\usepackage{braket}')

d1=pd.read_csv('../JLK.csv',comment='#')
d2=pd.read_csv('../data/BindingEnergyLam/ChiSquared.csv',comment='#')
d3=d2[d2["ChiSquare2"]<1.5]
#d3=d2.query("ChiSquare2 < 3")
#print(d3["index"])

#J=np.zeros(len(d3))
#L=np.zeros(len(d3))
#K=np.zeros(len(d3))
#ms_m=np.zeros(len(d3))

#for index in d3["index"]:
	#Delete KIDS\
	#if i<=20 and i>=26:
	#if d3["ChiSquare2"][i]<3:
	#J[i]=d1["J (MeV)"][index-1]
	#L[i]=d1["L (MeV)"][index-1]
	#K[i]=d1["K (MeV)"][index-1]
	#ms_m[i]=d1["m*/m"][index-1]

fig=plt.figure()
subplots_adjust(hspace=0.0,wspace=0.0,top=0.9,left=0.2,right=0.85)
ax = subplot(1,1,1)

#ax.scatter(J,ms_m,c=d3["ChiSquare2"])
#ax.scatter(J,ms_m,c=d2["ChiSquare3"])
plt.scatter(d1["J (MeV)"][d3["index"]-1],d1["K (MeV)"][d3["index"]-1],c=d3["ChiSquare2"],cmap=plt.cm.jet,s=30,edgecolor='k')
cbar=plt.colorbar(aspect=40,pad=0.08,orientation='vertical')
#cbar.set_label(r"$\braket{B_{\Lambda,exp}-B_{\Lambda,HF}}$",fontsize=14)
cbar.set_label(r"$< (B_{\Lambda,exp}-B_{\Lambda,HF})^2 >^{1/2}$",fontsize=14)
#cbar.set_label("mean deviation squared",fontsize=14)
#cax=plt.axes([0.85,0.1,0.075,0.8])
#plt.colorbar(cax=cax)

#ax.text(7.5,-20,'SLy4',{'color':'k','fontsize':14})
#ax.text(7,-25,r'GKW3(u<1.5)',{'color':'k','fontsize':14})
#ax.text(0.5,-15,'$^{208}_\Lambda$Pb',{'color':'k','fontsize':14})

ax.legend(loc='upper left',frameon=0,numpoints=1,fontsize=14)
ax.set_xlim(-33,-26)
ax.set_ylim(0,600)
#plt.yticks(arange(0.01,0.08,0.02), fontsize=14)
ax.set_xlabel(r'$J_\Lambda$ (MeV)',fontsize=16)
ax.set_ylabel(r'$K_\Lambda$ (MeV)',fontsize=16)
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
plt.savefig("JK_Chi.pdf",dpi=300)
plt.savefig("JK_Chi.png",dpi=300)
plt.show()

#指定可能なファイル形式は emf, eps, jpeg, jpg, pdf, png, ps, raw, rgba, svg,
#svgz, tif, tiff です。拡張子を指定すると勝手に判断されます。
