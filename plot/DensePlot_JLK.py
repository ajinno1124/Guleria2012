#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pylab as pl
from pylab import *
import cmath

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['figure.subplot.bottom'] = 0.15

#d1_dense=pd.read_csv('../data/Potential/Potential_GKW2/DensityDep_GKW2.csv')
#d2_dense=pd.read_csv('../data/Potential/Potential_GKW3/DensityDep_GKW3.csv')
#d3_dense=pd.read_csv('../data/Potential/Potential_GKW2+MD1/DensityDep_GKW2+MD1.csv')
#d4_dense=pd.read_csv('../data/Potential/Potential_GKW3+MD2/DensityDep_GKW3+MD2.csv')
#d5_dense=pd.read_csv('../data/Potential/Potential_GKW3+MD3/DensityDep_GKW3+MD3.csv')
#d6_dense=pd.read_csv('../data/Potential/Potential_15/DensityDep_15.csv') #HPL2
d7_dense=pd.read_csv('../data/Potential/Potential_10/DensityDep_10.csv') #LY1

#d8_dense=pd.read_csv('../data/Potential/Potential_32/DensityDep_32.csv') #GKW2_medium(rho1.5)
#d9_dense=pd.read_csv('../data/Potential/Potential_33/DensityDep_33.csv') #GKW3_medium(rho1.5)
#d10_dense=pd.read_csv('../data/Potential/Potential_42/DensityDep_42.csv') #GKW2_medium(rho1.5)+Kohno2(k1.5)
#d11_dense=pd.read_csv('../data/Potential/Potential_43/DensityDep_43.csv') #GKW3_medium(rho1.5)+Kohno3(k1.5)


fig=plt.figure()
subplots_adjust(hspace=0.0,wspace=0.0,top=0.9,bottom=0.1,left=0.2,right=0.85)
#subplots_adjust(hspace=1,wspace=0.0,top=1.1,left=0.2,right=0.7)
ax = subplot(1,1,1)

df1=pd.read_csv('../Fitting Parameters/Givendata/SNM2BF_lower.csv')
df2=pd.read_csv('../Fitting Parameters/Givendata/SNM2BF_upper.csv')
df3=pd.read_csv('../Fitting Parameters/Givendata/SNM3BF_lower.csv')
df4=pd.read_csv('../Fitting Parameters/Givendata/SNM3BF_upper.csv')
ax.fill(np.append(df1["density"],df2["density"][::-1]),np.append(df1["U"],df2["U"][::-1]),label="GKW2",color=(0,0,0,0.08),edgecolor='black')
ax.fill(np.append(df3["density"],df4["density"][::-1]),np.append(df3["U"],df4["U"][::-1]),label="GKW3",color=(1,0,0,0.08),edgecolor='black')

#for i in np.arange(51,1563):
for i in np.arange(51,3578):
	df=pd.read_csv(f'../data/Potential/Potential_{i}/DensityDep_{i}.csv')
	ax.plot(df["density"],df["U"],color='k',linestyle='-')	

#ax.plot(d1_dense["density"],d1_dense["U"],label="GKW2 Fit",linewidth=1.9,color='tab:brown',linestyle=':')
#ax.plot(d2_dense["density"],d2_dense["U"],label="GKW3 Fit",linewidth=2.2,color='red',linestyle='-')
#ax.plot(d3_dense["density"],d3_dense["U"],label="GKW2+MD1",linewidth=1.5,color='darkorange',linestyle='-.')
#ax.plot(d4_dense["density"],d4_dense["U"],label="GKW3+MD2",linewidth=2.2,color='m',linestyle='--')
#ax.plot(d5_dense["density"],d5_dense["U"],label="GKW3+MD3",linewidth=2.7,color='darkgreen',linestyle=':')
#ax.plot(d6_dense["density"],d6_dense["U"],label="HPL2",linewidth=2,color='m',linestyle='-.')
ax.plot(d7_dense["density"],d7_dense["U"],label="LY1",linewidth=2,color='k',linestyle=':')

dMS2_dense=pd.read_csv('../Fitting Parameters/Givendata/Density_Dependence/MS2_dense.csv')
ax.plot(dMS2_dense["density"],dMS2_dense["U"],label=r"MS2 $\times$ 2/3",linewidth=2.5,color="b",linestyle="-")

#df=pd.read_csv('SNM_mom_k=0.csv')
#xdata=df[df.density<3.5].density
#ydata=df[df.density<3.5].U
#ax.scatter(df["density"],df["U"],s=1,label='Kohno', alpha=1, linewidths=10, color="black")

#ax1.text(-0.9,0.06,'charged hadrons',{'color':'k','fontsize':16})

ax.legend(loc='best',frameon=0,numpoints=1,fontsize=13)
ax.set_xlim(0,2)
ax.set_ylim(-50,30)
plt.yticks(arange(-40,30.1,10), fontsize=14)
ax.set_ylabel(r'$U_\Lambda$ (MeV)',fontsize=16)
ax.set_xlabel(r'$\rho/\rho_0$',fontsize=16)
#ax3.tick_params(axis='both', which='both', direction='in',labelsize=14)
#plt.text(1.4,0.05,'MS+GKW3',{'color':'k','fontsize':12})
#plt.text(8,0.035,'mid-central Au + Au at 7.7 GeV',{'color':'k','fontsize':12})

#plt.plot((0.0,40),(0,0.0),'k:')

#xticklabels = ax1.get_xticklabels()+ax2.get_xticklabels()
#setp(xticklabels, visible=False)

#yticklabels = ax2.get_yticklabels()+ax4.get_yticklabels()
#setp(yticklabels, visible=False)

#ax.set_xticks([0,0.5,1,1.5,2])
#ax.set_yticks([-50,0,50,100,200,300])
ax.tick_params(axis='x',direction='in')
ax.tick_params(axis='y',which='both',left='true',right='true',direction='in')
ax.tick_params(labelsize=12)

#plt.tight_layout()

#plt.axis([0.0,1.0, 0.0,0.2])
#plot.ylabel(r'$v_2$',fontsize=20)
#plt.text(1.25,9,'Au+Au b<3.4 fm',{'color':'b','fontsize':20})
#plt.text(0.17,11.5,'$\mu=0.0$',{'color':'b','fontsize':20})
#plt.text(0.17,10,r'hadrons=$(n,\Delta,\pi,\rho$)',{'color':'b','fontsize':20})
#plt.xlim(0.0,15)
#plt.ylim(0.0,10)

#plt.savefig("v2ch.eps",format='eps',dpi=1000)
#plt.savefig("pxe895qmdmsv.eps",format='eps',bbox__inches='tight')
plt.savefig("DensityDependence_JLK.pdf",dpi=300)
plt.savefig("DensityDependence_JLK.png",dpi=300)
#plt.savefig("timeevolv1.png")
plt.show()

#指定可能なファイル形式は emf, eps, jpeg, jpg, pdf, png, ps, raw, rgba, svg,
#svgz, tif, tiff です。拡張子を指定すると勝手に判断されます。