#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pylab as pl
from pylab import *
import cmath

isMD=True

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['figure.subplot.bottom'] = 0.15

d1_dense=pd.read_csv('../../data/Potential/Potential_21/DensityDep_21.csv')
d2_dense=pd.read_csv('../../data/Potential/Potential_22/DensityDep_22.csv')
d3_dense=pd.read_csv('../../data/Potential/Potential_23/DensityDep_23.csv')
d4_dense=pd.read_csv('../../data/Potential/Potential_24/DensityDep_24.csv')
d5_dense=pd.read_csv('../../data/Potential/Potential_25/DensityDep_25.csv')

fig=plt.figure()
subplots_adjust(hspace=0.0,wspace=0.0,top=0.9,bottom=0.1,left=0.2,right=0.85)
#subplots_adjust(hspace=1,wspace=0.0,top=1.1,left=0.2,right=0.7)
ax = subplot(1,1,1)

df1=pd.read_csv('../../Fitting Parameters/Givendata/SNM2BF_lower.csv')
df2=pd.read_csv('../../Fitting Parameters/Givendata/SNM2BF_upper.csv')
df3=pd.read_csv('../../Fitting Parameters/Givendata/SNM3BF_lower.csv')
df4=pd.read_csv('../../Fitting Parameters/Givendata/SNM3BF_upper.csv')
ax.fill(np.append(df1["density"],df2["density"][::-1]),np.append(df1["U"],df2["U"][::-1]),label="GKW2",color=(0,0,0,0.08),edgecolor='black')
ax.fill(np.append(df3["density"],df4["density"][::-1]),np.append(df3["U"],df4["U"][::-1]),label="GKW3",color=(1,0,0,0.08),edgecolor='black')

ax.plot(d1_dense["density"],d1_dense["U"],label="KIDS0",linewidth=1,ls="-",c='b')
ax.plot(d2_dense["density"],d2_dense["U"],label="KIDSA",linewidth=1,ls="--",c='k')
ax.plot(d3_dense["density"],d3_dense["U"],label="KIDSB",linewidth=1,ls="-.",c='r')
ax.plot(d4_dense["density"],d4_dense["U"],label="KIDSC",linewidth=1,ls="-",c='darkgreen')
ax.plot(d5_dense["density"],d5_dense["U"],label="KIDSD",linewidth=1,ls="-",c='darkorange')

#ax.plot(RhoN,d_KIDS0,label="KIDS0",ls="-",c='b')
#ax.plot(RhoN,d_KIDSA,label="KIDSA",ls="--",c='k')
#ax.plot(RhoN,d_KIDSB,label="KIDSB",ls="-.",c='r')
#ax.plot(RhoN,d_KIDSC,label="KIDSC",ls="-",c='darkgreen')
#ax.plot(RhoN,d_KIDS0,label="KIDSD",ls="-",c='darkorange')

#ax.plot(d1_dense["density"],d1_dense["U"],label="GKW2 Fit",linewidth=1.9,color='tab:brown',linestyle=':')
#ax.plot(d2_dense["density"],d2_dense["U"],label="GKW3 Fit",linewidth=2.2,color='red',linestyle='-')
#ax.plot(d3_dense["density"],d3_dense["U"],label="GKW2+MD1",linewidth=1.5,color='darkorange',linestyle='-.')
#ax.plot(d4_dense["density"],d4_dense["U"],label="GKW3+MD2",linewidth=2.2,color='m',linestyle='--')
#ax.plot(d5_dense["density"],d5_dense["U"],label="GKW3+MD3",linewidth=2.7,color='darkgreen',linestyle=':')
#ax.plot(d6_dense["density"],d6_dense["U"],label="HPL2",linewidth=2,color='m',linestyle='-.')
#ax.plot(d7_dense["density"],d7_dense["U"],label="LY1",linewidth=2,color='k',linestyle=':')

#dMS2_dense=pd.read_csv('../Fitting Parameters/Givendata/Density_Dependence/MS2_dense.csv')
#ax.plot(dMS2_dense["density"],dMS2_dense["U"],label=r"MS2 $\times$ 2/3",linewidth=2.5,color="b",linestyle="-")

#df=pd.read_csv('SNM_mom_k=0.csv')
#xdata=df[df.density<3.5].density
#ydata=df[df.density<3.5].U
#ax.scatter(df["density"],df["U"],s=1,label='Kohno', alpha=1, linewidths=10, color="black")

#ax1.text(-0.9,0.06,'charged hadrons',{'color':'k','fontsize':16})

ax.legend(loc='upper center',frameon=0,numpoints=1,fontsize=13)
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
if isMD==False:
	plt.savefig("DensityDependence.pdf",dpi=300)
	plt.savefig("DensityDependence.png",dpi=300)
elif isMD==True:
	plt.savefig("DensityDependenceMD.pdf",dpi=300)
	plt.savefig("DensityDependenceMD.png",dpi=300)
#plt.savefig("timeevolv1.png")
plt.show()

#???????????????????????????????????? emf, eps, jpeg, jpg, pdf, png, ps, raw, rgba, svg,
#svgz, tif, tiff ??????????????????????????????????????????????????????????????????