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

fig=plt.figure()
subplots_adjust(hspace=0.0,wspace=0.0,top=0.9,left=0.2,right=0.85)
ax = subplot(1,1,1)

df_data=pd.read_csv("../../LamBindingEnergy.csv",comment="#",skiprows=[23])
df_data=df_data[df_data["Core A"]>=12]
#df_data=df_data[df_data["lLam"]==0]
ax.errorbar(df_data["Core A"]**(-2/3),df_data["B. E. (MeV)"],yerr=df_data["error(MeV)"],label="exp.",fmt='x',markersize=7,ecolor='k',markeredgecolor = "black",color='k',zorder=10,fillstyle='none')

def Plot_OnePot(NParamType,LParamType,label,color,fmt,linewidth):
	df=pd.read_csv(f'../../data/BindingEnergyLam/BindingEnergy{NParamType}{LParamType}.csv',comment='#')
	d0=df[df["lLam"]==0]
	d1=df[df["lLam"]==1]
	d2=df[df["lLam"]==2]
	d3=df[df["lLam"]==3]
	d4=df[df["lLam"]==4]
	#d5=df[df["lLam"]==5]
	#print(d1_pot.head)

	ax.plot(d0["A"]**(-2/3),d0["BE_Int(MeV)"],fmt,label=label,c=color,lw=linewidth,ms=7,fillstyle='none')
	ax.plot(d1["A"]**(-2/3),d1["BE_Int(MeV)"],fmt,c=color,lw=linewidth,ms=7,fillstyle='none')
	ax.plot(d2["A"]**(-2/3),d2["BE_Int(MeV)"],fmt,c=color,lw=linewidth,ms=7,fillstyle='none')
	ax.plot(d3["A"]**(-2/3),d3["BE_Int(MeV)"],fmt,c=color,lw=linewidth,ms=7,fillstyle='none')
	ax.plot(d4["A"]**(-2/3),d4["BE_Int(MeV)"],fmt,c=color,lw=linewidth,ms=7,fillstyle='none')
	#ax.plot(d5["A"]**(-2/3),d5["BE_Int(MeV)"],fmt,c=color,lw=linewidth,ms=7,fillstyle='none')

#Plot_OnePot("SLy4","GKW2_medium(rho1.5)+Kohno2(k1.5)_a3tuned",'GKW2+Kohno2','b','o--',1.5)
Plot_OnePot("SLy4","GKW3_medium(rho1.5)_a3tuned",'GKW3','r','s-',1.5)
Plot_OnePot("SLy4","GKW3_medium(rho1.5)+Kohno3(k1.5)_a3tuned",'GKW3+Kohno3','darkgreen','o--',1.5)
#Plot_OnePot("SLy4","LY1","LY1",'black','^-',1.5)
Plot_OnePot("SLy4","LY4","LY4",'black','^--',1.5)
Plot_OnePot("SLy4","HPL2","HPL2",'k','D:',1.5)

ax.set_xlim(0,0.25)
ax.set_ylim(-3,30)
ax.legend(loc='upper right',frameon=0,numpoints=1,fontsize=14)
#ax.legend(loc='lower left',frameon=0,numpoints=1,fontsize=14)
#ax.set_xlim(0.20,0.25)
#ax.set_ylim(11.0,13.5)
plt.yticks(arange(0,30.1,5), fontsize=14)
ax.set_ylabel('$B_\Lambda $ (MeV)',fontsize=16)
ax.set_xlabel(r'$A^{-2/3}$',fontsize=16)
ax.tick_params(axis='x', which='both',top='true',bottom='true', direction='in',labelsize=14)
ax.tick_params(axis='y', which='both',left='true',right='true', direction='in',labelsize=14)
#plt.tight_layout()

#ax.set_xticks([0,1,2,3,4,5])
#ax.set_yticks([-50,0,50,100,200,300])
#ax.tick_params(labelsize=12)

#plt.tight_layout()

plt.savefig("BLam_GKW3.pdf",dpi=300)
plt.savefig("BLam_GKW3.png",dpi=300)
plt.show()

#指定可能なファイル形式は emf, eps, jpeg, jpg, pdf, png, ps, raw, rgba, svg,
#svgz, tif, tiff です。拡張子を指定すると勝手に判断されます。
