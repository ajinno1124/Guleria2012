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

fig=plt.figure(figsize=(12,12))
subplots_adjust(hspace=0.0,wspace=0.0,top=0.9,left=0.2,right=0.85)

#df_data=pd.read_csv("../../LamBindingEnergy.csv",comment="#",skiprows=[23])
df_data=pd.read_csv("../../LamBindingEnergy.csv",comment="#")
df_data=df_data[df_data["Core A"]>=12]

def Plot_OnePot(ax,NParamType,LParamType,label,color,fmt,linewidth):
	#ax.errorbar(df_data["Core A"]**(-2/3),df_data["B. E. (MeV)"],yerr=df_data["error(MeV)"],label="exp.",fmt='x',markersize=7,ecolor='k',markeredgecolor = "black",color='k',zorder=10,fillstyle='none')
	df=pd.read_csv(f'../../data/BindingEnergyLam/BindingEnergy{NParamType}{LParamType}.csv',comment='#',skiprows=[3,4])
	d0=df[df["lLam"]==0]
	d1=df[df["lLam"]==1]
	d2=df[df["lLam"]==2]
	d3=df[df["lLam"]==3]
	d4=df[df["lLam"]==4]
	#d5=df[df["lLam"]==5]
	#print(d1_pot.head)

	ax.plot((d0["A"]+1)**(-2/3),d0["BE_Int(MeV)"],fmt,label=label,c=color,lw=linewidth,ms=7,fillstyle='none')
	ax.plot((d1["A"]+1)**(-2/3),d1["BE_Int(MeV)"],fmt,c=color,lw=linewidth,ms=7,fillstyle='none')
	ax.plot((d2["A"]+1)**(-2/3),d2["BE_Int(MeV)"],fmt,c=color,lw=linewidth,ms=7,fillstyle='none')
	ax.plot((d3["A"]+1)**(-2/3),d3["BE_Int(MeV)"],fmt,c=color,lw=linewidth,ms=7,fillstyle='none')
	ax.plot((d4["A"]+1)**(-2/3),d4["BE_Int(MeV)"],fmt,c=color,lw=linewidth,ms=7,fillstyle='none')
	#ax.plot(d5["A"]**(-2/3),d5["BE_Int(MeV)"],fmt,c=color,lw=linewidth,ms=7,fillstyle='none')

	ax.set_xlim(-0.025,0.25)
	ax.set_ylim(-10,30)
	#ax.set_xlim(0.20,0.25)
	#ax.set_ylim(11.0,13.5)
	ax.tick_params(axis='x', which='both',top='true',bottom='true', direction='in',labelsize=14,pad=10)
	ax.tick_params(axis='y', which='both',left='true',right='true', direction='in',labelsize=14,pad=10)
	#plt.tight_layout()

def Add_Nuclide(ax):
	ax.text(208**(-2/3)-0.04,-3,r"$^{208} _{\Lambda}{\rm Pb}$",fontsize=14)
	ax.text(139**(-2/3)-0.01,-9,r"$^{139} _{\Lambda}{\rm La}$",fontsize=14)
	ax.text(89**(-2/3)-0.008,-5,r"$^{89} _{\Lambda}{\rm Y}$",fontsize=14)
	ax.text(51**(-2/3)-0.003,-5,r"$^{51} _{\Lambda}{\rm V}$",fontsize=14)
	#ax.text(40**(-2/3)-0.01,-5,r"$^{40} _{\Lambda}{\rm Ca}$",fontsize=14)
	ax.text(28**(-2/3)-0.01,-5,r"$^{28} _{\Lambda}{\rm Si}$",fontsize=14)
	ax.text(16**(-2/3)-0.01,-5,r"$^{16} _{\Lambda}{\rm O}$",fontsize=14)
	ax.text(13**(-2/3)-0.01,-5,r"$^{13} _{\Lambda}{\rm C}$",fontsize=14)
	pointPb = {
        'start': [208**(-2/3), -1],
        'end': [208**(-2/3), 4]
    }
	pointLa = {
        'start': [139**(-2/3), -7],
        'end': [139**(-2/3), 0]
    }
	pointY = {
        'start': [89**(-2/3), -3],
        'end': [89**(-2/3), 0]
    }
	pointV = {
        'start': [51**(-2/3), -3],
        'end': [51**(-2/3), -0.1]
    }
	pointO = {
        'start': [16**(-2/3), -3],
        'end': [16**(-2/3), 1]
    }
	ax.annotate('', xy=pointPb['end'], xytext=pointPb['start'],
                arrowprops=dict(shrink=0, width=1, headwidth=8, 
                                headlength=10, connectionstyle='arc3',
                                facecolor='k', edgecolor='k'))
	ax.annotate('', xy=pointLa['end'], xytext=pointLa['start'],
                arrowprops=dict(shrink=0, width=1, headwidth=8, 
                                headlength=10, connectionstyle='arc3',
                                facecolor='k', edgecolor='k'))
	ax.annotate('', xy=pointY['end'], xytext=pointY['start'],
                arrowprops=dict(shrink=0, width=1, headwidth=8, 
                                headlength=10, connectionstyle='arc3',
                                facecolor='k', edgecolor='k'))
	ax.annotate('', xy=pointV['end'], xytext=pointV['start'],
                arrowprops=dict(shrink=0, width=1, headwidth=8, 
                                headlength=10, connectionstyle='arc3',
                                facecolor='k', edgecolor='k'))
	ax.annotate('', xy=pointO['end'], xytext=pointO['start'],
				arrowprops=dict(shrink=0, width=1, headwidth=8, 
					headlength=10, connectionstyle='arc3',
					facecolor='k', edgecolor='k'))


def Add_Orbital(ax):
	ax.text(-0.005,27,r"$1s$",fontsize=16)
	ax.text(-0.005,22,r"$1p$",fontsize=16)
	ax.text(-0.005,16.5,r"$1d$",fontsize=16)
	ax.text(-0.005,11.5,r"$1f$",fontsize=16)
	ax.text(-0.005,6.5,r"$1g$",fontsize=16)

ax1 = fig.add_subplot(2,2,1)
ax1.errorbar((df_data["Core A"]+1)**(-2/3),df_data["B. E. (MeV)"],yerr=df_data["error(MeV)"],label="exp.",fmt='x',markersize=7,ecolor='k',markeredgecolor = "black",color='k',zorder=10,fillstyle='none')
ax1.legend(loc='upper right',frameon=0,numpoints=1,fontsize=14)
Plot_OnePot(ax1,"SLy4","GKW3_medium(rho1.5)_a3tuned",'GKW3','r','-',1.5)
#Plot_OnePot("SLy4","GKW2_medium(rho1.5)+Kohno2(k1.5)_a3tuned",'GKW2+Kohno2','b','o--',1.5)
ax1.set_xticks(arange(0.00,0.21,0.05), fontsize=14)
ax1.set_yticks(arange(-10,30.1,5), fontsize=14)
ax1.set_ylabel('$B_\Lambda $ (MeV)',fontsize=16)
ax1.text(0.15,20,'GKW3',fontsize=14)
setp(ax1.get_xticklabels(), visible=False)
Add_Nuclide(ax1)
Add_Orbital(ax1)

ax2 = fig.add_subplot(2,2,2)
Plot_OnePot(ax2,"SLy4","GKW3_medium(rho1.5)+Kohno3(k1.5)_a3tuned",'GKW3+Kohno3','darkgreen','-',2)
#ax2.legend(loc='upper right',frameon=0,numpoints=1,fontsize=14)
ax2.errorbar((df_data["Core A"]+1)**(-2/3),df_data["B. E. (MeV)"],yerr=df_data["error(MeV)"],label="exp.",fmt='x',markersize=7,ecolor='k',markeredgecolor = "black",color='k',zorder=10,fillstyle='none')
ax2.set_xticks(arange(0.00,0.21,0.05), fontsize=0)
ax2.set_yticks(arange(-10,29,5), fontsize=0)
setp(ax2.get_xticklabels(), visible=False)
setp(ax2.get_yticklabels(), visible=False)
ax2.text(0.13,20,'GKW3+Kohno3',fontsize=14)
Add_Nuclide(ax2)
Add_Orbital(ax2)

ax3 = fig.add_subplot(2,2,3)
Plot_OnePot(ax3,"SLy4","LY4","LY4",'black','-',1.5)
#ax3.legend(loc='upper right',frameon=0,numpoints=1,fontsize=14)
ax3.errorbar((df_data["Core A"]+1)**(-2/3),df_data["B. E. (MeV)"],yerr=df_data["error(MeV)"],label="exp.",fmt='x',markersize=7,ecolor='k',markeredgecolor = "black",color='k',zorder=10,fillstyle='none')
ax3.set_xticks(arange(0.00,0.21,0.05), fontsize=14)
ax3.set_yticks(arange(-10,29,5), fontsize=14)
ax3.set_xlabel('$A^{-2/3}$',fontsize=16)
ax3.set_ylabel('$B_\Lambda $ (MeV)',fontsize=16)
ax3.text(0.15,20,r'LY${\rm~IV}$',fontsize=14)
Add_Nuclide(ax3)
Add_Orbital(ax3)

ax4 = fig.add_subplot(2,2,4)
Plot_OnePot(ax4,"SLy4","HPL2","HPL2",'k','-',1.5)
#ax4.legend(loc='upper right',frameon=0,numpoints=1,fontsize=14)
ax4.errorbar((df_data["Core A"]+1)**(-2/3),df_data["B. E. (MeV)"],yerr=df_data["error(MeV)"],label="exp.",fmt='x',markersize=7,ecolor='k',markeredgecolor = "black",color='k',zorder=10,fillstyle='none')
ax4.set_xticks(arange(0.00,0.26,0.05), fontsize=14)
ax4.set_yticks(arange(-10,29,5), fontsize=14)
ax4.set_xlabel('$A^{-2/3}$',fontsize=16)
setp(ax4.get_yticklabels(), visible=False)
ax4.text(0.15,20,r'HP$\Lambda$2',fontsize=14)
Add_Nuclide(ax4)
Add_Orbital(ax4)

plt.savefig("BLam_GKW3.pdf",dpi=300)
plt.savefig("BLam_GKW3.png",dpi=300)
plt.show()

#???????????????????????????????????? emf, eps, jpeg, jpg, pdf, png, ps, raw, rgba, svg,
#svgz, tif, tiff ??????????????????????????????????????????????????????????????????
