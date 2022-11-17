#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pylab as pl
from pylab import *
import cmath

# create BindingEnergyLY.csv
def BindingEnergyL(NParamType,LParamType):
	ZN=np.array([
		[6,6],
		#[6,7],
		#[7,8],# not included in Ohnishi-san's File
        [8,7],
        [14,13],
		[16,15],
		[20,19],
        [23,27],
		[26,29],
        [39,49],
        [57,81],
        [82,125]
	])

	f=open(f'BindingEnergy{NParamType}{LParamType}.csv','w',encoding="utf-8")
	f.write('#Elcheck = e_Lam - (e_Lam using Rearrangement Energy)\n')
	f.write('A,Z,N,jLam,lLam,B.E Lambda(MeV),EL_check(MeV)\n')

	for i in range(len(ZN)):
		df1=pd.read_csv(f"../../data/{NParamType}NaN/Z{ZN[i][0]}N{ZN[i][1]}L0_{NParamType}NaN/Energy2.csv",comment="#")
		df2=pd.read_csv(f"../../data/{NParamType}{LParamType}/Z{ZN[i][0]}N{ZN[i][1]}L1_{NParamType}{LParamType}/Energy.csv",comment="#")
		l=0
		for n in range(len(df2)):
			if df2["lLam"][n]==l:
				f.write(f'{ZN[i][0]+ZN[i][1]}')
				f.write(f',{ZN[i][0]}')
				f.write(f',{ZN[i][1]}')
				f.write(f',{df2["jLam"][n]}')
				f.write(f',{df2["lLam"][n]}')
				f.write(f',{-df2["Etot(MeV)"][n]+df1["Etot2(MeV)"][0]}')
				f.write(f',{df2["El_Check(MeV)"][n]}\n')
				l+=1

	f.close()

#BindingEnergyL("SK3","LY1")
#BindingEnergyL("SK3","LY1_a3zero")
BindingEnergyL("SLy4","LY1")
BindingEnergyL("SK3","LY1")

################################################3
# main plot process

#plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['text.usetex'] = True
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['figure.subplot.bottom'] = 0.15

fig=plt.figure()
subplots_adjust(hspace=0.0,wspace=0.0,top=0.9,left=0.2,right=0.85)
ax = subplot(1,1,1)

df_data=pd.read_csv("../../LamBindingEnergy.csv",comment="#")
df_data=df_data[df_data["Core A"]>=12]
ax.errorbar(df_data["Core A"]**(-2/3),df_data["B. E. (MeV)"],yerr=df_data["error(MeV)"],label="exp.",fmt='o',markersize=5,ecolor='k',markeredgecolor = "black",color='k',zorder=10)

def Plot_OnePot(NParamType,LParamType,color,fmt,linewidth):
	df=pd.read_csv(f'BindingEnergy{NParamType}{LParamType}.csv',comment='#')
	d0=df[df["lLam"]==0]
	d1=df[df["lLam"]==1]
	d2=df[df["lLam"]==2]
	d3=df[df["lLam"]==3]
	d4=df[df["lLam"]==4]
	d5=df[df["lLam"]==5]

	#print(d1_pot.head)

	ax.plot(d0["A"]**(-2/3),d0["B.E Lambda(MeV)"],fmt,label=f"{NParamType}, {LParamType}",c=color,lw=linewidth,ms=5,fillstyle='none')
	ax.plot(d1["A"]**(-2/3),d1["B.E Lambda(MeV)"],fmt,c=color,lw=linewidth,ms=5,fillstyle='none')
	ax.plot(d2["A"]**(-2/3),d2["B.E Lambda(MeV)"],fmt,c=color,lw=linewidth,ms=5,fillstyle='none')
	ax.plot(d3["A"]**(-2/3),d3["B.E Lambda(MeV)"],fmt,c=color,lw=linewidth,ms=5,fillstyle='none')
	ax.plot(d4["A"]**(-2/3),d4["B.E Lambda(MeV)"],fmt,c=color,lw=linewidth,ms=5,fillstyle='none')
	ax.plot(d5["A"]**(-2/3),d5["B.E Lambda(MeV)"],fmt,c=color,lw=linewidth,ms=5,fillstyle='none')

Plot_OnePot("SLy4","LY1",'b','o-',1)
Plot_OnePot("SK3","LY1",'r','s-',1)


#ax.text(0.1,28,r'SK3',{'color':'k','fontsize':14})
#ax.text(0.1,26,r'LY1',{'color':'k','fontsize':14})

ax.legend(loc='upper right',frameon=0,numpoints=1,fontsize=14)
ax.set_xlim(0.0,0.25)
ax.set_ylim(-2,30)
#plt.yticks(arange(0.01,0.08,0.02), fontsize=14)
ax.set_ylabel('$B_\Lambda $ (MeV)',fontsize=16)
ax.set_xlabel(r'$A^{-2/3}$',fontsize=16)
ax.tick_params(axis='x', which='both',top='true',bottom='true', direction='in',labelsize=14)
ax.tick_params(axis='y', which='both',left='true',right='true', direction='in',labelsize=14)
#plt.tight_layout()

#ax.set_xticks([0,1,2,3,4,5])
#ax.set_yticks([-50,0,50,100,200,300])
#ax.tick_params(labelsize=12)

plt.savefig("CompareSLy4SK3.pdf",dpi=300)
plt.savefig("CompareSLy4SK3.png",dpi=300)
plt.show()

#指定可能なファイル形式は emf, eps, jpeg, jpg, pdf, png, ps, raw, rgba, svg,
#svgz, tif, tiff です。拡張子を指定すると勝手に判断されます。
