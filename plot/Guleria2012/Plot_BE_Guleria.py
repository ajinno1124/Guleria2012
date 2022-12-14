#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pylab as pl
from pylab import *
import cmath

# create BindingEnergyLY.csv
'''
def BindingEnergyL(NParamType,LParamType):
	ZN=[
		[2,5],
		[4,3],
		[3,5],
		[4,5],
		[5,4],
		[5,5],
		[5,6],
        [6,5],
		[6,6],
		[7,8],
        [8,7],
        [14,13],
		[16,15],
		[20,19],
        [23,27],
        [39,49],
        [57,81],
        [82,125]
    ]

	f=open('BindingEnergyHPL2.csv','w',encoding="utf-8")
	f.write('#Elcheck = e_Lam - (e_Lam using Rearrangement Energy)\n')
	f.write('A,Z,N,jLam,lLam,B.E Lambda(MeV),EL_check(MeV)\n')

	for i in range(len(ZN)):
		df1=pd.read_csv(f"../data/{NParamType}{LParamType}/data/Z{ZN[i][0]}N{ZN[i][1]}L0_{NParamType}NaN/Energy2.csv",comment="#")
		df2=pd.read_csv(f"../data//{NParamType}{LParamType}/data/Z{ZN[i][0]}N{ZN[i][1]}L1_{NParamType}{LParamType}/Energy.csv",comment="#")
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

BindingEnergyL("SLy4","HPL2")
'''

################################################3
# main plot process

#plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['text.usetex'] = True
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['figure.subplot.bottom'] = 0.15

df=pd.read_csv('../BindingEnergyHPL2.csv',comment='#')
d0=df[df["lLam"]==0]
d1=df[df["lLam"]==1]
d2=df[df["lLam"]==2]
d3=df[df["lLam"]==3]
d4=df[df["lLam"]==4]
d5=df[df["lLam"]==5]
#print(d1_pot.head)


fig=plt.figure()
subplots_adjust(hspace=0.0,wspace=0.0,top=0.9,left=0.2,right=0.85)
ax = subplot(1,1,1)

#ax.errorbar(df["A^(-2/3)"],df["exp(MeV)"],yerr=df["exp_error(MeV)"],label="exp.",fmt='o',markersize=5,ecolor='k',markeredgecolor = "black",color='k',zorder=10)
ax.plot(d0["A"]**(-2/3),d0["B.E Lambda(MeV)"],label="$s$ orbital")
ax.plot(d1["A"]**(-2/3),d1["B.E Lambda(MeV)"],label="$p$ orbital")
ax.plot(d2["A"]**(-2/3),d2["B.E Lambda(MeV)"],label="$d$ orbital")
ax.plot(d3["A"]**(-2/3),d3["B.E Lambda(MeV)"],label="$f$ orbital")
ax.plot(d4["A"]**(-2/3),d4["B.E Lambda(MeV)"],label="$g$ orbital")
ax.plot(d5["A"]**(-2/3),d5["B.E Lambda(MeV)"],label="$h$ orbital")

ax.text(0.1,28,r'SLy4',{'color':'k','fontsize':14})
ax.text(0.1,26,r'HPL2',{'color':'k','fontsize':14})

ax.legend(loc='upper right',frameon=0,numpoints=1,fontsize=14)
ax.set_xlim(0.0,0.3)
ax.set_ylim(0,30)
#plt.yticks(arange(0.01,0.08,0.02), fontsize=14)
ax.set_ylabel('$B_\Lambda $ (MeV)',fontsize=16)
ax.set_xlabel(r'$A^{-2/3}$',fontsize=16)
ax.tick_params(axis='x', which='both',top='true',bottom='true', direction='in',labelsize=14)
ax.tick_params(axis='y', which='both',left='true',right='true', direction='in',labelsize=14)
#plt.tight_layout()

#ax.set_xticks([0,1,2,3,4,5])
#ax.set_yticks([-50,0,50,100,200,300])
#ax.tick_params(labelsize=12)

plt.savefig("BElam_Guleria.pdf",dpi=300)
plt.savefig("BElam_Guleria.png",dpi=300)
plt.show()

#???????????????????????????????????? emf, eps, jpeg, jpg, pdf, png, ps, raw, rgba, svg,
#svgz, tif, tiff ??????????????????????????????????????????????????????????????????
