#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pylab as pl
from pylab import *
from matplotlib import cm
import cmath


#plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.family'] = 'Times New Roman'
#plt.rcParams['text.usetex'] = True
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['figure.subplot.bottom'] = 0.15
#plt.rc('text.latex', preamble=r'\usepackage{braket}')
#plt.rc('text.latex', preamble=r'\usepackage{physics}')


#d1=pd.read_csv('../JLK.csv',comment='#')
d1=pd.read_csv('../../data/BindingEnergyLam/ChiSquared.csv',comment='#')
df=d1[d1["Number of Data"]==25.0]
df=df[df["index"]>50]
df=df[df["ChiSquare6"]<0.8]

fig=plt.figure()
subplots_adjust(hspace=0.0,wspace=0.0,top=0.9,bottom=0.1,left=0.2,right=0.85)
#subplots_adjust(hspace=1,wspace=0.0,top=1.1,left=0.2,right=0.7)
ax = subplot(1,1,1)

df1=pd.read_csv('../../Fitting Parameters/Givendata/SNM2BF_lower.csv')
df2=pd.read_csv('../../Fitting Parameters/Givendata/SNM2BF_upper.csv')
df3=pd.read_csv('../../Fitting Parameters/Givendata/SNM3BF_lower.csv')
df4=pd.read_csv('../../Fitting Parameters/Givendata/SNM3BF_upper.csv')
ax.fill(np.append(df1["density"],df2["density"][::-1]),np.append(df1["U"],df2["U"][::-1]),label="GKW2",color=(0,0,1,0.08),edgecolor='black')
ax.fill(np.append(df3["density"],df4["density"][::-1]),np.append(df3["U"],df4["U"][::-1]),label="GKW3",color=(1,0,0,0.08),edgecolor='black')

for index in df["index"]:
	#print(index,df["K (MeV)"][index-1])
	df_pot=pd.read_csv(f'../../data/Potential/Potential_{index}/DensityDep_{index}.csv',comment='#')
	if abs(df["K (MeV)"][index-1]-300)<0.1:
		ax.plot(df_pot["density"],df_pot["U"],color="r",linewidth=1)
	elif abs(df["K (MeV)"][index-1]-600)<0.1:
		ax.plot(df_pot["density"],df_pot["U"],color="b",linewidth=1)
	else:
		aother=ax.plot(df_pot["density"],df_pot["U"],color="darkgreen",linewidth=1)

nanvec=np.zeros(2)
nanvec[:]=np.nan
ax.plot(nanvec,nanvec,color="r",linewidth=1,label=r"$K_\Lambda = 300$ (MeV)")
ax.plot(nanvec,nanvec,color="b",linewidth=1,label=r"$K_\Lambda = 300$ (MeV)")
ax.plot(nanvec,nanvec,color="darkgreen",linewidth=1,label=r"others")

ax.legend(loc='best',frameon=0,numpoints=1,fontsize=13)
ax.set_xlim(0,2)
ax.set_ylim(-50,30)
plt.xticks(arange(0,2.1,0.5), fontsize=14)
plt.yticks(arange(-40,30.1,10), fontsize=14)
ax.set_ylabel(r'$U_\Lambda$ (MeV)',fontsize=16)
ax.set_xlabel(r'$\rho/\rho_0$',fontsize=16)

ax.tick_params(axis='x',direction='in')
ax.tick_params(axis='y',which='both',left='true',right='true',direction='in')
ax.tick_params(labelsize=12)

plt.tight_layout()
plt.savefig("DensePot08MeV.pdf",dpi=300)
plt.savefig("DensePot08MeV.png",dpi=300)
plt.show()

#指定可能なファイル形式は emf, eps, jpeg, jpg, pdf, png, ps, raw, rgba, svg,
#svgz, tif, tiff です。拡張子を指定すると勝手に判断されます。
