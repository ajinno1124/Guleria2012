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

df1_dat=pd.read_csv('../../Fitting Parameters/Givendata/SNM_mom.csv')
xdata=df1_dat[df1_dat.k<2.5].k
ydata=df1_dat[df1_dat.k<2.5].Um
a1=ax.scatter(xdata,ydata-ydata[0],marker='o',s=20,label='Kohno3', edgecolor="k",facecolor='none')
#a1=ax.plot(xdata,ydata-ydata[0],label=r'Kohno3',linewidth=5, color="black")

df2_dat=pd.read_csv('../../Fitting Parameters/Givendata/SNM_mom_2BF.csv')
xdata=df2_dat[df2_dat.k<2.5].k
ydata=df2_dat[df2_dat.k<2.5].Um
a2=ax.scatter(xdata,ydata-ydata[0],marker='s',s=20,label='Kohno2', color="black")
#a2=ax.plot(xdata,ydata-ydata[0],label=r'Kohno2',linewidth=3.5, color="black",ls='--')

c300='r'
lw300=3
zo300=1
ls300='-'

c600="b"
lw600=1
zo600=10
ls600='--'

cother="darkgreen"
lwother=1
zoother=15
lsother='--'


for index in df["index"]:
	#print(index,df["K (MeV)"][index-1])
	df_pot=pd.read_csv(f'../../data/Potential/Potential_{index}/MomentumDep_{index}.csv',comment='#')
	if abs(df["K (MeV)"][index-1]-300)<0.1:
		ax.plot(df_pot["k"],df_pot["Um"]-df_pot["Um"][0],color=c300,linewidth=lw300,zorder=zo300,linestyle=ls300)
	elif abs(df["K (MeV)"][index-1]-600)<0.1:
		ax.plot(df_pot["k"],df_pot["Um"]-df_pot["Um"][0],color=c600,linewidth=lw600,zorder=zo600,linestyle=ls600)
	else:
		aother=ax.plot(df_pot["k"],df_pot["Um"]-df_pot["Um"][0],color=cother,linewidth=lwother,zorder=zoother,linestyle=lsother)

nanvec=np.zeros(2)
nanvec[:]=np.nan
b1=ax.plot(nanvec,nanvec,label=r"$K_\Lambda = 300$ (MeV)",color=c300,linewidth=lw300,zorder=zo300,linestyle=ls300)
b2=ax.plot(nanvec,nanvec,label=r"$K_\Lambda = 600$ (MeV)",color=c600,linewidth=lw600,zorder=zo600,linestyle=ls600)
b3=ax.plot(nanvec,nanvec,label=r"others",color=cother,linewidth=lwother,zorder=zoother,linestyle=lsother)

#leg1=ax.legend(handles=[b1[0],b2[0],b3[0]],loc='lower right',frameon=0,numpoints=1,fontsize=16)
leg2=ax.legend(handles=[a1,a2,b1[0],b2[0],b3[0]],loc='upper left',frameon=0,numpoints=1,fontsize=16)
#plt.gca().add_artist(leg1)
ax.set_xlim(0,2)
ax.set_ylim(-5,30)
#plt.yticks(arange(0.01,0.08,0.02), fontsize=14)
#ax.set_ylabel(r'$U_\Lambda(\rho=\rho_0)$ (MeV)',fontsize=16)
ax.set_ylabel(r'$U_{\Lambda}(k)$ (MeV)',fontsize=16)
ax.set_xlabel(r'$k~(\mathrm{fm}^{-1})$',fontsize=16)
ax.tick_params(axis='both', which='both', direction='in',labelsize=16)
ax.tick_params(axis='y',right='true',left='true',labelsize=16)
#plt.text(1.4,0.05,'MS+GKW3',{'color':'k','fontsize':12})
#plt.text(8,0.035,'mid-central Au + Au at 7.7 GeV',{'color':'k','fontsize':12})
ax.tick_params(labelsize=12)

plt.tight_layout()
plt.savefig("MomPot08MeV.pdf",dpi=300)
plt.savefig("MomPot08MeV.png",dpi=300)
plt.show()

#指定可能なファイル形式は emf, eps, jpeg, jpg, pdf, png, ps, raw, rgba, svg,
#svgz, tif, tiff です。拡張子を指定すると勝手に判断されます。
