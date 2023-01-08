#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pylab as pl
from pylab import *
import cmath

#d1_mom=pd.read_csv('../../data/Potential/Potential_GKW2/MomentumDep_GKW2.csv')
#d2_mom=pd.read_csv('../../data/Potential/Potential_GKW3/MomentumDep_GKW3.csv')
#d3_mom=pd.read_csv('../../data/Potential/Potential_GKW2+MD1/MomentumDep_GKW2+MD1.csv')
#d4_mom=pd.read_csv('../../data/Potential/Potential_GKW3+MD2/MomentumDep_GKW3+MD2.csv')
#d5_mom=pd.read_csv('../../data/Potential/Potential_GKW3+MD3/MomentumDep_GKW3+MD3.csv')
#d6_mom=pd.read_csv('../../data/Potential/Potential_HPL2/MomentumDep_HPL2.csv')
#d7_mom=pd.read_csv('../../data/Potential/Potential_9/MomentumDep_9.csv')

#d8_mom=pd.read_csv('../../data/Potential/Potential_GKW2_1.5/MomentumDep_GKW2_1.5.csv')
#d9_mom=pd.read_csv('../../data/Potential/Potential_GKW3_1.5/MomentumDep_GKW3_1.5.csv')
#d10_mom=pd.read_csv('../../data/Potential/Potential_GKW2_1.5+Kohno2/MomentumDep_GKW2_1.5+Kohno2.csv')
#d11_mom=pd.read_csv('../../data/Potential/Potential_GKW3_1.5+Kohno3/MomentumDep_GKW3_1.5+Kohno3.csv')

#d12_mom=pd.read_csv('../../data/Potential/Potential_42/MomentumDep_42.csv') #GKW2_medium(rho1.5)+Kohno2(k1.5)
#d13_mom=pd.read_csv('../../data/Potential/Potential_43/MomentumDep_43.csv') #GKW3_medium(rho1.5)+Kohno3(k1.5)

fig=plt.figure()
subplots_adjust(hspace=0.0,wspace=0.0,top=0.97,left=0.24,right=0.96)
ax = subplot(1,1,1)

df1_dat=pd.read_csv('../../Fitting Parameters/Givendata/SNM_mom.csv')
xdata=df1_dat[df1_dat.k<2.5].k
ydata=df1_dat[df1_dat.k<2.5].Um
a1=ax.scatter(xdata,ydata-ydata[0],marker='o',s=15,label='Kohno3', edgecolor="black",facecolor='None')
#a1=ax.plot(xdata,ydata-ydata[0],label=r'Kohno3',linewidth=5, color="black")

df2_dat=pd.read_csv('../../Fitting Parameters/Givendata/SNM_mom_2BF.csv')
xdata=df2_dat[df2_dat.k<2.5].k
ydata=df2_dat[df2_dat.k<2.5].Um
a2=ax.scatter(xdata,ydata-ydata[0],marker='s',s=15,label='Kohno2', color="black")
#a2=ax.plot(xdata,ydata-ydata[0],label=r'Kohno2',linewidth=3.5, color="black",ls='--')

#for i in np.arange(51,1563):
for i in np.arange(51,3578):
    df=pd.read_csv(f'../../data/Potential/Potential_{i}/MomentumDep_{i}.csv')
    ax.plot(df["k"],df["Um"]-df["Um"][0],color=(0,0,1,0.05),linestyle='-',linewidth=0.5,rasterized=True)


#b1=ax.plot(d3_mom["k"],d3_mom["Um"],label=r"MD1",linewidth=2,color='darkorange',ls='-.')
#b2=ax.plot(d4_mom["k"],d4_mom["Um"],label=r"MD2",linewidth=2,color='m',ls='--')
#b3=ax.plot(d5_mom["k"],d5_mom["Um"],label=r"MD3",linewidth=2,color='darkgreen')
#b1=ax.plot(d12_mom["k"],d12_mom["Um"]-d12_mom["Um"][0],label=r"Kohno2 ($k<1.5$/fm) Fit",linewidth=1,color='darkgreen',ls='-')
#b2=ax.plot(d13_mom["k"],d13_mom["Um"]-d13_mom["Um"][0],label=r"Kohno3 ($k<1.5$/fm) Fit",linewidth=1.9,color='r',ls='-')
#b4=ax.plot(d6_mom["k"],d6_mom["Um"]-d6_mom["Um"][0],label=r"HPL2",linewidth=2,color='tab:cyan',ls=":")
#b5=ax.plot(d7_mom["k"],d7_mom["Um"]-d7_mom["Um"][0],label=r"LY1",linewidth=2,color='k',ls=":")

dMS2_mom=pd.read_csv('../../Fitting Parameters/Givendata/Momentum_Dependence/MS2_mom.csv')
a3=ax.plot(dMS2_mom["k"],dMS2_mom["Um"]-dMS2_mom["Um"][0],label=r"MS2 $\times$ 2/3",linewidth=2.5,color="b",linestyle="-")

#ax.legend(loc='upper left',frameon=0,numpoints=1,fontsize=13)
#ax.legend(loc='upper left',frameon=0,fontsize=16)
#leg1=ax.legend(handles=[b1[0],b2[0],b3[0],b4[0],b5[0]],loc='upper left',frameon=0,numpoints=1,fontsize=16)
#leg1=ax.legend(handles=[b1[0],b2[0],b5[0]],loc='upper left',frameon=0,numpoints=1,fontsize=16)
leg2=ax.legend(handles=[a1,a2,a3[0]],loc='upper left',frameon=0,numpoints=1,fontsize=16)

ax.set_xlim(0,3)
ax.set_ylim(-3,40)
#plt.yticks(arange(0.01,0.08,0.02), fontsize=14)
#ax.set_ylabel(r'$U_\Lambda(\rho=\rho_0)$ (MeV)',fontsize=16)
ax.set_ylabel(r'$U_{\mathrm{opt}}(k)-U_{\mathrm{opt}}(k=0)$ (MeV)',fontsize=16)
ax.set_xlabel(r'$k~(\mathrm{fm}^{-1})$',fontsize=16)
ax.tick_params(axis='both', which='both', direction='in',labelsize=16)
ax.tick_params(axis='y',right='true',left='true',labelsize=16)
#plt.text(1.4,0.05,'MS+GKW3',{'color':'k','fontsize':12})
#plt.text(8,0.035,'mid-central Au + Au at 7.7 GeV',{'color':'k','fontsize':12})

plt.tight_layout()

#plt.savefig("v2ch.eps",format='eps',dpi=1000)
#plt.savefig("pxe895qmdmsv.eps",format='eps',bbox__inches='tight')
plt.savefig("MomentumDependence_JLK.pdf",dpi=300)
plt.savefig("MomentumDependence_JLK.png",dpi=300)
#plt.savefig("timeevolv1.png")
plt.show()

#指定可能なファイル形式は emf, eps, jpeg, jpg, pdf, png, ps, raw, rgba, svg,
#svgz, tif, tiff です。拡張子を指定すると勝手に判断されます。
