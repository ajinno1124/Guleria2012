#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pylab as pl
from pylab import *
import cmath

#d1_mom=pd.read_csv('../data/Potential_GKW2/MomentumDep_GKW2.csv')
#d2_mom=pd.read_csv('../data/Potential_GKW3/MomentumDep_GKW3.csv')
#d3_mom=pd.read_csv('../data/Potential_GKW2+MD1/MomentumDep_GKW2+MD1.csv')
#d4_mom=pd.read_csv('../data/Potential_GKW3+MD2/MomentumDep_GKW3+MD2.csv')
#d5_mom=pd.read_csv('../data/Potential_GKW3+MD3/MomentumDep_GKW3+MD3.csv')
#d6_mom=pd.read_csv('../data/Potential_HPL2/MomentumDep_HPL2.csv')
#d7_mom=pd.read_csv('../data/Potential_LY1/MomentumDep_LY1.csv')

#d8_mom=pd.read_csv('../data/Potential_GKW2_1.5/MomentumDep_GKW2_1.5.csv')
#d9_mom=pd.read_csv('../data/Potential_GKW3_1.5/MomentumDep_GKW3_1.5.csv')
#d10_mom=pd.read_csv('../data/Potential_GKW2_1.5+Kohno2/MomentumDep_GKW2_1.5+Kohno2.csv')
#d11_mom=pd.read_csv('../data/Potential_GKW3_1.5+Kohno3/MomentumDep_GKW3_1.5+Kohno3.csv')

d12_mom=pd.read_csv('../data/Potential_GKW3_1.5+Kohno3/MomentumDep_GKW3_1.5+Kohno3.csv')
d13_mom=pd.read_csv('../data/Potential_GKW3_1.5+Kohno3/MomentumDep_GKW3_1.5+Kohno3.csv')

fig=plt.figure()
subplots_adjust(hspace=0.0,wspace=0.0,top=0.97,left=0.24,right=0.96)
ax = subplot(1,1,1)

df1_dat=pd.read_csv('../Fitting Parameters/Givendata/SNM_mom.csv')
xdata=df1_dat[df1_dat.k<2.5].k
ydata=df1_dat[df1_dat.k<2.5].Um
#ax.scatter(xdata,ydata,s=1,label='Kohno3', alpha=1, linewidths=10, color="black")
a1=ax.plot(xdata,ydata-ydata[0],label=r'Kohno3',linewidth=5, color="black")

df2_dat=pd.read_csv('../Fitting Parameters/Givendata/SNM_mom_2BF.csv')
xdata=df2_dat[df2_dat.k<2.5].k
ydata=df2_dat[df2_dat.k<2.5].Um
#ax.scatter(xdata,ydata,s=1,label='Kohno3', alpha=1, linewidths=10, color="black")
a2=ax.plot(xdata,ydata-ydata[0],label=r'Kohno2',linewidth=3.5, color="black",ls='--')


#b1=ax.plot(d3_mom["k"],d3_mom["Um"],label=r"MD1",linewidth=2,color='darkorange',ls='-.')
#b2=ax.plot(d4_mom["k"],d4_mom["Um"],label=r"MD2",linewidth=2,color='m',ls='--')
#b3=ax.plot(d5_mom["k"],d5_mom["Um"],label=r"MD3",linewidth=2,color='darkgreen')
b1=ax.plot(d10_mom["k"],d10_mom["Um"]-d10_mom["Um"][0],label=r"Kohno2",linewidth=2,color='darkorange',ls='-.')
b2=ax.plot(d11_mom["k"],d11_mom["Um"]-d11_mom["Um"][0],label=r"Kohno3",linewidth=2,color='m',ls='--')
b4=ax.plot(d6_mom["k"],d6_mom["Um"]-d6_mom["Um"][0],label=r"HPL2",linewidth=2,color='tab:cyan',ls=":")
b5=ax.plot(d7_mom["k"],d7_mom["Um"]-d7_mom["Um"][0],label=r"LY1",linewidth=2,color='k',ls=":")

#ax.legend(loc='upper left',frameon=0,numpoints=1,fontsize=13)
#ax.legend(loc='upper left',frameon=0,fontsize=16)
#leg1=ax.legend(handles=[b1[0],b2[0],b3[0],b4[0],b5[0]],loc='upper left',frameon=0,numpoints=1,fontsize=16)
leg1=ax.legend(handles=[b1[0],b2[0],b4[0],b5[0]],loc='upper left',frameon=0,numpoints=1,fontsize=16)
leg2=ax.legend(handles=[a1[0],a2[0]],loc='lower right',frameon=0,numpoints=1,fontsize=16)
plt.gca().add_artist(leg1)

ax.set_xlim(0,3)
ax.set_ylim(-3,40)
#plt.yticks(arange(0.01,0.08,0.02), fontsize=14)
#ax.set_ylabel(r'$U_\Lambda(\rho=\rho_0)$ (MeV)',fontsize=16)
ax.set_ylabel(r'$U_{\mathrm{opt}}(\rho=\rho_0)$ (MeV)',fontsize=16)
ax.set_xlabel(r'$p~(\mathrm{fm}^{-1})$',fontsize=16)
ax.tick_params(axis='both', which='both', direction='in',labelsize=16)
ax.tick_params(axis='y',right='true',left='true',labelsize=16)
#plt.text(1.4,0.05,'MS+GKW3',{'color':'k','fontsize':12})
#plt.text(8,0.035,'mid-central Au + Au at 7.7 GeV',{'color':'k','fontsize':12})

plt.tight_layout()

#plt.savefig("v2ch.eps",format='eps',dpi=1000)
#plt.savefig("pxe895qmdmsv.eps",format='eps',bbox__inches='tight')
plt.savefig("MomentumDependence.pdf",dpi=300)
plt.savefig("MomentumDependence.png",dpi=300)
#plt.savefig("timeevolv1.png")
plt.show()

#指定可能なファイル形式は emf, eps, jpeg, jpg, pdf, png, ps, raw, rgba, svg,
#svgz, tif, tiff です。拡張子を指定すると勝手に判断されます。
