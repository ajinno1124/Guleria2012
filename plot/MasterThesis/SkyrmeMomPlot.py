#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pylab as pl
from pylab import *
import cmath


d1_mom=pd.read_csv('../../data/Potential/Potential_49/MomentumDep_49.csv') #GKW2+Kohno2
d2_mom=pd.read_csv('../../data/Potential/Potential_50/MomentumDep_50.csv') #GKW3+Kohno3
d3_mom=pd.read_csv('../../data/Potential/Potential_14/MomentumDep_14.csv') #HPL2
d4_mom=pd.read_csv('../../data/Potential/Potential_9/MomentumDep_9.csv') #LY1
d5_mom=pd.read_csv('../../data/Potential/Potential_12/MomentumDep_12.csv') #LY4

fig=plt.figure()
subplots_adjust(hspace=0.0,wspace=0.0,top=0.97,left=0.24,right=0.96)
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


b1=ax.plot(d1_mom["k"],d1_mom["Um"]-d1_mom["Um"][0],label=r"MD1",linewidth=2,color='darkorange',ls='-.')
b2=ax.plot(d2_mom["k"],d2_mom["Um"]-d2_mom["Um"][0],label=r"MD2",linewidth=2,color='darkgreen',ls='--')
b3=ax.plot(d3_mom["k"],d3_mom["Um"]-d3_mom["Um"][0],label=r"HPL2",linewidth=1.5,color='k',ls="-")
b4=ax.plot(d4_mom["k"],d4_mom["Um"]-d4_mom["Um"][0],label=r"LY1",linewidth=2,color='k',ls="--")
b5=ax.plot(d5_mom["k"],d5_mom["Um"]-d5_mom["Um"][0],label=r"LY4",linewidth=2,color='k',ls=":")

#dMS2_mom=pd.read_csv('../../Fitting Parameters/Givendata/Momentum_Dependence/MS2_mom.csv')
#a3=ax.plot(dMS2_mom["k"],dMS2_mom["Um"]-dMS2_mom["Um"][0],label=r"MS2 $\times$ 2/3",linewidth=2.5,color="b",linestyle="-")

#ax.legend(loc='upper left',frameon=0,numpoints=1,fontsize=13)
#ax.legend(loc='upper left',frameon=0,fontsize=16)
#leg1=ax.legend(handles=[b1[0],b2[0],b3[0],b4[0],b5[0]],loc='upper left',frameon=0,numpoints=1,fontsize=16)
#leg1=ax.legend(handles=[b1[0],b2[0],b4[0],b5[0],b6[0]],loc='upper left',frameon=0,numpoints=1,fontsize=16)
leg1=ax.legend(handles=[b1[0],b2[0],b3[0],b4[0],b5[0]],loc='upper left',frameon=0,numpoints=1,fontsize=16)
leg2=ax.legend(handles=[a1,a2],loc='lower right',frameon=0,numpoints=1,fontsize=16)
plt.gca().add_artist(leg1)

ax.set_xlim(0,2)
ax.set_ylim(-5,30)
#plt.yticks(arange(0.01,0.08,0.02), fontsize=14)
#ax.set_ylabel(r'$U_\Lambda(\rho=\rho_0)$ (MeV)',fontsize=16)
ax.set_ylabel(r'$U_{\Lambda}(k)-U_{\Lambda}(k=0)$ (MeV)',fontsize=16)
ax.set_xlabel(r'$k~(\mathrm{fm}^{-1})$',fontsize=16)
ax.tick_params(axis='both', which='both', direction='in',labelsize=16)
ax.tick_params(axis='y',right='true',left='true',labelsize=16)
#plt.text(1.4,0.05,'MS+GKW3',{'color':'k','fontsize':12})
#plt.text(8,0.035,'mid-central Au + Au at 7.7 GeV',{'color':'k','fontsize':12})

plt.tight_layout()

#plt.savefig("v2ch.eps",format='eps',dpi=1000)
#plt.savefig("pxe895qmdmsv.eps",format='eps',bbox__inches='tight')
plt.savefig("SkyrmeMomentumDependence.pdf",dpi=300)
plt.savefig("SkyrmeMomentumDependence.png",dpi=300)
#plt.savefig("timeevolv1.png")
plt.show()

#???????????????????????????????????? emf, eps, jpeg, jpg, pdf, png, ps, raw, rgba, svg,
#svgz, tif, tiff ??????????????????????????????????????????????????????????????????
