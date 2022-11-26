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
#print(df.head(),len(df))


J_mesh=np.array([-33,-32,-31,-30,-29,-28,-27])
L_mesh=np.array([-50,-40,-30,-20,-10,0,10,20])
K_mesh=np.array([0,100,200,300,400,500,600])
ms_m_mesh=np.array([0.6,0.65,0.7,0.75,0.8,0.85,0.90,0.95,1.00])


f=open('Jm_best.csv','w',encoding="utf-8")
f.write("#ChiSquare1 = 1/Nd*sum(B_exp-B_th)^2/sigma\n")
f.write("#ChiSquare2 = sqrt(1/Nd*sum(B_exp-B_th)^2)^2) MeV\n")
f.write("#ChiSquare3 = sum(B_exp-B_th)^2/B_exp^2\n")
f.write("#ChiSquare4 = ChiSquare2 using only s-wave")
f.write("#ChiSquare5 = ChiSquare2 using only s-p splitting\n")
f.write("#ChiSquare6 = ChiSquare2 using only heavyer than 13C_Lam\n")
f.write('index,NParamType,LParameterType,J (MeV),L (MeV),K (MeV),m*/m,ChiSquare6\n')

for j in J_mesh:
	for mi in ms_m_mesh:
		df1=df[abs(df["J (MeV)"]-j)<0.01]
		df1=df1[abs(df1["m*/m"]-mi)<0.01]
		MaxId=df1["ChiSquare6"].idxmin()
		f.write(f'{df1["index"][MaxId]}')
		f.write(f',{df1["NParamType"][MaxId]}')
		f.write(f',{df1["LParameterType"][MaxId]}')
		f.write(f',{df1["J (MeV)"][MaxId]}')
		f.write(f',{df1["L (MeV)"][MaxId]}')
		f.write(f',{df1["K (MeV)"][MaxId]}')
		f.write(f',{df1["m*/m"][MaxId]}')
		f.write(f',{df1["ChiSquare6"][MaxId]}')
		f.write('\n')

f.close()

#d3=d2[d2["ChiSquare2"]<1.5]
#d3=d2.query("ChiSquare2 < 3")
#print(d3["index"])

#J=np.zeros(len(d3))
#L=np.zeros(len(d3))
#K=np.zeros(len(d3))
#ms_m=np.zeros(len(d3))

#for index in d3["index"]:
	#Delete KIDS\
	#if i<=20 and i>=26:
	#if d3["ChiSquare2"][i]<3:
	#J[i]=d1["J (MeV)"][index-1]
	#L[i]=d1["L (MeV)"][index-1]
	#K[i]=d1["K (MeV)"][index-1]
	#ms_m[i]=d1["m*/m"][index-1]
	

#fig=plt.figure()
#subplots_adjust(hspace=0.0,wspace=0.0,top=0.9,left=0.2,right=0.85)
#ax = subplot(1,1,1)
fig, ax=plt.subplots(subplot_kw={"projection":"3d"})

df_best=pd.read_csv('Jm_best.csv',comment='#')
print(df_best.head())

x_data,y_data=np.meshgrid(J_mesh,ms_m_mesh,indexing='ij')
z_data=np.zeros([len(J_mesh),len(ms_m_mesh)])
for i in range(0,len(J_mesh)):
	for j in range(0,len(ms_m_mesh)):
		df_temp=df_best[abs(df_best["J (MeV)"]-J_mesh[i])<0.01]
		df_temp=df_temp[abs(df_best["m*/m"]-ms_m_mesh[j])<0.01]
		print(df_temp)
		if len(df_temp)==0:
			z_data[i][j]=np.nan
		else:
			z_data[i][j]=df_temp["ChiSquare6"][df_temp["ChiSquare6"].idxmin()]


#ax.plot_surface(x_data,y_data,z_data,cmap=cm.jet)
ax.plot_wireframe(x_data,y_data,z_data,color="k")
levels=[0.8,1.0,1.5]
colors=['r','b','darkgreen']
linestyles=['solid','dashed','dashdot']
CS=ax.contour(x_data,y_data,z_data,levels=levels,linestyles=linestyles, colors=colors,labels=levels)

#plt.scatter(df_best['J (MeV)'],df_best["m*/m"],c=df_best["ChiSquare6"],cmap=plt.cm.jet,s=30,edgecolor='k')
#plt.scatter(d1["J (MeV)"][d3["index"]-1],d1["m*/m"][d3["index"]-1],c=d3["ChiSquare5"],cmap=plt.cm.jet,s=30,edgecolor='k')
#print("ChiSquare5")
#cbar=plt.colorbar(aspect=40,pad=0.08,orientation='vertical')
#cbar.set_label(r"$\Braket{B_{\Lambda,exp}-B_{\Lambda,HF}}$".format('B'))
#cbar.set_label(r"$< (B_{\Lambda,exp}-B_{\Lambda,HF})^2 >^{1/2}$",fontsize=14)
#cbar.set_label("mean deviation squared")
#cax=plt.axes([0.85,0.1,0.075,0.8])
#plt.colorbar(cax=cax)

#ax.text(7.5,-20,'SLy4',{'color':'k','fontsize':14})
#ax.text(7,-25,r'GKW3(u<1.5)',{'color':'k','fontsize':14})
#ax.text(0.5,-15,'$^{208}_\Lambda$Pb',{'color':'k','fontsize':14})

plt.clabel(CS, inline=1, fontsize=10)
#plt.title('Simplest default with labels')

labels = [f'{levels[0]}',f'{levels[1]}',f'{levels[2]}']
for i in range(len(labels)):
    CS.collections[i].set_label(labels[i])

ax.legend(loc='upper right',frameon=0,numpoints=1,fontsize=14)
ax.set_xlim(-34,-26)
ax.set_ylim(0.55,1.05)
ax.set_zlim(0,2)
#plt.yticks(arange(0.01,0.08,0.02), fontsize=14)
ax.set_xlabel(r'$J_\Lambda$ (MeV)',fontsize=16)
ax.set_ylabel(r'$m_\Lambda^*/m_\Lambda$',fontsize=16)
#ax.set_zlabel(r'$\ev{(B_{\Lambda, {\rm exp}} - B_{\Lambda {\rm exp})^2}}^{1/2}$',fontsize=16)
ax.tick_params(axis='x', which='both', direction='in',labelsize=14)
ax.tick_params(axis='y', which='both', direction='in',labelsize=14)
ax.tick_params(axis='z', which='both', direction='in',labelsize=14)
plt.tight_layout()
ax.view_init(azim=-126,elev=34)
#ax.set_xticks([0,1,2,3,4,5])
#ax.set_yticks([-50,0,50,100,200,300])
#ax.tick_params(labelsize=12)

#plt.axis([0.0,1.0, 0.0,0.2])
#plot.ylabel(r'$v_2$',fontsize=20)
#plt.text(1.25,9,'Au+Au b<3.4 fm',{'color':'b','fontsize':20})
#plt.text(0.17,11.5,'$\mu=0.0$',{'color':'b','fontsize':20})
#plt.text(0.17,10,r'hadrons=$(n,\Delta,\pi,\rho$)',{'color':'b','fontsize':20})

#plt.savefig("v2ch.eps",format='eps',dpi=1000)
#plt.savefig("pxe895qmdmsv.eps",format='eps',bbox__inches='tight')
plt.savefig("Jm_best.pdf",dpi=300)
plt.savefig("Jm_best.png",dpi=300)
#plt.show()

#指定可能なファイル形式は emf, eps, jpeg, jpg, pdf, png, ps, raw, rgba, svg,
#svgz, tif, tiff です。拡張子を指定すると勝手に判断されます。
