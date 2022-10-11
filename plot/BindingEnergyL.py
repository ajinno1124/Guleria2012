#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd

def BindingEnergyL(NParamType,LParamType):
	ZN=[
        [6,5],
        [8,7],
        [14,13],
        [23,27],
        [39,49],
        [57,81],
        [82,125]
    ]

	f=open('BindingEnergyLY.csv','w',encoding="utf-8")
	f.write('#Elcheck = e_Lam - (e_Lam using Rearrangement Energy)\n')
	f.write('A,Z,N,jLam,lLam,B.E Lambda(MeV),EL_check(MeV)\n')

	for i in range(len(AN)):
		NParamType="SK3"
		LParamType="LY1"
		df1=pd.read_csv(f"../data/Z{AN[i][0]}N{AN[i][1]}L0_{NParamType}NaN/Energy2.csv",comment="#")
		df2=pd.read_csv(f"../data/Z{AN[i][0]}N{AN[i][1]}L1_{NParamType}{LParamType}/Energy.csv",comment="#")
		l=0
		for n in range(len(df2)):
			if df2["lLam"][n]==l:
				f.write(f'{AN[i][0]+AN[i][1]}')
				f.write(f',{AN[i][0]}')
				f.write(f',{AN[i][1]}')
				f.write(f',{df2["jLam"][n]}')
				f.write(f',{df2["lLam"][n]}')
				f.write(f',{-df2["Etot(MeV)"][n]+df1["Etot2(MeV)"][0]}')
				f.write(f',{df2["El_Check(MeV)"][n]}\n')
				l+=1

	f.close()

BindingEnergyL("SK3","LY1")
