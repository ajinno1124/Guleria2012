#!/usr/bin/env python
# -*- coding: utf-8 -*-

#check program in ipynb, then 実装
import numpy as np
import pandas as pd

def BindingEnergyL(NParamType,LParamType):
    AN=np.array([
		[6,5],
		[8,7],
		[14,13],
		[23,27],
		[39,49],
		[52,81],
		[82,125]
    ])

    f=open('BindingEnergyLY.csv','w',encoding="utf-8")
    f.write('This is a test\n')
    f.close()

    for i in range(len(AN)):
        NParamType="SK3"
        LParamType="LY1"
        df1=pd.read_csv(f"../data/Z{AN[i][0]}N{AN[i][1]}L0_{NParamType}NaN/Energy2.csv",comment="#")
        df2=pd.read_csv(f"../data/Z{AN[i][0]}N{AN[i][1]}L1_{NParamType}{LParamType}/Energy.csv",comment="#")
        
    
BindingEnergyL("SK3","LY1")