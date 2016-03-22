#!/usr/bin/python
import sys,os
from numpy import *
sys.path.append('./pyCode/')
import LikeDM

try:
   import minuit
except:
   import iminuit as minuit




from Read_parameters import *

'''
###############################################################
#
# LikeDM python interface version 
# Please check the function by using 
#        "print LikeDM.__doc__"
#                                                          
# Author: Xiaoyuan Huang, Yue-Lin Sming Tsai, Qiang Yuan
# Email: huangxiaoyuan@gmail.com, smingtsai@gmail.com, yuanq@pmo.ac.cn
# Date: 2013-10-21    
# Date: 2015-08-30           
###############################################################

'''
if len(sys.argv) < 1:  
    print "Need, at least, a name of ini-file!!!"
    print 'For example, "./pyLikeDM.py likedm.ini dnde.txt"'
    exit()

ini=read_setting(sys.argv[1])




LikeDM.readtable.sets__dmdecay=ini.decayDM

mx=ini.mx
sv=ini.sv

LikeDM.readtable.sets__seebug=int(ini.seebug)
#LikeDM.readtable.sets__seebug=3

LikeDM.readtable.sets__halo=ini.WhatHalo
WhatHalo=ini.WhatHalo
WhatGALPROP=ini.WhatGALPROP
LikeDM.readtable.sets__use_dsphs=ini.use_dSphs
LikeDM.readtable.sets__use_ep=ini.use_ep
LikeDM.readtable.sets__use_ap=ini.use_ap
LikeDM.readtable.smod__ep=ini.epmod
LikeDM.readtable.smod__ap=ini.apmod


foutput= ini.output_name
dsphs_map=ini.dsphs_map

#------------------------------------------------------------------#
LikeDM.readtable.game_initializing(dsphs_map,foutput)


xr=linspace(log10(10.0),3.0, num=50, endpoint=True)
mx_list=10.0**xr[:]

xr=linspace(-30.0,-23.0, num=50, endpoint=True)
sv_list=10.0**xr[:]

LikeDM.readtable.sets__seebug=0

for mx in mx_list:
   for sv in sv_list:

     if ini.use_pppc4:
       LikeDM.pythia_pppc4.readpppc4(ini.BR)
       if ini.decayDM:
         LikeDM.pythia_pppc4.dnde_from_pppc4(mx,mx)
       else:
         LikeDM.pythia_pppc4.dnde_from_pppc4(mx,2.0*mx)

       LikeDM.readtable.load_dnde('copy from PPPC4')
     else:
       LikeDM.readtable.load_dnde(sys.argv[2]) 


     LikeDM.readtable.load_dsphs_likemap()
     #--- Gamma ray likelihood for dSphs ---
     dSphs=LikeDM.dsphs_gamma.dsphs_calclike(mx,sv)
     print mx,sv,dSphs




