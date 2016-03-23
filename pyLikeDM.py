#!/usr/bin/python
import sys,os
from numpy import *
sys.path.append('./src/')
import LikeDM

try:
   import iminuit as minuit
except:
   import minuit
#   sys.path.append('./pyminuit/build/lib.linux-x86_64-2.7/') # perhaps depend on your computer



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
# Date: 2016-03-22           
###############################################################

'''
if len(sys.argv) < 1:  
    print "Need, at least, a name of ini-file!!!"
    print 'For example, "./pyLikeDM.py likedm.ini dnde.txt"'
    exit()

ini=read_setting(sys.argv[1])



TS_summary=[0.0,0.0,0.0,0.0,0.0,0.0,0.0]

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


if ini.use_pppc4:
  LikeDM.pythia_pppc4.readpppc4(ini.BR)
  if ini.decayDM:
    LikeDM.pythia_pppc4.dnde_from_pppc4(mx,mx)
  else:
    LikeDM.pythia_pppc4.dnde_from_pppc4(mx,2.0*mx)

  LikeDM.readtable.load_dnde('copy from PPPC4')
else:
  LikeDM.readtable.load_dnde(sys.argv[2]) 



#------------------------------------------------------------------#
LikeDM.readtable.game_initializing(dsphs_map,foutput)

if ini.use_dSphs:
  LikeDM.readtable.load_dsphs_likemap()
  #--- Gamma ray likelihood for dSphs ---
  dSphs=LikeDM.dsphs_gamma.dsphs_calclike(mx,sv)
  TS_summary[0]=dSphs



if ini.use_ep or ini.use_ap:

   LikeDM.readtable.loadgreen(WhatHalo,WhatGALPROP,ini.decayDM)

if ini.use_ep:
   LikeDM.charge_lepton.gen_epem_full_spectrum(mx,sv)

   def ep(a1,a2,a3,b1,b2,b3):
      chisq=0.0
      alpha=[a1,a2,a3]
      beta=[b1,b2,b3] 
      chisq+=LikeDM.charge_lepton.epemchisq('AMS02efr',alpha,beta)
      chisq+=LikeDM.charge_lepton.epemchisq('AMS02e+',alpha,beta)
      chisq+=LikeDM.charge_lepton.epemchisq('AMS02e-',alpha,beta)
      chisq+=LikeDM.charge_lepton.epemchisq('AMS02e+e-',alpha,beta)
      return chisq

   m = minuit.Minuit(ep, a1=1.0,limit_a1=(0.1,10.0), 
                     a2=1.0,limit_a2=(0.1,10.0),
                     a3=1.0,limit_a3=(0.1,10.0),
                     b1=0.0,limit_b1=(-0.5,0.5),
                     b2=0.0,limit_b2=(-0.5,0.5),
                     b3=0.0,limit_b3=(-0.5,0.5))
   m.migrad()  
   TSep=m.fval

   aep=[m.values["a1"], m.values["a2"],m.values["a3"]]
   bep=[m.values["b1"], m.values["b2"],m.values["b3"]]
   TS_summary[2]=LikeDM.charge_lepton.epemchisq('AMS02efr',aep,bep)
   TS_summary[3]=LikeDM.charge_lepton.epemchisq('AMS02e+',aep,bep)
   TS_summary[4]=LikeDM.charge_lepton.epemchisq('AMS02e-',aep,bep)
   TS_summary[5]=LikeDM.charge_lepton.epemchisq('AMS02e+e-',aep,bep)



if ini.use_ap:
   LikeDM.charge_antip.gen_ap_full_spectrum(mx,sv)

   def ap(alpha,beta):
      chisq=0.0
      chisq+=LikeDM.charge_antip.apchisq(alpha,beta)
      return chisq


   m = minuit.Minuit(ap, alpha=1.0,limit_alpha=(0.1,10.0), 
                     beta=0.0,limit_beta=(-0.5,0.5))
   m.migrad()   
   TSap=m.fval
   aap=m.values["alpha"]
   bap=m.values["beta"]
   TS_summary[6]=LikeDM.charge_antip.apchisq(aap,bap)


LikeDM.print_debug_info(mx,sv,aep,bep,aap,bap,TS_summary)

