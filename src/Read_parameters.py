#!/usr/bin/python

###############################################################
#
# LikeDM python interface to read parameters 
# Please check the function by using 
#        "print LikeDM.__doc__"
#                                                          
# Author: Yue-Lin Sming Tsai
# Email: smingtsai@gmail.com
# Date: 2013-10-21    
# Date: 2014-08-13           
###############################################################

class read_setting:
  def __init__(self,ini_file):

     lines = open(ini_file).readlines()
     self.BR = [0.0 for i in range(28)]

     for line_ in lines:
       line=line_.strip()
       if not line.startswith("#"):
           if line.startswith("decayDM"):
             thisline=line.split('=')[1].split('#')[0].strip() 
             if thisline in ('Y', 'T', 'yes', 'true', 'True'):
               self.decayDM=True
             else:
               self.decayDM=False

           if line.startswith("DMmass"):
             self.mx=float(line.split('=')[1].split('#')[0])

           if line.startswith("sigmav"):
             self.sv=float(line.split('=')[1].split('#')[0])

           if line.startswith("decay_time"):
             self.tau=float(line.split('=')[1].split('#')[0])

           if line.startswith("seebug"):
             self.seebug=line.split('=')[1].split('#')[0].strip()

           if line.startswith("WhatHalo"):
             self.WhatHalo=int(line.split('=')[1].split('#')[0])

           if line.startswith('WhatGALPROP'):
             self.WhatGALPROP=int(line.split('=')[1].split('#')[0])

           if line.startswith("use_dSphs"):
             thisline=line.split('=')[1].split('#')[0].strip() 
             if thisline in ('Y', 'T', 'yes', 'true', 'True'):
               self.use_dSphs=True
             else:
               self.use_dSphs=False

           if line.startswith("use_ep"):
             thisline=line.split('=')[1].split('#')[0].strip() 
             if thisline in ('Y', 'T', 'yes', 'true', 'True'):
               self.use_ap=True
             else:
               self.use_ap=False

           if line.startswith("use_ep"):
             thisline=line.split('=')[1].split('#')[0].strip() 
             if thisline in ('Y', 'T', 'yes', 'true', 'True'):
               self.use_ep=True
             else:
               self.use_ep=False

           if line.startswith("use_DD"):
             thisline=line.split('=')[1].split('#')[0].strip() 
             if thisline in ('Y', 'T', 'yes', 'true', 'True'):
               self.use_DD=True
             else:
               self.use_DD=False

           if line.startswith("use_pppc4"):
             thisline=line.split('=')[1].split('#')[0].strip()
             if thisline in ('Y', 'T', 'yes', 'true', 'True'):
               self.use_pppc4=True
             else:
               self.use_pppc4=False

           if line.startswith("epmod"):
             self.epmod=float(line.split('=')[1].split('#')[0])

           if line.startswith("apmod"):
             self.apmod=float(line.split('=')[1].split('#')[0])

           if line.startswith("dsphs_map"):
             self.dsphs_map=line.split('=')[1].split('#')[0].strip()

           if line.startswith("output_name"):
             self.output_name=line.split('=')[1].split('#')[0].strip()


           if line.startswith("BR_"):
             for i in xrange(28):
               if 'BR_'+str(i+1) == line.split('=')[0].strip(): 
                 self.BR[i]=float(line.split('=')[1].split('#')[0])

           if line.startswith("sigsip"):
             self.sigsip=float(line.split('=')[1].split('#')[0])

           if line.startswith("sigsin"):
             self.sigsin=float(line.split('=')[1].split('#')[0])

           if line.startswith("sigsdp"):
             self.sigsdp=float(line.split('=')[1].split('#')[0])

           if line.startswith("sigsdn"):
             self.sigsdn=float(line.split('=')[1].split('#')[0])

     if self.decayDM: self.sv=1.0/self.tau

