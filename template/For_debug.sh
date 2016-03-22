#!/bin/bash

this_dir=`pwd`
fc=gfortran
#... Installation of pyLikeDM ...
echo '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
echo 'Start installing pyLikeDM'
echo '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
cd ${this_dir}/src
rm -f likedm.pyf
f2py -m LikeDM -h likedm.pyf MathLib.f90 PYTHIA_PPPC4.f90 charge_data.f90 ReadTable.f90 gamma_dSphs.f90 charge_bkg.f90 charge_antip.f90 charge_lepton.f90 monitorLikeDM.f90
f2py -c --fcompiler=${fc} likedm.pyf MathLib.f90 PYTHIA_PPPC4.f90 charge_data.f90 ReadTable.f90 gamma_dSphs.f90 charge_bkg.f90 charge_antip.f90 charge_lepton.f90 monitorLikeDM.f90


echo 'End installing pyLikeDM'
echo 'Enjoy use!'
echo '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
