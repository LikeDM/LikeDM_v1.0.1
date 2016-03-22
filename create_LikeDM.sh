#!/bin/bash

this_dir=`pwd`
fco=gfortran
#... Installation of pyminuit ...
echo '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
echo 'Start installing pyminuit'
echo '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
echo 'Two ways to install pyminuit: enter use_pip or local'
echo 'Other keys for doing nothing.'
read set_pyminuit
case ${set_pyminuit} in
      use_pip)
      # If you have "pip" in your computer:
              sudo pip install iminuit 
              ;;
      local)
      # if you don't have pip: 
              tar -zxf CPPMinuit.tar.gz
              cd ${this_dir}/CPPMinuit
              ./configure
              make

              if [ -e "${this_dir}/CPPMinuit/src/.libs/liblcg_Minuit.so.0" ]; then
                 sudo cp ${this_dir}/CPPMinuit/src/.libs/liblcg_Minuit.so.0 /usr/lib
                 cd ${this_dir}/pyminuit
                 python setup.py clean
                 python setup.py install --home=${this_dir}/pyminuit --with-minuit=${this_dir}/CPPMinuit
                 echo "CPPMinuit is installed."
              else
                echo "CPPMinuit is failed to be installed."
                exit 1;
              fi


              ;;

	*)
              echo "User does not want to install pyminuit now."
esac
echo 'End installing pyminuit'
echo '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

#... Installation of pyLikeDM ...
echo '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
echo 'Start installing pyLikeDM'
echo '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

cd ${this_dir}/src
rm -f likedm.pyf
f2py -m LikeDM -h likedm.pyf MathLib.f90 PYTHIA_PPPC4.f90 charge_data.f90 ReadTable.f90 gamma_dSphs.f90 charge_bkg.f90 charge_antip.f90 charge_lepton.f90 monitorLikeDM.f90
f2py -c --fcompiler=${fco} likedm.pyf MathLib.f90 PYTHIA_PPPC4.f90 charge_data.f90 ReadTable.f90 gamma_dSphs.f90 charge_bkg.f90 charge_antip.f90 charge_lepton.f90 monitorLikeDM.f90


if [ -e "${this_dir}/src/LikeDM.so" ]; then
  echo "pyLikeDM is installed."
else
  echo "pyLikeDM is NOT installed."
  exit 1;
fi


echo 'End installing pyLikeDM'
echo 'Enjoy use!'
echo '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
