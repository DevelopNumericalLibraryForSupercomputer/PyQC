#import aces2py as a2
import numpy as np
#import scipy.linalg as sp
import math
import time
import sys
import os
import EOMCCSD_Guess as guess
import EOMCCSD_Davidson as dav
import EOMCCSD_Hbar as hbar
import EOMCCSD_Util as util
import EOMCCSD_AcesPy as ap


def read_inp():
    # this will be moved to a separte file
    EnvVal={} 
    EnvVal['Irrep']=1
    EnvVal['Nbas']=11
    EnvVal['Nocc']=2
    EnvVal['Nvrt']=9
    EnvVal['RefOrb']='RHF'     
    EnvVal['MaxIter']=20
    EnvVal['Nroot']=1
    EnvVal['EngTol']=1.0E-5    #Energy Tolerence
    EnvVal['DavSubSpDim']=5    #Maximum subspace dimension
#   EnvVal['GuessType']='CIS_FILE' 
    EnvVal['GuessType']='Hdiag' 
    EnvVal['HbarType']='FILE'  
    
    return EnvVal 


def main():
    laces=False

    # -----------------
    # acespy initialize 
    # -----------------
    if (laces):
#      if aces2py not in sys.modules:
       import aces2py as a2
       f=a2.init()
       #a2.buildhbar()
    
    util.make_header('EOM-CCSD program')

    # setting variables up 
    EnvVal=read_inp()
    if (laces):
       EnvVal=ap.read_mol_info(EnvVal)

    #Smat,Hcore,Jmat,Kmat=ap.read_ints(Nbas)  # get integrals
    lprint=False
    F,W,T = hbar.get_Hbar(EnvVal,lprint)

    R=guess.driver(EnvVal,F,W)

    dav.Diag_Davidson(EnvVal,F,W,T,R) 

    return

main()



