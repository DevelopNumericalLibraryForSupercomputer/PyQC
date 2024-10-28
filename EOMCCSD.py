#import aces2py as a2
import numpy as np
#import scipy.linalg as sp
import math
import time
import sys
import os
import Base_Input as inp
import Base_Util as util
import EOMCCSD_Guess as guess
import EOMCCSD_Davidson as dav
import EOMCCSD_Hbar as hbar
import EOMCCSD_AcesPy as ap



def main(argv):
    laces=False

    # -----------------
    # acespy initialize 
    # -----------------
    if (laces):
       #if aces2py not in sys.modules:
       import aces2py as a2
       f=a2.init()
       #a2.buildhbar()
    
    util.make_header('EOM-CCSD program')

    # setting variables up 
    EnvVal=inp.driver(argv)
    if (laces):
       EnvVal=ap.read_mol_info(EnvVal)

    #Smat,Hcore,Jmat,Kmat=ap.read_ints(Nbas)  # get integrals
    lprint=False
    F,W,T,L = hbar.get_Hbar(EnvVal,lprint)

    R=guess.driver(EnvVal,F,W)

    dav.Diag_Davidson(EnvVal,F,W,T,L,R) 

    print('\n * End of the program')
    return


if __name__ == '__main__' :
    argv = sys.argv
    main(argv)


