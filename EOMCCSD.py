#import aces2py as a2
import numpy as np
#import scipy.linalg as sp
import os
import sys
import time
import math
import Base_Input as inp
import Base_Util as util
import Base_AcesPy as ap
import EOMCCSD_Hbar as hbar
import EOMCCSD_Guess as guess
import EOMCCSD_Davidson as dav



def main(argv):
    laces=False
    util.make_header('EOM-CCSD program')

    # setting variables up 
    EnvVal=inp.driver(argv)
    if (laces):
       EnvVal=ap.read_mol_info(EnvVal)

    #Smat,Hcore,Jmat,Kmat=ap.read_ints(Nbas)  # get integrals
    F,W,T = hbar.driver(EnvVal)

    R=guess.driver(EnvVal,F,W)

    dav.Diag_Davidson(EnvVal,F,W,T,R) 

    print('\n * End of the program')
    return


if __name__ == '__main__' :
    argv = sys.argv
    main(argv)


