import os
import numpy as np
import Base_Util as util
import EOMCCSD_Sigma as sig
import EOMCCSD_Guess as guess
#np.set_printoptions(precision=5,suppress=True)
np.set_printoptions(precision=5)

def Diag_Davidson(EnvVal,F,W,T,L,R):
    Nocc=EnvVal['NOCC']
    Nvrt=EnvVal['NVRT']
    RefOrb=EnvVal['REF_ORB']
    MaxIter=EnvVal['NMAXITER']
    Nroot=EnvVal['NROOT']  #for now, Nroot=1
    EngTol=EnvVal['VTOL_ENG']
    DavSubSpDim=EnvVal['NDIM_SUBSP'] #Maximum subspace dimension
    Nov=Nocc*Nvrt
    DimR=Nov+Nov*Nov
    NGuessSp=2

    Hdiag=guess.make_Hdiag(EnvVal,F)

    Eng=[0.0]*Nroot*NGuessSp

    print('\n * EOM-CCSD iteration, Max='+str(MaxIter))
    for Iter in range(MaxIter): 

        # make each R vector orthogonal with QR decomp. 
        Qval,Rval=np.linalg.qr(R)
        R=Qval

        # save energies
        Eng_old=Eng[:Nroot] 

        # make the subspace and diagonalize
        DimG=R.shape[1]
        G=np.zeros((DimG,DimG))
        Z=np.zeros((DimR,DimG))
        for i in range(DimG):
            Z[:,i]=sig.make_Sigma(EnvVal,F,W,T,L,R[:,i])
            for j in range(DimG):
                G[j,i]=np.dot(R[:,j],Z[:,i])
        #print('Subspace matrix (G)')
        #print(G)
        Eng,Evec=np.linalg.eig(G)

        # select root (by energy) 
        idx=Eng.argsort()[:Nroot]
        Eng=np.real(Eng[idx])
        Evec=np.real(Evec[:,idx])

        # get new vector
        Enorm=0.0
        for i in range(Nroot):
            ResVec = np.dot(Z-Eng[i]*R, Evec[:,i])
            CorrVec=ResVec/(Eng[i]-Hdiag)
            dE=abs(Eng[i]-Eng_old[i])
            print('   Iter.%3d :  E[%d] = %.10f    dE = %.10f ' % (Iter,i+1,Eng[i],dE))
            Enorm += dE*dE
        Enorm=np.sqrt(Enorm)

        # check conv. / if not, update R
        if (Enorm < EngTol) and (Iter>1): 
           print('\n * EOM-CCSD Converged')
           for i in range(Nroot):
               print(' - Energy for the state %d = %.10f ' % (i+1,Eng[i]))
           return
        else:   
           if DimG >= DavSubSpDim:
              R=np.dot(R,Evec)    # start the new subspace with the last estimate
              Eng=Eng_old
           else: 
              R=np.c_[R,CorrVec]  # update R 

    return

