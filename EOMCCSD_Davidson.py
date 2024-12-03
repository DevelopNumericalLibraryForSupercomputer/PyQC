import os
import numpy as np
import Base_Util as util
import Base_Param as param
import EOMCCSD_Sigma as sig
import EOMCCSD_Guess as guess
#np.set_printoptions(precision=5,suppress=True)
np.set_printoptions(precision=5)
au2eV=param.unit_conversion('energy','au','ev')

def Diag_Davidson(EnvVal,F,W,T,R):
    Nocc=EnvVal['NOCC']
    Nvrt=EnvVal['NVRT']
    Nroot=EnvVal['NROOT']  
    Nov=Nocc*Nvrt
    Rrow=Nov+Nov*Nov
    Rcol=1
    R0=np.zeros((Rrow,Rcol))
    for iroot in range(Nroot):
        print('\n * Solving root = '+str(iroot+1))
        R0[:,0]=R[:,iroot]
        Davidson_Iteration(EnvVal,F,W,T,R0)

def Davidson_Iteration(EnvVal,F,W,T,R):
    Nocc=EnvVal['NOCC']
    Nvrt=EnvVal['NVRT']
    RefOrb=EnvVal['REF_ORB']
    MaxIter=EnvVal['NMAX_ITER']
    EngTol=EnvVal['VTOL_ENG']
    NDimSubSp=EnvVal['NDIM_SUBSP'] #Maximum subspace dimension
    Nov=Nocc*Nvrt
    DimR=Nov+Nov*Nov
    NGuessSp=1
    NrootSub=1

    Hdiag=guess.make_Hdiag(EnvVal,F)

    Eng=[0.0]*NrootSub*NGuessSp

    print('\n * EOM-CCSD iteration, Max='+str(MaxIter))
    for Iter in range(MaxIter): 

        # make each R vector orthogonal with QR decomp. 
        Qval,Rval=np.linalg.qr(R)
        R=Qval

        # save energies
        EngOld=Eng[:NrootSub] 

        # make the subspace and diagonalize
        DimG=R.shape[1]
        G=np.zeros((DimG,DimG))
        Z=np.zeros((DimR,DimG))
        for i in range(DimG):
            Z[:,i]=sig.make_Sigma(EnvVal,F,W,T,R[:,i])
            for j in range(DimG):
                G[j,i]=np.dot(R[:,j],Z[:,i])
        #print('Subspace matrix (G)')
        #print(G)
        Eng,Evec=np.linalg.eig(G)

        # select root (by energy) 
        idx=Eng.argsort()[:NrootSub]
        Eng=np.real(Eng[idx])
        Evec=np.real(Evec[:,idx])

        # get new vector
        EngNorm=0.0
        for i in range(NrootSub):
            ResVec = np.dot(Z-Eng[i]*R, Evec[:,i])
            CorrVec=ResVec/(Eng[i]-Hdiag)
            dE=abs(Eng[i]-EngOld[i])
            print('   Iter.%3d :  E[%d] = %.10f    dE = %.10f ' % (Iter,i+1,Eng[i],dE))
            EngNorm += dE*dE
        EngNorm=np.sqrt(EngNorm)

        # check conv. / if not, update R
        if (EngNorm < EngTol) and (Iter>1): 
           print('\n * EOM-CCSD Converged')
           Eng2eV=Eng*au2eV
           for i in range(NrootSub):
               print(' - Energy for the state %d = %.10f  (%.5f eV)' % (i+1,Eng[i],Eng2eV[i]))
           return
        else:   
           if DimG >= NDimSubSp:
              R=np.dot(R,Evec)    # start the new subspace with the last estimate
              Eng=EngOld
           else: 
              R=np.c_[R,CorrVec]  # update R 

    return

