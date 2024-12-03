import os
import sys
import numpy as np
import Base_Util as util
import EOMCCSD_Davidson as dav
#sys.path.append("/mnt/c/Ubuntu/Workspace/Code/KISTI/Einsum")
#import einsum as es
#np.set_printoptions(precision=5)

# Equations follow the Chem. Phys. Lett. 248, 189 (1996) by Gwaltney et al.
# with the spin integration for RHF.


def make_Sigma(EnvVal,F,W,T,R0):
    Nocc=EnvVal['NOCC']
    Nvrt=EnvVal['NVRT']

    R=expand_Vec1D(R0,Nocc,Nvrt)
    Z1=HR_Z1(EnvVal,F,W,R)
    Z2=HR_Z2(EnvVal,F,W,T,R)
    Z=np.concatenate((Z1,Z2),axis=0)
    return Z

def HR_Z1(EnvVal,F,W,R):
    Nocc=EnvVal['NOCC']
    Nvrt=EnvVal['NVRT']
    Nov=Nocc*Nvrt

    R1a=R['1a']
    R2ab=R['2ab']

    Foo=F['oo'] #91 
    Fvv=F['vv'] #92 
    Fov=F['ov'] #93
    Wovoo=W['ovoo'] 
    Wooov=W['ooov'] 
    Wvovv=W['vovv'] 
    Wovov=W['ovov'] 
    Wovvo=W['ovvo'] 

    ## Hss 
    ## (DFT1INT1/DT1INT1 in ACES2)
    Z1a = 0.0
    Z1a += np.einsum("ae,ie->ia",   Fvv,   R1a)          
    Z1a -= np.einsum("mi,ma->ia",   Foo,   R1a)          
    Z1a += np.einsum("maei,me->ia", Wovvo, R1a)*2.0  
    Z1a -= np.einsum("maie,me->ia", Wovov, R1a)      

    ## Hsd 
    ## (DT2INT1A/DT2INT1B in ACES2)
    Z1a += np.einsum("me,imae->ia",   Fov,   R2ab)*2.0  #[93]
    Z1a -= np.einsum("me,miae->ia",   Fov,   R2ab)
    Z1a += np.einsum("amef,imef->ia", Wvovv, R2ab)*2.0  #[27] 
    Z1a -= np.einsum("amfe,imef->ia", Wvovv, R2ab)      #[30] 
    Z1a -= np.einsum("mnie,mnae->ia", Wooov, R2ab)*2.0  #[7]  
    Z1a += np.einsum("nmie,mnae->ia", Wooov, R2ab)      #[10] 

#   util.check_sum('Z1 - final',Z1a,2)
    Z1=Z1a.flatten()
    return Z1

def HR_Z2(EnvVal,F,W,T,R):
    Nocc=EnvVal['NOCC']
    Nvrt=EnvVal['NVRT']
    Nov=Nocc*Nvrt

    R1a=R['1a']
    R2ab=R['2ab']
    T2ab=T['2ab']
    Foo=F['oo']  #9a
    Fvv=F['vv']  #92

    Wovoo=W['ovoo'] 
    Wooov=W['ooov'] 
    Wvovv=W['vovv'] 
    Woooo=W['oooo'] 
    Wovov=W['ovov'] 
    Wovvo=W['ovvo'] 
    Wvvvv=W['vvvv'] 
    Wvvvo=W['vvvo'] 
    Woovv=W['oovv'] 

    ## Hds (wo/three-body terms)
    ## (DT1INT2A/DT1INT2B in ACES2)
    Z2ab  = 0.0
    Z2ab -= np.einsum("mbij,ma->ijab", Wovoo, R1a)         # P(ab)[W(maij)R(mb)]
    Z2ab += np.einsum("abej,ie->ijab", Wvvvo, R1a)         # P(ij)[W(abej)R(ie)]

    ## Hdd (wo/three-body terms)
    ## (DT2INT2 in ACES2)
    Z2ab += np.einsum("be,ijae->ijab", Fvv, R2ab)          # P(ab)[F(be)R(ijae)]
    Z2ab -= np.einsum("mj,imab->ijab", Foo, R2ab)          #-P(ij)[F(mj)R(imab)]
    Z2ab += np.einsum("abef,ijef->ijab", Wvvvv, R2ab)*0.5  # W(abef)R(ijef)*0.5
    Z2ab += np.einsum("mnij,mnab->ijab", Woooo, R2ab)*0.5  # W(mnij)R(mnab)*0.5
    # P(ab)P(ij)[W(mbej)R(imae)] 
    Z2ab += np.einsum("mbej,imae->ijab", Wovvo, R2ab)*2.0    
    Z2ab -= np.einsum("mbej,miae->ijab", Wovvo, R2ab)        
    Z2ab -= np.einsum("mbje,imae->ijab", Wovov, R2ab)          
    Z2ab -= np.einsum("maje,imeb->ijab", Wovov, R2ab)       

    ## Hds (Three body terms)
    ## (FORMQ1/GFORMG2 in ACES2)
    Yvv   = 0.0
    Yvv  += np.einsum("amfe,me->af",   Wvovv, R1a)*2.0 
    Yvv  -= np.einsum("amef,me->af",   Wvovv, R1a)     
    Z2ab += np.einsum("af,ijfb->ijab", Yvv,   T2ab)    

    Yoo   = 0.0
    Yoo  -= np.einsum("nmie,me->ni",   Wooov, R1a)*2.0
    Yoo  += np.einsum("mnie,me->ni",   Wooov, R1a)    
    Z2ab += np.einsum("ni,njab->ijab", Yoo,   T2ab)   

    ## Hdd (Three body terms) 
    ## (FORMQ1/GFORMG2 in ACES2)
    Yvv   = 0.0
    Yvv  -= np.einsum("nmfe,nmae->af", Woovv, R2ab)    
    Z2ab += np.einsum("af,ijfb->ijab", Yvv,   T2ab)    

    Yoo   = 0.0
    Yoo  -= np.einsum("nmfe,imfe->ni", Woovv, R2ab)
    Z2ab += np.einsum("ni,njab->ijab", Yoo,   T2ab)   

#   util.check_sum('Z2 - final',Z2ab,4)
    Z2ab += Z2ab.transpose([1,0,3,2])

    Z2=Z2ab.flatten()
    return Z2

def compress_Vec1D(V):
    V1a=V['1a'].flatten()
    V2ab=V['2ab'].flatten()
    V0=np.concatenate((V1a,V2ab),axis=0)
    return V0

def expand_Vec1D(V0,Nocc,Nvrt):
    Nov=Nocc*Nvrt
    Dim1=Nov
    Dim2=Nov*Nov

    V={}
    V['1a']=V0[:Dim1].reshape(Nocc,Nvrt)
    V['2ab']=V0[Dim1:Dim1+Dim2].reshape(Nocc,Nocc,Nvrt,Nvrt)
    return V

def print_vec(string,Z):
    print('\n - '+string)
    print(Z)

