import os
import numpy as np
import EOMCCSD_Util as util
import EOMCCSD_Davidson as dav
import einsum as es

#np.set_printoptions(precision=5)

def make_Sigma(EnvVal,F,W,T,R0):
    Nocc=EnvVal['Nocc']
    Nvrt=EnvVal['Nvrt']

    R=expand_Vec1D(R0,Nocc,Nvrt)
    Z1=HR_Z1(EnvVal,F,W,R)
    Z2=HR_Z2(EnvVal,F,W,T,R)
    Z=np.concatenate((Z1,Z2),axis=0)
    return Z

def HR_Z1(EnvVal,F,W,R):
    Nocc=EnvVal['Nocc']
    Nvrt=EnvVal['Nvrt']
    Nov=Nocc*Nvrt

    R1a=R['1a']
    R2aa=R['2aa']

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
    Z1a += es.c_einsum("ae,ie->ia",   Fvv,   R1a)          
    Z1a -= es.c_einsum("mi,ma->ia",   Foo,   R1a)          
    Z1a += es.c_einsum("maei,me->ia", Wovvo, R1a)*2.0  
    Z1a += es.c_einsum("maie,me->ia",-Wovov, R1a)      

    ## Hsd 
    ## (DT2INT1A/DT2INT1B in ACES2)
    Z1a += es.c_einsum("me,miea->ia",   Fov,   R2aa)*2.0  #[93]
    Z1a += es.c_einsum("me,imea->ia",  -Fov,   R2aa)
    Z1a += es.c_einsum("amef,imef->ia", Wvovv, R2aa)*2.0  #[27] 
    Z1a += es.c_einsum("amfe,imef->ia",-Wvovv, R2aa)      #[30] 
    Z1a -= es.c_einsum("mnie,mnae->ia", Wooov, R2aa)*2.0  #[7]  
    Z1a -= es.c_einsum("nmie,mnae->ia",-Wooov, R2aa)      #[10] 

#   util.check_sum('Z1 - final',Z1a,2)
    Z1=Z1a.flatten()
    return Z1

def HR_Z2(EnvVal,F,W,T,R):
    Nocc=EnvVal['Nocc']
    Nvrt=EnvVal['Nvrt']
    Nov=Nocc*Nvrt

    R1a=R['1a']
    R2aa=R['2aa']
    T2aa=T['2aa']
    Foo=F['oo']  #9a
    Fvv=F['vv']  #92

    Wovoo=W['ovoo'] 
    Wooov=W['ooov'] 
    Wvovv=W['vovv'] 
    Woooo=W['oooo'] 
    Wovov=W['ovov'] 
    Wovvo=W['ovvo'] 
    Wvvvv=W['vvvv'] 
    Woovv=W['oovv'] 
    Wvvvo=W['vvvo'] 

    ## Hds (wo/three-body terms)
    ## (DT1INT2A/DT1INT2B in ACES2)
    Z2aa  = 0.0
    Z2aa += es.c_einsum("abej,ie->ijab", Wvvvo, R1a) 
    Z2aa -= es.c_einsum("mbij,ma->ijab", Wovoo, R1a) 

    ## Hdd (wo/three-body terms)
    ## (DT2INT2 in ACES2)
    Z2aa += es.c_einsum("be,ijae->ijab", Fvv, R2aa)
    Z2aa -= es.c_einsum("mj,imab->ijab", Foo, R2aa) 
    Z2aa += es.c_einsum("abef,ijef->ijab", Wvvvv, R2aa)*0.5  #[231]
    Z2aa += es.c_einsum("mnij,mnab->ijab", Woooo, R2aa)*0.5  #[51]
    Z2aa += es.c_einsum("mbej,imae->ijab", Wovvo, R2aa)*2.0  #[14] 
    Z2aa += es.c_einsum("mbej,miae->ijab",-Wovvo, R2aa)      #[14]  
    Z2aa -= es.c_einsum("mbje,miae->ijab", Wovov, R2aa)      #[16]    
    Z2aa -= es.c_einsum("maje,imeb->ijab", Wovov, R2aa)      #[16] 

    ## Three body terms in Hds and Hdd
    ## (FORMQ1/GFORMG2 in ACES2)
    Yvv   = 0.0
    Yvv  += es.c_einsum("bmfe,me->bf",   Wvovv, R1a)*2.0 
    Yvv  += es.c_einsum("bmef,me->bf",  -Wvovv, R1a)     
    Yvv  -= es.c_einsum("nmfe,nmbe->bf", Woovv, R2aa)    
    Z2aa += es.c_einsum("bf,ijaf->ijab", Yvv,   T2aa)    
    Yoo   = 0.0
    Yoo  -= es.c_einsum("nmje,me->nj",   Wooov, R1a)*2.0
    Yoo  -= es.c_einsum("mnje,me->nj",  -Wooov, R1a)    
    Yoo  -= es.c_einsum("nmef,jmef->nj", Woovv, R2aa)   
    Z2aa += es.c_einsum("nj,inab->ijab", Yoo,   T2aa)   

#   util.check_sum('Z2 - final',Z2aa,4)
    ## Some permutations are already applied in Hbar
    Z2aa += Z2aa.transpose([1,0,3,2])

    Z2=Z2aa.flatten()
    return Z2

def compress_Vec1D(V):
    V1a=V['1a'].flatten()
    V2aa=V['2aa'].flatten()

    V0=np.concatenate((V1a,V2aa),axis=0)
    return V0

def expand_Vec1D(V0,Nocc,Nvrt):
    Nov=Nocc*Nvrt
    Dim1=Nov
    Dim2=Nov*Nov

    V={}
    V['1a']=V0[:Dim1].reshape(Nocc,Nvrt)
    V['2aa']=V0[Dim1:Dim1+Dim2].reshape(Nocc,Nocc,Nvrt,Nvrt)
    return V

def print_vec(string,Z):
    print('\n - '+string)
    print(Z)

