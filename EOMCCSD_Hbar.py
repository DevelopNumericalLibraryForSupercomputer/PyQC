import numpy as np
import os
import Base_Util as util
import Base_AcesPy as ap

def get_Hbar(EnvVal,lprint):
    HbarType=EnvVal['HBAR_TYPE']

    if HbarType=='ACES2':
       print('This needs to be fixed.')
       #F,W,T=read_from_ACES2(EnvVal,lprint)
       exit(1)
    elif HbarType=='FILE':
       F,W,T,L=get_from_file(EnvVal)

    return F,W,T,L


def get_from_file(EnvVal):
    print('\n * Getting Hbar from file')
    Nocc=EnvVal['NOCC']
    Nvrt=EnvVal['NVRT']
    
    F={}
    W={}
    T={}
    L={}

    T['2ab']=util.read_data('T2',EnvVal)

    F['oo']=util.read_data('Foo',EnvVal)
    F['vv']=util.read_data('Fvv',EnvVal)
    F['ov']=util.read_data('Fov',EnvVal)

    W['ovoo'] = util.read_data('Wovoo',EnvVal) # [110] * ovoo 
    W['ooov'] = util.read_data('Wooov',EnvVal) # [10] * ooov 
    W['vovv'] = util.read_data('Wvovv',EnvVal) # [30] * vovv 
    W['oooo'] = util.read_data('Woooo',EnvVal) # [53] * oooo 
    W['ovov'] = util.read_data('Wovov',EnvVal) # [54] * ovov 
    W['ovvo'] = util.read_data('Wovvo',EnvVal) # [56] * ovvo 
    W['vvvv'] = util.read_data('Wvvvv',EnvVal) # [233] * vvvv
    W['vvvo'] = util.read_data('Wvvvo',EnvVal) # [?] * vvvo
    L['oovv'] = util.read_data('Loovv',EnvVal) # [?] * oovv

    return F,W,T,L


def read_from_ACES2(EnvVal,lprint):
    Nocc=EnvVal['NOCC']
    Nvrt=EnvVal['NVRT']

    F={}
    print('here')
    F['oo']=ap.get_Vec2D(Nocc,Nvrt,'Foo',lprint)
    F['vv']=ap.get_Vec2D(Nocc,Nvrt,'Fvv',lprint)
    F['ov']=ap.get_Vec2D(Nocc,Nvrt,'Fov',lprint)

    # T2                     
    T={}
    T['ov']=ap.get_Vec2D(Nocc,Nvrt,'Tov',lprint)
    T['ovov']=ap.get_Vec4D(Nocc,Nvrt,'Tovov',lprint)
    

    AddDiag=True
    if (AddDiag):
       OrbEng=util.read_data('SCFEVALA',True)
       print('OrbEng')
       print(OrbEng)
       # add diagonal terms (orbital energies)
       Foo_aa=F['oo_aa']
       for i in range(Nocc):
           Foo_aa[i][i]+=OrbEng[i]
           
       Fvv_aa=F['vv_aa']
       for a in range(Nvrt):
           Fvv_aa[a][a]+=OrbEng[Nocc+a]

    return F,W,T

def get_tau(T2,T1):
    # Tau(RHF) = T2 + T1*T1
    # TildeTau(RHF)+0.5*T1*T1 = (T2 + 0.5*T1*T1) + (0.5*T1*T1) = Tau(RHF)
    tau  = T2 
    tau += np.einsum('ia,jb->ijab',T1,T1)
    return tau

def get_Loovv(Goovv):
    Loovv  = 2.0*Goovv
    Loovv -= np.einsum("pqrs->pqsr",Goovv)  #swap 3-4
    return Loovv

def get_Lvovv(Gvovv):
    Lvovv  = 2.0*Gvovv
    Lvovv -= np.einsum("pqrs->pqsr",Gvovv)  #swap 3-4
    return Lvovv

def get_Looov(Gooov):
    Looov  = 2.0*Gooov
    Looov -= np.einsum("pqrs->qprs",Gooov)  #swap 1-2
    return Looov


def get_Fvv(fvv,fov,Lvovv,Loovv,T1,Tau):
    # F(ae) = f(ae) - t(ma)f(me) + t(mf)<ma||fe> - tau(mnaf)<mn|ef> 
    Fvv  = fvv 
    Fvv -= np.einsum("ma,me->ae",     T1, fov)
    Fvv += np.einsum("mf,amef->ae",   T1, Lvovv)
    Fvv -= np.einsum("mnaf,mnef->ae", Tau, Loovv) 
    return Fvv

def get_Foo(foo,fov,Looov,Loovv,T1):
    # F(mi) = f(mi) + t(ei)f(me) + t(en)<mn||ie> + 1/2 tau(inef)<mn|ef> 
    Foo  = foo
    Foo += np.einsum('ie,me->mi',     T1, fov)
    Foo += np.einsum('ne,mnie->mi',   T1, Looov)
    Foo += np.einsum('inef,mnef->mi', Tau, Loovv)
    return Foo

def get_Fov(fov,Loovv,T1):
    # F(me) = f(me) + t(fn)<mn||ef>
    Fov = fov
    Fov += np.einsum('nf,mnef->me',   T1, Loovv)
    return Fov

def get_Woooo(Goooo,Gooov,Goovv,T1,Tau):
    Woooo  = Goooo  #G(mnij)
    Woooo += np.einsum('je.mnie->mnij', T1, Gooov)
    Woooo += np.einsum('ie.mnej->mnij', T1, Gooov) 
    Woooo += np.einsum('ijef,mnef->mnij', Tau, Goovv)
    return Woooo

def get_Wvvvv(Gvovv,Goovv,T1,Tau):
    Wvvvv  = Gvvvv  #G(abef)
    Wvvvv -= np.einsum('mb.amef->abef', T1, Gvovv)
    Wvvvv -= np.einsum('ma.bmfe->abef', T1, Gvovv) 
    Wvvvv += np.einsum('mnab,mnef->abef', Tau, Goovv)
    return Wvvvv

def get_Wovvo(Govvo,Govvv,Goovo,Goovv):
    # W(mbej) = <mb||ej> + t(jf)<mb||ef> - t(nb)<mn||ej> - (t(jnfb)+t(jf)t(nb))<mn||ef>
    Wovvo  = Govvo  #G(mbej)
    Wovvo += np.einsum('jf,mbef->mbej', T1, Govvv)
    Wovvo -= np.einsum('nb,mnej->mbej', T1, Goovo)
    Wovvo -= np.einsum('jnfb,mnef->mbej', Tau, Goovv)
    Wovvo += np.einsum('njfb,mnef->mbej', T2, Goovv)
    return Wovvo

def get_Wovov(Wovvo):
    # W(mbje) = -W(mbej)
    return -Wovvo.swapaxes(2,3) 
 
def get_Wooov(Gooov,Goovv,T1):
    Wooov  = Gooov  #G(mnie)
    Wooov += np.einsum('if,mnfe->mnie', T1, Goovv)

def get_Wvovv():
    Wvovv  = Gvovv  #G(amef)
    Wvovv += np.einsum('na,nmef->amef', T1, Goovv)
    return Wvovv

def get_Wovoo(Govoo):
    Wovoo  = Govoo  #G(mbij)
    Wovoo += np.einsum('ijeb,me->mbij', T2, Fov)
    Wovoo -= np.einsum('nb,mnij->mbij', T2, Woooo)
    Wovoo += np.einsum('ijef,mbef->mbij', Tau, Govvv)

    Wovoo += np.einsum('jnbe,mnie->mbij', T2, Looov)  #P(ij)t(jnbe)<mn||ie>
    Wovoo -= np.einsum('njbe,mnie->mbij', T2, Gooov)  
    Wovoo -= np.einsum('ineb,mnej->mbij', T2, Gooov)

    Wovoo += np.einsum('ie,mbej->mbij',T1, Govov)     #P(ij)t(ie)<mb||ej>
    Wovoo += np.einsum('je,mbie->mbij',T1, Govov)

    tmp    = np.einsum('njbf,mnef->mbej',T2, Goovv)   #P(ij)t(ie){...}
    tmp   -= np.einsum('jnbf,mnef->mbej',T2, Goovv)
    Wovoo += np.einsum('ie,mbej->mbij',T1, tmp)
    tmp    = np.einsum('infb,mnfe->mbie',T2, Goovv)
    Wovoo += np.einsum('je,mbie->mbij',T1, tmp)
    return Wovoo

def get_Wvvvo():
    Wvvvo  = Gvvvo  #G(abei)
    Wvvvo -= np.einsum('miab,me->abei',T2, Fov)       #t(miab)F(me)
    Wvvvo += np.einsum('if,abef->abei',T1, Wvvvv)     #t(if)W(abef)
    Wvvvo += np.einsum('mnab,mnei->abei',Tau, Goovo)  #tau(mnab)<mn|ei>
                                                      
    Wvvvo -= np.einsum('miaf,mbef->abei',T2, Wovvv)   #P(ab)t(miaf)<mb||ef>
    Wvvvo -= np.einsum('mibf,amef->abei',T2. Wvovv) 
    Wvvvo += np.einsum('imbf,amef->abei',T2. Wvovv) 

    Wvvvo -= np.einsum('ma,mbei->abei',T1, Wovvo)     #P(ab)t(ma)<mb||ei>
    Wvvvo -= np.einsum('mb,amei->abei',T1, Wvovo)

    tmp    = np.einsum('nibf,mnef->mbei',T2, Woovv)   #P(ab)t(ma)t(nibf)<mn||ef>
    tmp   -= np.einsum('inbf,mnef->mbei',T2, Woovv)
    Wvvvo += np.einsum('ma,mbei->abei',T1, tmp)
    tmp    = np.einsum('niaf,nmef->amei',T2, Woovv)
    Wvvvo  = np.einsum('mb,amei->abei',T1, tmp)
    return Wvvvo


