import numpy as np
import os
import EOMCCSD_Util as util
import EOMCCSD_AcesPy as ap

def get_Hbar(EnvVal,lprint):
    HbarType=EnvVal['HbarType']

    if HbarType=='ACES2':
       F,W,T=get_Hbar_ACES2(EnvVal,lprint)
    elif HbarType=='FILE':
       #F,W,T=get_Hbar_file(EnvVal)
       F,W,T=get_Hbar_file2(EnvVal)

    return F,W,T

def get_Hbar_file2(EnvVal):
    print('\n * Getting Hbar fron file')
    Nocc=EnvVal['Nocc']
    Nvrt=EnvVal['Nvrt']

    lread=True

    F={}
    W={}
    T={}

    T['2aa']=util.read_data('T2',lread)

    F['oo']=util.read_data('Foo',lread)
    F['vv']=util.read_data('Fvv',lread)
    F['ov']=util.read_data('Fov',lread)

    W['ovoo'] = util.read_data('Wovoo',lread) # [110] * ovoo 
    W['ooov'] = util.read_data('Wooov',lread) # [10] * ooov 
    W['vovv'] = util.read_data('Wvovv',lread) # [30] * vovv 
    W['oooo'] = util.read_data('Woooo',lread) # [53] * oooo 
    W['ovov'] = util.read_data('Wovov',lread) # [54] * ovov 
    W['ovvo'] = util.read_data('Wovvo',lread) # [56] * ovvo 
    W['vvvv'] = util.read_data('Wvvvv',lread) # [233] * vvvv
    W['oovv'] = util.read_data('Woovv',lread) # [?] * oovv
    W['vvvo'] = util.read_data('Wvvvo',lread) # [?] * vvvo

    return F,W,T

def get_Hbar_file(EnvVal):
    print('\n * Getting Hbar fron file')
    Nocc=EnvVal['Nocc']
    Nvrt=EnvVal['Nvrt']

    lread=True

    F={}
    W={}
    T={}

    F['oo_aa']=util.read_data('Foo_aa',lread)
    F['oo_bb']=util.read_data('Foo_bb',lread)
    F['vv_aa']=util.read_data('Fvv_aa',lread)
    F['vv_bb']=util.read_data('Fvv_bb',lread)
    F['vo_aa']=util.read_data('Fvo_aa',lread)
    F['vo_bb']=util.read_data('Fvo_bb',lread)

    # H(IJ,KA)
    W['5ooov_aaaa'] = util.read_data2('Wooov_aaaa',lread) #[7]   W(IJKA) 
    W['5ooov_bbbb'] = util.read_data2('Wooov_bbbb',lread) #[8]   W(ijka) 
    W['5oovo_abab'] = util.read_data2('Wooov_abba',lread) #[9]   W(IjAk) *
    W['5ooov_abab'] = util.read_data2('Wooov_abab',lread) #[10]  W(IjKa) 

    # H(IJ,AB)                                                 
    W['4vvoo_aaaa'] = util.read_data2('Woovv_aaaa',lread) #[14]  W(IJAB) << ABIJ
    W['4vvoo_bbbb'] = util.read_data2('Woovv_bbbb',lread) #[15]  W(ijab) << abij
    W['4vvoo_abab'] = util.read_data2('Woovv_abab',lread) #[16]  W(IjAb) << AbIj
                                                               
    # H(AI,BC) #FI X ME ordering H(CI,AB)?                     
    W['5vvvo_aaaa'] = util.read_data2('Wvovv_aaaa',lread) #[27]  W(ABCI)
    W['5vvvo_bbbb'] = util.read_data2('Wvovv_bbbb',lread) #[28]  W(abci)
    W['5vvov_abab'] = util.read_data2('Wvovv_baab',lread) #[29]  W(AbIc) *
    W['5vvvo_abab'] = util.read_data2('Wvovv_abab',lread) #[30]  W(AbCi)
                                                               
    # H(MN,IJ)                                                 
    W['4oooo_aaaa'] = util.read_data2('Woooo_aaaa',lread) #[51]  W(MNIJ) <- MNIJ
    W['4oooo_bbbb'] = util.read_data2('Woooo_bbbb',lread) #[52]  W(mnij) <- mnij
    W['4oooo_abab'] = util.read_data2('Woooo_abab',lread) #[53]  W(MnIj) <- MnIj 
                                                               
    # H(MB,JE)                                                 
    W['1vovo_aaaa'] = util.read_data2('Wovov_aaaa',lread) #[54]  W(MBJE) << EMBJ
    W['1vovo_bbbb'] = util.read_data2('Wovov_bbbb',lread) #[55]  W(mbje) << embj
    W['1vovo_aabb'] = util.read_data2('Wovov_abba',lread) #[56]  W(EMAI) << EMBJ
    W['1vovo_bbaa'] = util.read_data2('Wovov_baab',lread) #[57]  W(emAI) << emBJ [*]
    W['4vovo_abab'] = util.read_data2('Wovov_baba',lread) #[58]  W(mBjE) << EmBj
    W['4vovo_baba'] = util.read_data2('Wovov_abab',lread) #[59]  W(MbJe) << eMbJ

    # H(IA,JK) or H(KA,IJ)
    W['3ooov_aaaa'] = util.read_data2('Wovoo_aaaa',lread) #[107]  W(IJKA) 
    W['3ooov_bbbb'] = util.read_data2('Wovoo_bbbb',lread) #[108]  W(ijka) 
    W['3oovo_abab'] = util.read_data2('Wovoo_baab',lread) #[109]  W(IjAk) 
    W['3ooov_abab'] = util.read_data2('Wovoo_abab',lread) #[110]  W(IjKa) 
                                                                
    # H(AB,CI)                                                  
    W['3vvvo_aaaa'] = util.read_data2('Wvvvo_aaaa',lread) #[127]  W(ABCI)
    W['3vvvo_bbbb'] = util.read_data2('Wvvvo_bbbb',lread) #[128]  W(abci)
    W['3vvov_abab'] = util.read_data2('Wvvvo_abba',lread) #[129]  W(AbIc) *
    W['3vvvo_abab'] = util.read_data2('Wvvvo_abab',lread) #[130]  W(AbCi)
                                                                
    # H(AB,CD)                                                  
    W['4vvvv_aaaa'] = util.read_data2('Wvvvv_aaaa',lread) #[231]  W(ABCD) <- ABCD
    W['4vvvv_bbbb'] = util.read_data2('Wvvvv_bbbb',lread) #[232]  W(abcd) <- abcd
    W['4vvvv_abab'] = util.read_data2('Wvvvv_abab',lread) #[233]  W(AbCd) <- AbCd

    # T2
    T['vovo_aaaa'] = util.read_data2('Tvovo_aaaa',lread)
    T['vovo_bbbb'] = util.read_data2('Tvovo_bbbb',lread)
    T['vovo_abab'] = util.read_data2('Tvovo_abab',lread)
    T['ovov_bbaa'] = util.read_data('Tovov_bbaa',lread)


    AddDiag=True
    if (AddDiag):
       OrbEngA=util.read_data('SCFEVALA',True)
       OrbEngB=util.read_data('SCFEVALB',True)
       print('\n - Orbital energy')
       print(OrbEngA)
       print(OrbEngB)
       # add diagonal terms (orbital energies)
       Foo_aa=F['oo_aa']
       Foo_bb=F['oo_bb']
       for i in range(Nocc):
           Foo_aa[i][i]+=OrbEngA[i]
           Foo_bb[i][i]+=OrbEngB[i]

       Fvv_aa=F['vv_aa']
       Fvv_bb=F['vv_bb']
       for a in range(Nvrt):
           Fvv_aa[a][a]+=OrbEngA[Nocc+a]
           Fvv_bb[a][a]+=OrbEngB[Nocc+a]

    return F,W,T

def get_Hbar_ACES2(EnvVal,lprint):
    Nocc=EnvVal['Nocc']
    Nvrt=EnvVal['Nvrt']


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
    Wvvvo -= np.einsum('miab,me->abei',T2, Fov)      #t(miab)F(me)
    Wvvvo += np.einsum('if,abef->abei',T1, Wvvvv)    #t(if)W(abef)
    Wvvvo += np.einsum('mnab,mnei->abei',Tau, Goovo) #tau(mnab)<mn|ei>
     
    Wvvvo -= np.einsum('miaf') 
    Wvvvo 

def make_Hbar():
    # make Hbar 




