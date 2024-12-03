import numpy as np
import os
import sys
import Base_Input as inp
import Base_Util as util
import Base_AcesPy as ap

def driver(EnvVal):
    HbarType=EnvVal['HBAR_TYPE']

    if HbarType=='ACES2':
       print('This routine needs to be fixed.')
       #F,W,T=read_from_ACES2(EnvVal)
       exit(1)
    elif HbarType=='FILE':
       F,W,T=get_from_file(EnvVal)
    elif HbarType=='CALC':
       calc_Hbar(EnvVal) 
       return
    elif HbarType=='INTS':
       get_Ints(EnvVal)
       return

    return F,W,T

def test_results(EnvVal):
    util.compare_two_vectors(EnvVal,'Data-T1.txt','ACES2-T1.txt')
    util.compare_two_vectors(EnvVal,'Data-T2.txt','ACES2-T2.txt')
    return

def get_from_file(EnvVal):
    print('\n * Getting Hbar from file')
    Nocc=EnvVal['NOCC']
    Nvrt=EnvVal['NVRT']
    
    F={}
    W={}
    T={}

    T['2ab']=util.read_data('Data-T2.txt',EnvVal)

    F['oo']=util.read_data('Data-Foo.txt',EnvVal)
    F['vv']=util.read_data('Data-Fvv.txt',EnvVal)
    F['ov']=util.read_data('Data-Fov.txt',EnvVal)

    W['ovoo'] = util.read_data('Data-Wovoo.txt',EnvVal) # [110] * ovoo 
    W['ooov'] = util.read_data('Data-Wooov.txt',EnvVal) # [10] * ooov 
    W['vovv'] = util.read_data('Data-Wvovv.txt',EnvVal) # [30] * vovv 
    W['oooo'] = util.read_data('Data-Woooo.txt',EnvVal) # [53] * oooo 
    W['ovov'] = util.read_data('Data-Wovov.txt',EnvVal) # [54] * ovov 
    W['ovvo'] = util.read_data('Data-Wovvo.txt',EnvVal) # [56] * ovvo 
    W['oovv'] = util.read_data('Data-Woovv.txt',EnvVal) # [?] * oovv
    W['vvvv'] = util.read_data('Data-Wvvvv.txt',EnvVal) # [233] * vvvv
    W['vvvo'] = util.read_data('Data-Wvvvo.txt',EnvVal) # [?] * vvvo

    return F,W,T


def read_from_ACES2(EnvVal):
    Nocc=EnvVal['NOCC']
    Nvrt=EnvVal['NVRT']
    lprint=False

    F={}
    print('here')
    F['oo']=ap.get_Vec2D(Nocc,Nvrt,'Foo',EnvVal,lprint)
    F['vv']=ap.get_Vec2D(Nocc,Nvrt,'Fvv',EnvVal,lprint)
    F['ov']=ap.get_Vec2D(Nocc,Nvrt,'Fov',EnvVal,lprint)

    # T2                     
    T={}
    T['ov']=ap.get_Vec2D(Nocc,Nvrt,'Tov',EnvVal,lprint)
    T['ovov']=ap.get_Vec4D(Nocc,Nvrt,'Tovov',EnvVal,lprint)
    

    AddDiag=True
    if (AddDiag):
       OrbEng=util.read_data('Data-SCFEVALA.txt',True)
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


def get_Ints(EnvVal):
    EnvVal['Hbar_SubMO']='FILE'
    util.get_SubMO('oooo',EnvVal)
    util.compare_two_vectors(EnvVal,'Int-Goooo.txt','Int-Goooo.txt')
    util.get_SubMO('ooov',EnvVal)
    util.compare_two_vectors(EnvVal,'Int-Gooov.txt','Int-Gooov.txt')
    util.get_SubMO('ovov',EnvVal)
    util.compare_two_vectors(EnvVal,'Int-Govov.txt','Int-Govov.txt')


def calc_Hbar(EnvVal):
    print('\n * Calculate Hbar elements')
    Nocc=EnvVal['NOCC']
    Nvrt=EnvVal['NVRT']

    get_tau(EnvVal)

    lget_Lpqrs=True
    if (lget_Lpqrs):
       Gooov=util.get_SubMO('ooov',EnvVal)
       get_antisymmetrization(Gooov,'pq','Int-Looov.txt')
       Goovv=util.get_SubMO('oovv',EnvVal)
       get_antisymmetrization(Goovv,'rs','Int-Loovv.txt')
       Gvovv=util.get_SubMO('vovv',EnvVal)
       get_antisymmetrization(Gvovv,'rs','Int-Lvovv.txt')

    lget_Fpq=True
    if (lget_Fpq):
       get_Fvv(EnvVal)
       get_Foo(EnvVal)
       get_Fov(EnvVal)

    lget_Wpqrs=True
    if (lget_Wpqrs):
       get_Woooo(EnvVal)
       get_Wvvvv(EnvVal)
       get_Wovvo(EnvVal)
       get_Wovov(EnvVal) 
       get_Woovv(EnvVal)  
       get_Wooov(EnvVal)
       get_Wvovv(EnvVal)
       get_Wovoo(EnvVal)
       get_Wvvvo(EnvVal)

    return

def get_tau(EnvVal):
    # Tau(RHF) = T2 + T1*T1
    # TildeTau(RHF)+0.5*T1*T1 = (T2 + 0.5*T1*T1) + (0.5*T1*T1) = Tau(RHF)
    T1=util.read_data('Data-T1.txt',EnvVal)
    T2=util.read_data('Data-T2.txt',EnvVal)

    tau  = T2 
    tau += np.einsum('ia,jb->ijab',T1,T1)

    if (EnvVal['HBAR_OUT']=='FILE'):
       util.write_data('Data-tau.txt',tau,EnvVal)
       return
    else:
       return tau

def get_antisymmetrization(Gpqrs,*string):
    Lpqrs  = 2.0*Gpqrs
    if string[0]=='pq':
       Lpqrs -= np.einsum("pqrs->qprs",Gpqrs)  #swap p,q
    elif string[0]=='rs':
       Lpqrs -= np.einsum("pqrs->pqsr",Gpqrs)  #swap r,s

    if (EnvVal['HBAR_OUT']=='FILE'):
       util.write_data(string[1],Lpqrs,EnvVal)
       return
    else:
       return Lpqrs


def get_Fvv(EnvVal):
    # F(ae) = f(ae) - t(ma)f(me) + t(mf)<ma||fe> - tau(mnaf)<mn|ef> 
    T1=util.read_data('Data-T1.txt',EnvVal)
    tau=util.read_data('Data-tau.txt',EnvVal)
    fvv=util.read_data('Int-fvv.txt',EnvVal)
    fov=util.read_data('Int-fov.txt',EnvVal)
    Loovv=util.read_data('Int-Loovv.txt',EnvVal)
    Lvovv=util.read_data('Int-Lvovv.txt',EnvVal)

    Fvv  = fvv 
    Fvv -= np.einsum("ma,me->ae",     T1, fov)
    Fvv += np.einsum("mf,amef->ae",   T1, Lvovv)
    Fvv -= np.einsum("mnaf,mnef->ae", tau, Loovv) 

    if (EnvVal['HBAR_OUT']=='FILE'):
       util.write_data('Data-Fvv.txt',Fvv,EnvVal)
       return
    else:
       return Fvv 


def get_Foo(EnvVal):
    # F(mi) = f(mi) + t(ei)f(me) + t(en)<mn||ie> + 1/2 tau(inef)<mn|ef> 
    T1=util.read_data('Data-T1.txt',EnvVal)
    tau=util.read_data('Data-tau.txt',EnvVal)
    foo=util.read_data('Int-foo.txt',EnvVal)
    fov=util.read_data('Int-fov.txt',EnvVal)
    Looov=util.read_data('Int-Looov.txt',EnvVal)
    Loovv=util.read_data('Int-Loovv.txt',EnvVal)

    Foo  = foo
    Foo += np.einsum('ie,me->mi',     T1, fov)
    Foo += np.einsum('ne,mnie->mi',   T1, Looov)
    Foo += np.einsum('inef,mnef->mi', tau, Loovv)

    if (EnvVal['HBAR_OUT']=='FILE'):
       util.write_data('Data-Foo.txt',Foo,EnvVal)
       return
    else:
       return Foo


def get_Fov(EnvVal):
    # F(me) = f(me) + t(fn)<mn||ef>
    T1=util.read_data('Data-T1.txt',EnvVal)
    fov=util.read_data('Int-fov.txt',EnvVal)
    Loovv=util.read_data('Int-Loovv.txt',EnvVal)

    Fov = fov
    Fov += np.einsum('nf,mnef->me',   T1, Loovv)

    if (EnvVal['HBAR_OUT']=='FILE'):
       util.write_data('Data-Fov.txt',Fov,EnvVal)
       return
    else:
       return Fov


def get_Woooo(EnvVal):
    T1=util.read_data('Data-T1.txt',EnvVal)
    tau=util.read_data('Data-tau.txt',EnvVal)
    Goooo=util.get_SubMO('oooo',EnvVal)
    Gooov=util.get_SubMO('ooov',EnvVal)
    Goovo=util.get_SubMO('oovo',EnvVal)
    Goovv=util.get_SubMO('oovv',EnvVal)

    Woooo  = Goooo  #G(mnij)
    Woooo += np.einsum('je,mnie->mnij', T1, Gooov)
    Woooo += np.einsum('ie,mnej->mnij', T1, Goovo) #check 
    Woooo += np.einsum('ijef,mnef->mnij', tau, Goovv)

    if (EnvVal['HBAR_OUT']=='FILE'):
       util.write_data('Data-Woooo.txt',Woooo,EnvVal)
       return
    else:
       return Woooo


def get_Wvvvv(EnvVal):
    T1=util.read_data('Data-T1.txt',EnvVal)
    tau=util.read_data('Data-tau.txt',EnvVal)
    Goovv=util.get_SubMO('oovv',EnvVal)
    Gvovv=util.get_SubMO('vovv',EnvVal)
    Gvvvv=util.get_SubMO('vvvv',EnvVal)

    Wvvvv  = Gvvvv  #G(abef)
    Wvvvv -= np.einsum('mb,amef->abef', T1, Gvovv)
    Wvvvv -= np.einsum('ma,bmfe->abef', T1, Gvovv) 
    Wvvvv += np.einsum('mnab,mnef->abef', tau, Goovv)

    if (EnvVal['HBAR_OUT']=='FILE'):
       util.write_data('Data-Wvvvv.txt',Wvvvv,EnvVal)
       return
    else:
       return Wvvvv


def get_Wovvo(EnvVal):
    # W(mbej) = <mb||ej> + t(jf)<mb||ef> - t(nb)<mn||ej> - (t(jnfb)+t(jf)t(nb))<mn||ef>
    T1=util.read_data('Data-T1.txt',EnvVal)
    T2=util.read_data('Data-T2.txt',EnvVal)
    tau=util.read_data('Data-tau.txt',EnvVal)
    Goovo=util.get_SubMO('oovo',EnvVal)
    Goovv=util.get_SubMO('oovv',EnvVal)
    Govvo=util.get_SubMO('ovvo',EnvVal)
    Govvv=util.get_SubMO('ovvv',EnvVal)
    Loovv=util.read_data('Int-Loovv.txt',EnvVal)

    Wovvo  = Govvo  #G(mbej)
    Wovvo += np.einsum('jf,mbef->mbej', T1, Govvv)
    Wovvo -= np.einsum('nb,mnej->mbej', T1, Goovo)
    Wovvo -= np.einsum('jnfb,mnef->mbej', tau, Goovv)
    Wovvo += np.einsum('njfb,mnef->mbej', T2, Loovv)

    if (EnvVal['HBAR_OUT']=='FILE'):
       util.write_data('Data-Wovvo.txt',Wovvo,EnvVal)
       return
    else:
       return Wovvo


def get_Wovov(EnvVal):   #Check
    # W(mbje) = -W(mbej)
    T1=util.read_data('Data-T1.txt',EnvVal)
    tau=util.read_data('Data-tau.txt',EnvVal)
    Gooov=util.get_SubMO('ooov',EnvVal)
    Goovv=util.get_SubMO('oovv',EnvVal)
    Govov=util.get_SubMO('ovov',EnvVal)
    Gvovv=util.get_SubMO('vovv',EnvVal)
    Govvv=util.get_SubMO('ovvv',EnvVal)

    Wovov  = Govov  #G(mbje)
    Wovov += np.einsum('jf,mbfe->mbje', T1, Govvv)
    Wovov -= np.einsum('nb,mnje->mbje', T1, Gooov)
    Wovov -= np.einsum('jnfb,mnfe->mbje', tau, Goovv)

    if (EnvVal['HBAR_OUT']=='FILE'):
       util.write_data('Data-Wovov.txt',Wovov,EnvVal)
       return
    else:
       return Wovov


def get_Woovv(EnvVal):   #Check
    Loovv=util.read_data('Int-Loovv.txt',EnvVal)
    
    Woovv=Loovv

    if (EnvVal['HBAR_OUT']=='FILE'):
       util.write_data('Data-Woovv.txt',Woovv,EnvVal)
       return
    else:
       return Woovv

 
def get_Wooov(EnvVal):
    T1=util.read_data('Data-T1.txt',EnvVal)
    Gooov=util.get_SubMO('ooov',EnvVal)
    Goovv=util.get_SubMO('oovv',EnvVal)

    Wooov  = Gooov  #G(mnie)
    Wooov += np.einsum('if,mnfe->mnie', T1, Goovv)

    if (EnvVal['HBAR_OUT']=='FILE'):
       util.write_data('Data-Wooov.txt',Wooov,EnvVal)
       return
    else:
       return Wooov


def get_Wvovv(EnvVal):
    T1=util.read_data('Data-T1.txt',EnvVal)
    Goovv=util.get_SubMO('oovv',EnvVal)
    Gvovv=util.get_SubMO('vovv',EnvVal)

    Wvovv  = Gvovv  #G(amef)
    Wvovv -= np.einsum('na,nmef->amef', T1, Goovv)

    if (EnvVal['HBAR_OUT']=='FILE'):
       util.write_data('Data-Wvovv.txt',Wvovv,EnvVal)
       return
    else:
       return Wvovv

def get_Wovoo(EnvVal):  
    # need Fov, Woooo
    T1=util.read_data('Data-T1.txt',EnvVal)
    T2=util.read_data('Data-T2.txt',EnvVal)
    tau=util.read_data('Data-tau.txt',EnvVal)
    Gooov=util.get_SubMO('ooov',EnvVal)
    Goovo=util.get_SubMO('oovo',EnvVal)
    Goovv=util.get_SubMO('oovv',EnvVal)
    Govoo=util.get_SubMO('ovoo',EnvVal)
    Govov=util.get_SubMO('ovov',EnvVal)
    Govvo=util.get_SubMO('ovvo',EnvVal)
    Govvv=util.get_SubMO('ovvv',EnvVal)
    Looov=util.read_data('Int-Looov.txt',EnvVal)
    Loovv=util.read_data('Int-Loovv.txt',EnvVal)
    Fov=util.read_data('Data-Fov.txt',EnvVal)
    Woooo=util.read_data('Data-Woooo.txt',EnvVal)

    Wovoo  = Govoo  #G(mbij)
    Wovoo += np.einsum('ijeb,me->mbij', T2, Fov)
    Wovoo -= np.einsum('nb,mnij->mbij', T1, Woooo)
    Wovoo += np.einsum('ijef,mbef->mbij', tau, Govvv)

    Wovoo += np.einsum('jnbe,mnie->mbij', T2, Looov)  #P(ij)t(jnbe)<mn||ie>
    Wovoo -= np.einsum('njbe,mnie->mbij', T2, Gooov)  
    Wovoo -= np.einsum('ineb,mnej->mbij', T2, Goovo)

    Wovoo += np.einsum('ie,mbej->mbij',T1, Govvo)     #P(ij)t(ie)<mb||ej>
    Wovoo += np.einsum('je,mbie->mbij',T1, Govov)

    tmp    = np.einsum('njbf,mnef->mbej',T2, Goovv)   #P(ij)t(ie){...}
    tmp   -= np.einsum('jnbf,mnef->mbej',T2, Loovv)
    Wovoo -= np.einsum('ie,mbej->mbij',T1, tmp)
    tmp    = np.einsum('infb,mnfe->mbie',T2, Goovv)
    Wovoo -= np.einsum('je,mbie->mbij',T1, tmp)

    if (EnvVal['HBAR_OUT']=='FILE'):
       util.write_data('Data-Wovoo.txt',Wovoo,EnvVal)
       return
    else:
       return Wovoo


def get_Wvvvo(EnvVal):  #FIXME
    T1=util.read_data('Data-T1.txt',EnvVal)
    T2=util.read_data('Data-T2.txt',EnvVal)
    tau=util.read_data('Data-tau.txt',EnvVal)
    Goovo=util.get_SubMO('oovo',EnvVal)
    Goovv=util.get_SubMO('oovv',EnvVal)
    Govvo=util.get_SubMO('ovvo',EnvVal)
    Govvv=util.get_SubMO('ovvv',EnvVal)
    Gvovo=util.get_SubMO('vovo',EnvVal)
    Gvovv=util.get_SubMO('vovv',EnvVal)
    Gvvvo=util.get_SubMO('vvvo',EnvVal)
    Loovv=util.read_data('Int-Loovv.txt',EnvVal)
    Lvovv=util.read_data('Int-Lvovv.txt',EnvVal)

    Fov=util.read_data('Data-Fov.txt',EnvVal)
    Wvvvv=util.read_data('Data-Wvvvv.txt',EnvVal)

    Wvvvo  = Gvvvo  #G(abei)
    Wvvvo -= np.einsum('miab,me->abei',T2, Fov)       #t(miab)F(me)
    Wvvvo += np.einsum('if,abef->abei',T1, Wvvvv)     #t(if)W(abef)
    Wvvvo += np.einsum('mnab,mnei->abei',tau, Goovo)  #tau(mnab)<mn|ei>
                                                      
    Wvvvo -= np.einsum('miaf,mbef->abei',T2, Govvv)   #P(ab)t(miaf)<mb||ef>
    Wvvvo -= np.einsum('mibf,amef->abei',T2, Gvovv) 
    Wvvvo += np.einsum('imbf,amef->abei',T2, Lvovv) 

    Wvvvo -= np.einsum('ma,mbei->abei',T1, Govvo)     #P(ab)t(ma)<mb||ei>
    Wvvvo -= np.einsum('mb,amei->abei',T1, Gvovo)

    tmp    = np.einsum('nibf,mnef->mbei',T2, Goovv)   #P(ab)t(ma)t(nibf)<mn||ef>
    tmp   -= np.einsum('inbf,mnef->mbei',T2, Loovv)
    Wvvvo += np.einsum('ma,mbei->abei',T1, tmp)
    tmp    = np.einsum('niaf,nmef->amei',T2, Goovv)
    Wvvvo += np.einsum('mb,amei->abei',T1, tmp)

    if (EnvVal['HBAR_OUT']=='FILE'):
       util.write_data('Data-Wvvvo.txt',Wvvvo,EnvVal)
       return
    else:
       return Wvvvo


if __name__ == '__main__' :

    util.make_header('Hbar program')
    argv = sys.argv
    EnvVal=inp.driver(argv)

    EnvVal['HBAR_TYPE']='CALC'
#   EnvVal['HBAR_TYPE']='INTS'
#   EnvVal['HBAR_TYPE']='TEST'

    if  EnvVal['HBAR_TYPE']=='TEST':
        test_results(EnvVal)
    else:
        driver(EnvVal)



