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
    F['oo_aa']=ap.get_Vec2D(Nocc,Nvrt,'Foo_aa',lprint)
    F['oo_bb']=ap.get_Vec2D(Nocc,Nvrt,'Foo_bb',lprint)
    F['vv_aa']=ap.get_Vec2D(Nocc,Nvrt,'Fvv_aa',lprint)
    F['vv_bb']=ap.get_Vec2D(Nocc,Nvrt,'Fvv_bb',lprint)
    F['ov_aa']=ap.get_Vec2D(Nocc,Nvrt,'Fov_aa',lprint)
    F['ov_bb']=ap.get_Vec2D(Nocc,Nvrt,'Fov_bb',lprint)

    W={}
    # H(MN,IJ) OOOO
    W['oooo_aaaa']=ap.get_Vec4D(Nocc,Nvrt,'Woooo_aaaa',lprint) #[51] (M<N I<J)
    W['oooo_bbbb']=ap.get_Vec4D(Nocc,Nvrt,'Woooo_bbbb',lprint) #[52] (m<n i<j)
    W['oooo_abab']=ap.get_Vec4D(Nocc,Nvrt,'Woooo_abab',lprint) #[53] (M,n I,j)
                                                     
    # H(AI,BC) VOVV                         
    W['vovv_aaaa']=ap.get_Vec4D(Nocc,Nvrt,'Wvovv_aaaa',lprint) #[27] (A<B C,I)
    W['vovv_bbbb']=ap.get_Vec4D(Nocc,Nvrt,'Wvovv_bbbb',lprint) #[28] (a<b c,i)
    W['vovv_abab']=ap.get_Vec4D(Nocc,Nvrt,'Wovvv_abab',lprint) #[29] (A,b I,c)
    W['vovv_abba']=ap.get_Vec4D(Nocc,Nvrt,'Wvovv_abab',lprint) #[30] (A,b C,i)
                                                     
    # H(IJ,KA)                         
    W['vooo_aaaa']=ap.get_Vec4D(Nocc,Nvrt,'Wvooo_aaaa',lprint) #
    W['vooo_bbbb']=ap.get_Vec4D(Nocc,Nvrt,'Wvooo_bbbb',lprint)
    W['vooo_abab']=ap.get_Vec4D(Nocc,Nvrt,'Wvooo_abab',lprint)
    W['vooo_abba']=ap.get_Vec4D(Nocc,Nvrt,'Wvooo_abba',lprint)
                                                     
    # H(MB,EJ)                        
    W['ovvo_aaaa']=ap.get_Vec4D(Nocc,Nvrt,'Wovvo_aaaa',lprint)
    W['ovvo_bbbb']=ap.get_Vec4D(Nocc,Nvrt,'Wovvo_bbbb',lprint)
    W['ovvo_aabb']=ap.get_Vec4D(Nocc,Nvrt,'Wovvo_aabb',lprint)
    W['ovvo_bbaa']=ap.get_Vec4D(Nocc,Nvrt,'Wovvo_bbaa',lprint)
    W['ovvo_abab']=ap.get_Vec4D(Nocc,Nvrt,'Wovvo_abab',lprint)
    W['ovvo_abba']=ap.get_Vec4D(Nocc,Nvrt,'Wovvo_baba',lprint)
                                                     
    # H(AB,CI)                       
    W['ovvv_aaaa']=ap.get_Vec4D(Nocc,Nvrt,'Wovvv_aaaa',lprint)
    W['ovvv_bbbb']=ap.get_Vec4D(Nocc,Nvrt,'Wovvv_bbbb',lprint)
    W['ovvv_abab']=ap.get_Vec4D(Nocc,Nvrt,'Wovvv_abab',lprint)
    W['ovvv_abba']=ap.get_Vec4D(Nocc,Nvrt,'Wovvv_abba',lprint)
                                                     
    # H(IA,JK)                      
    W['oovo_aaaa']=ap.get_Vec4D(Nocc,Nvrt,'Woovo_aaaa',lprint)
    W['oovo_bbbb']=ap.get_Vec4D(Nocc,Nvrt,'Woovo_bbbb',lprint)
    W['oovo_abab']=ap.get_Vec4D(Nocc,Nvrt,'Woovo_abab',lprint)
    W['oovo_abba']=ap.get_Vec4D(Nocc,Nvrt,'Woovo_abba',lprint)
                                                     
    # H(AB,CD)                     
    W['vvvv_aaaa']=ap.get_Vec4D(Nocc,Nvrt,'Wvvvv_aaaa',lprint)
    W['vvvv_bbbb']=ap.get_Vec4D(Nocc,Nvrt,'Wvvvv_bbbb',lprint)
    W['vvvv_abab']=ap.get_Vec4D(Nocc,Nvrt,'Wvvvv_abab',lprint)

    # T2                     
    T={}
    T['ovov_aaaa']=ap.get_Vec4D(Nocc,Nvrt,'Tovov_aaaa',lprint)
    T['ovov_bbbb']=ap.get_Vec4D(Nocc,Nvrt,'Tovov_bbbb',lprint)
    T['ovov_aabb']=ap.get_Vec4D(Nocc,Nvrt,'Tovov_aabb',lprint)
    T['ovov_bbaa']=ap.get_Vec4D(Nocc,Nvrt,'Tovov_bbaa',lprint)

    

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
    tau  = T2 
    tau += np.einsum('ia,jb->ijab',T1,T1)
    return

def get_Fvv():
    # F(ae) = f(ae) - t(am)f(me) + t(fm)<am||ef> - 1/2 tau(afmn)<mn|ef> 
    # F(ab) = f(ab) - t(ma)f(mb) + t(mf)<am||ef> - tau(mnfa)<mn||fe>    
    Fvv_aa  = fvv_a 
    Fvv_aa -= np.einsum("am,me->ae",     T_a, fov_a)
    Fvv_aa += np.einsum("fm,amef->ae",   T_a, Gvovv_aa)
    Fvv_aa -= np.einsum("afmn,mnef->ae", Tau_aa, Goovv_aa) * 0.5
#   Fvv_aa += np.einsum("ambe,em->ab",   Gvovv_ab, T_b)
#   Fvv_aa -= np.einsum("mnbf,afmn->ab", Goovv_ab, T_ab)

def get_Foo():
    # F(mi) = f(mi) + t(ei)f(me) + t(en)<mn||ie> + 1/2 tau(efin)<mn|ef> 
    Foo_aa  = foo_a
    Foo_aa += np.einsum('ei,me->mi',     T_a, fov_a)
    Foo_aa += np.einsum('en,mnie->mi',   T_a, Gooov_aa)
    Foo_aa += np.einsum('efin,mnef->mi', Tau_aa, Goovv_aa)
    return Foo_aa

def get_Fov(f,G,T):
    fov_a = f['ov_a']
    Goovv_aa = G['oovv_aa']
    T_a = T['a']

    # F(me) = f(me) + t(fn)<mn||ef>
    Fov_aa  = fov_a 
    Fov_aa += np.einsum('fn,mnef->me',   T_a, Goovv_aa) 
    return Fov_aa

def make_Hbar():
    
    Foo_aa  = np.einsum("je,ei->ji",     H.a.ov,     T.a)
    Foo_aa =+ np.einsum("jmie,em->ji",   H0.aa.ooov, T.a)
    Foo_aa =+ np.einsum("jmie,em->ji",   H0.ab.ooov, T.b)
    Foo_aa =+ np.einsum("jnef,efin->ji", H0.aa.oovv, T.aa) * 0.5
    Foo_aa =+ np.einsum("jnef,efin->ji", H0.ab.oovv, T.ab)





