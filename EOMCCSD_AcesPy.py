import sys
import numpy as np
import EOMCCSD_Util as util


def read_mol_info(EnvVal):
    import aces2py as a2

    Irrep=1
    Nbas = np.array([0])  # Number of basis functions
    Nocc = np.array([0])  # Number of electron pairs
    Enuc = np.zeros([1])  # Nuclear Repulsion Energy
    charge = np.array([0])  # Charge of Species

    a2.igetrecpy('NBASTOT', Nbas, 1)
    a2.igetrecpy('NMPROTON', Nocc, 1)
    a2.dgetrecpy('NUCREP', Enuc, 1)
    charge = a2.flags.iflags[27]

    Nbas = Nbas[0]
    print('- No. basis function = '+str(Nbas))
    Nocc = int((Nocc[0] - charge) / 2)
    print('- No, occupied orbital = '+str(Nocc))
    Enuc = Enuc[0]
    print('- Nuclear repulsion energy = '+str(Enuc))

    Nvrt=Nbas-Nocc

    EnvVal['Irrep']=Irrep
    EnvVal['Nbas']=Nbas
    EnvVal['Nocc']=Nocc
    EnvVal['Nvrt']=Nvrt

    return EnvVal


def read_ints(Nbas):
    import aces2py as a2
    Smat = np.zeros(Nbas**2, dtype=float)  # Overlap Matrix
    a2.get1ehpy('overlap', Smat, Nbas)  # Overlap Matrix
    Smat = Smat.reshape(Nbas,Nbas)
    print ('\n * S matrix')
    print (Smat)
    #print(Smat.nbytes)

    Int1 = np.zeros(Nbas**2)  # One-Electron Integrals
    a2.get1ehpy('oneh',Int1,Nbas)  # One-Electron Integrals
    Int1 = Int1.reshape(Nbas,Nbas)
    print ('\n *H core matrix')
    print (Int1)

    Int2 = np.zeros(Nbas**4)  # Two-Electron Integrals
    a2.get2ehpy('2elints',Int2,Nbas)  # Two-Electron Integrals
    Int2 = Int2.reshape(Nbas,Nbas,Nbas,Nbas)
#   print ('H rep matrix')
#   print (Int2)

#   Int2 = Int2.swapaxes(2, 3)
#   Jmat = Int2.reshape(Nbas**2,Nbas**2)
    Jmat = Int2.swapaxes(2, 3)
    Jmat = Jmat.reshape(Nbas**2,Nbas**2)
    Kmat = Int2.swapaxes(1, 2)  # Why not swapaxes(1,3)?
    Kmat = Kmat.reshape(Nbas**2,Nbas**2)

    # np.savetxt('Int1.csv', Int1, delimiter=',')
    # np.savetxt('Int2.csv', Int2, delimiter=',')
    # np.savetxt('Jmat.csv', Jmat, delimiter=',')
    # np.savetxt('Kmat.csv', Kmat, delimiter=',')

    return Smat,Int1,Jmat,Kmat



def get_Vec2D(Nocc,Nvrt,String,lprint):
    import aces2py as a2
    # See checkhbar.F
    Irrep=1

    dim=[];Ndim=1;ind=0
    for i in range(1,3):
        if (String[i]=='o'):
           dim.append(int(Nocc))
           Ndim=Ndim*Nocc
        elif (String[i]=='v'):
           dim.append(int(Nvrt))
           Ndim=Ndim*Nvrt
        else:
           print('read_Hbar1: Error!!!')

    if (String=='Foo_aa'):  ind=91; ispin=1 #[91] Foo(MI)
    if (String=='Foo_bb'):  ind=91; ispin=2 #[91] Foo(MI)
    if (String=='Fvv_aa'):  ind=92; ispin=1 #[92] Fvv(EA)
    if (String=='Fvv_bb'):  ind=92; ispin=2 #[92] Fvv(EA)
    if (String=='Fov_aa'):  ind=93; ispin=1 #[93] Fov(ME)
    if (String=='Fov_bb'):  ind=93; ispin=2 #[93] Fov(ME)

    if (ind==0):
       print('Error, String is not right!')

    Mat=np.zeros(Ndim)
    a2.getlstpy(1,1,Irrep,ispin,ind,Ndim,Mat)
    if (String=='Fov_aa'): 
       Mat=Mat.reshape(dim[0],dim[1])
    else:
       Mat=Mat.reshape(dim[0],dim[1])
#      Mat=Mat.swapaxes(0,1)

    if (lprint):
       print(String+', Ndim='+str(Ndim))
       print(dim)
       print(Mat)

    lwrite=True
    util.write_data(String,Mat,dim,lwrite)

    return Mat


def get_Vec4D(Nocc,Nvrt,String,lprint):
    import aces2py as a2
    # See checkhbar.F
    Irrep=1

    dim=[];Ndim=1;ind=0
    for i in range(1,5):
        if (String[i]=='o'):
           dim.append(int(Nocc))
           Ndim=Ndim*Nocc
        elif (String[i]=='v'):
           dim.append(int(Nvrt))
           Ndim=Ndim*Nvrt
        else:
           print('read_Vec4D: Error!!!')

    # H(MN,IJ) 0
    if (String=='Woooo_aaaa'):  ind=51;   # W(MN,IJ) !  M<N I<J
    if (String=='Woooo_bbbb'):  ind=52;   # W(MN,IJ) !  m<n i<j
    if (String=='Woooo_abab'):  ind=53;   # W(MN,IJ) !  M,n I,j

    # H(AI,BC) 3
    if (String=='Wvovv_aaaa'):  ind=27;   # W(AI,BC) !  A<B C,I  
    if (String=='Wvovv_bbbb'):  ind=28;   # W(AI,BC) !  a<b c,i
    if (String=='Wvovv_aabb'):  ind=29;   # W(AI,BC) !  A,b I,c *
    if (String=='Wvovv_abba'):  ind=30;   # W(AI,BC) !  A,b C,i

    # H(IJ,KA) 1
    if (String=='Wooov_aaaa'):  ind=7;    # W(IJ,KA) !  I<J K,A 
    if (String=='Wooov_bbbb'):  ind=8;    # W(ij,ka) !  i<j k,a
    if (String=='Wooov_abba'):  ind=9;    # W(Ij,kA) !  I,j A,k *
    if (String=='Wooov_abab'):  ind=10;   # W(Ij,Ka) !  I,j K,a 

    # H(MB,EJ) 2
    if (String=='Wovov_aaaa'):  ind=54;   # W(MB,JE) !  E,M B,j
    if (String=='Wovov_bbbb'):  ind=55;   # W(MB,JE) !  e,m b,j 
    if (String=='Wovov_aabb'):  ind=56;   # W(MB,JE) !  E,M b,j
    if (String=='Wovov_bbaa'):  ind=57;   # W(MB,JE) !  e,m B,J
    if (String=='Wovov_abab'):  ind=58;   # W(MB,JE) !  E,m B,j
    if (String=='Wovov_baba'):  ind=59;   # W(MB,JE) !  e,M b,J

    # H(AB,CI) 3
    if (String=='Wvvvo_aaaa'):  ind=127;  # W(AB,CI) !  A<B C,I 
    if (String=='Wvvvo_bbbb'):  ind=128;  # W(AB,CI) !  a<b c,i
    if (String=='Wvvvo_abba'):  ind=129;  # W(AB,CI) !  A,b I,c *
    if (String=='Wvvvo_abba'):  ind=130;  # W(AB,CI) !  A,b C,i

    # H(IA,JK) 1
    if (String=='Wovoo_aaaa'):  ind=107;  # W(IA,JK) !  I<J K,A 
    if (String=='Wovoo_bbbb'):  ind=108;  # W(IA,JK) !  i<j k,a
    if (String=='Wovoo_aabb'):  ind=109;  # W(IA,JK) !  I,j A,k *
    if (String=='Wovoo_abba'):  ind=110;  # W(IA,JK) !  I,j K,a

    # H(AB,CD) 4
    if (String=='Wvvvv_aaaa'):  ind=231;  # W(AB,CD) !  A<B C<D 
    if (String=='Wvvvv_bbbb'):  ind=232;  # W(AB,CD) !  a<b c<d
    if (String=='Wvvvv_abab'):  ind=233;  # W(AB,CD) !  A,b C,d

    # T2 (34-37 or 44-46)? 
    if (String=='Tovov_aaaa'):  ind=34;   # T 
    if (String=='Tovov_bbbb'):  ind=35;   
    if (String=='Tovov_aabb'):  ind=36;  
    if (String=='Tovov_bbaa'):  ind=37;  

    if (ind==0):
       print('Error, String is not right! '+String)

    Mat=np.zeros(Ndim)
    print(String)
    a2.getallpy(Mat,Ndim,Irrep,ind)
    Mat=Mat.reshape(dim[0],dim[1],dim[2],dim[3])

    if (lprint):
       print(String+', Ndim='+str(Ndim))
       print(dim)
       #print(Mat)

    lwrite=True
    util.write_data(String,Mat,dim,lwrite)

    return Mat


def get_Tamp(Nocc,Nvrt,String,lprint):
    import aces2py as a2

    print(String)


def get_CIS_vector(Nocc,Nvrt):
    import aces2py as a2
    Irrep=1
    Nov=Nocc*Nvrt
    Ndim=Nov*Nov
    dim=[Nov,Nov]
    Mat=np.zeros(Ndim)
    a2.getlstpy(1,1,Irrep,1,94,Ndim,Mat)
    Mat=Mat.reshape(dim[0],dim[1])

    lprint=True
    String='CIS vector'
    if (lprint):
       print(String+', Ndim='+str(Ndim))
       print(dim)
       print(Mat)


def get_CIS_matrix(Nocc,Nvrt):
    import aces2py as a2
    # See makess.f

    lprint=False
    Irrep=1

    Nbas=Nocc+Nvrt
    Nov=Nocc*Nvrt
    OrbEng=np.zeros(Nbas)
    a2.dgetrecpy('SCFEVALA',OrbEng,Nbas)
    print(OrbEng)

    if (False):
       # True when NON-HF reference is used.
       #Foo
       Ndim=Nocc*Nocc
       dim=[Nocc,Nocc]
       Foo=np.zeros(Ndim)
       a2.getlstpy(1,1,Irrep,1,91,Ndim,Foo)
       Foo=Foo.reshape(dim[0],dim[1])
       String='Foo matrix'
       lprint=True
       if (lprint):
          print(String+', Ndim='+str(Ndim))
          print(dim)
          print(Foo)

       #Fvv
       Ndim=Nvrt*Nvrt
       dim=[Nvrt,Nvrt]
       Fvv=np.zeros(Ndim)
       a2.getlstpy(1,1,Irrep,1,92,Ndim,Fvv)
       Fvv=Fvv.reshape(dim[0],dim[1])
       String='Fvv matrix'
       lprint=True
       if (lprint):
          print(String+', Ndim='+str(Ndim))
          print(dim)
          print(Fvv)

    #<IA||JB>
    Ndim=Nov*Nov
    dim=[Nov,Nov]
    Wovov_aaaa=np.zeros(Ndim)
    a2.getallpy(Wovov_aaaa,Ndim,Irrep,23)
    Wovov_aaaa=Wovov_aaaa.reshape(dim[0],dim[1])
    String='TEI matrix (aaaa)'
    if (lprint):
       print(String+', Ndim='+str(Ndim))
       print(dim)
#      print(Wovov_aaaa)
       print(Wovov_aaaa[0,:])

    #<Ab||Ij>
    Ndim=Nov*Nov
    dim=[Nov,Nov]
    Wovov_abab=np.zeros(Ndim)
    a2.getallpy(Wovov_abab,Ndim,Irrep,18)
    Wovov_abab=Wovov_abab.reshape(dim[0],dim[1])
    String='TEI matrix (abab)'
    if (lprint):
       print(String+', Ndim='+str(Ndim))
       print(dim)
#      print(Wovov_abab)
       print(Wovov_abab[0,:])


    Wovov = -Wovov_aaaa + Wovov_abab
    String='TEI matrix (aaaa+abab)'
    if (lprint):
       print(String+', Ndim='+str(Ndim))
       print(dim)
#      print(Wovov_abab)
       print(Wovov[0,:])

    # H(ai,bj) = {F(ab)-F(ij)}d(i==j)d(a==b) - W(ai|bj)
    Mat = -Wovov_aaaa + Wovov_abab
    ia=0
    for i in range(Nocc):
        for a in range(Nvrt):
            ia +=1
            jb =0
            for j in range(Nocc):
                for b in range(Nvrt):
                    jb +=1
                    if (i==j and a==b):
                       Mat[ia-1,jb-1] += OrbEng[a+Nocc]-OrbEng[i]

    String='CIS matrix'
    if (lprint):
       print(String+', Ndim='+str(Ndim))
       print(dim)
       print(Mat)

    lwrite=True
    util.write_data('CIS-Matrix',Mat,dim,lwrite)

    RunDiag=True
    if (RunDiag):
       EigVal,EigVec = np.linalg.eig(Mat)
       util.write_data('CIS-EigVec',EigVec,dim,lwrite)
       dim=[Nov]
       util.write_data('CIS-EigVal',EigVal,dim,lwrite)
       dim=[Nbas]
       util.write_data('SCFEVALA',OrbEng,dim,lwrite)
       print('EigVal')
       print(EigVal)

    return EigVal,EigVec



