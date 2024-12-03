import sys
import numpy as np
import Base_Util as util

def read_dgetrecpy(string,dim1):
    import aces2py as a2
    Vec=np.zeros([dim1]) 
    a2.dgetrecpy(string,Vec,dim1)
    return Vec

def read_getlstpy(idx,ispin,dim1,dim2):
    import aces2py as a2
    Irrep=1
    dim=dim1*dim2
    Vec=np.zeros(dim)
    a2.getlstpy(1,1,Irrep,ispin,idx,dim,Vec)
    Vec=Vec.reshape(dim1,dim2)
    return Vec

def read_getallpy(idx,Dim):
    import aces2py as a2
    Xdim=len(Dim)
    Irrep=1

    Ndim=1
    for i in range (Xdim):
        Ndim=Dim[i]
    if Ndim==0:
       write(' * Error: in reading getallpy, idx='+str(idx)) 
       write('          Ndim = 0')

    Vec=np.zeros(Ndim)
    a2.getallpy(Vec,Ndim,Irrep,idx)
    if Xdim==2:
       Vec=Vec.reshape(Dim[0],Dim[1])
    elif Xdim==4:
       Vec=Vec.reshape(Dim[0],Dim[1],Dim[2],Dim[3])

    return Vec

def read_get12ehpy(string,Nbas,Ndim):
    import aces2py as a2
    
    Mat=np.zeros(Nbas**Ndim, dtype=float) 
    if (Ndim==2):
       a2.get1ehpy(string,Mat,Nbas) 
       Mat=Mat.reshape(Nbas,Nbas)
    elif (Ndim==4):
       a2.get2ehpy(string,Mat,Nbas) 
       Mat=Mat.reshape(Nbas,Nbas,Nbas,Nbas)
    #print('\n * matrix('+string+')')
    #print(Mat)
    return Mat

def read_mol_info(EnvVal):
    import aces2py as a2

    Irrep=1
    Nbas = np.array([0])  # Number of basis functions
    Nocc = np.array([0])  # Number of electron pairs
    charge = np.array([0])  # Charge of Species

    a2.igetrecpy('NBASTOT', Nbas, 1)
    a2.igetrecpy('NMPROTON', Nocc, 1)
    Enuc=read_dgetrecpy('NUCREP',1) # Nuclear Repulsion Energy
    charge = a2.flags.iflags[27]

    Nbas = Nbas[0]
    print('- No. basis function = '+str(Nbas))
    Nocc = int((Nocc[0] - charge) / 2)
    print('- No, occupied orbital = '+str(Nocc))
    print('- Nuclear repulsion energy = '+str(Enuc[0]))

    Nvrt=Nbas-Nocc

    EnvVal['NIRREP']=Irrep
    EnvVal['NBAS']=Nbas
    EnvVal['NOCC']=Nocc
    EnvVal['NVRT']=Nvrt

    return EnvVal


def read_ints(Nbas):
     
    Smat = read_get12ehpy('overlap',Nbas,2)  # overlap Matrix
    Int1 = read_get12ehpy('oneh',Nbas,2)     # H core Matrix
    Int2 = read_get12ehpy('2elints',Nbas,4)  # 2e integral 

#   Int2 = Int2.swapaxes(2, 3)
#   Jmat = Int2.reshape(Nbas**2,Nbas**2)
    Jmat = Int2.swapaxes(2, 3)
    Jmat = Jmat.reshape(Nbas**2,Nbas**2)
    Kmat = Int2.swapaxes(1, 2)  
    Kmat = Kmat.reshape(Nbas**2,Nbas**2)

    return Smat,Int1,Jmat,Kmat

def write_aces2_ints_file():
    import aces2py as a2
    f=a2.init()
    a2.run('ints')
    #a2.run('scf')
    #a2.run('1props')
    #a2.run('cc')
    #a2.run('ee')
    #a2.run('hbar')

    return

def get_Vec2D(Nocc,Nvrt,String,EnvVal,lprint):
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

    String=String+'.txt'
    util.write_data(String,Mat,EnvVal)

    return Mat


def get_Vec4D(Nocc,Nvrt,String,EnvVal,lprint):
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

#   if (String=='Tvvoo'):

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

    String=String+'.txt'
    util.write_data(String,Mat,EnvVal)

    return Mat


def get_T(Nocc,Nvrt,String,lprint):
    import aces2py as a2
    print(String)

def get_F(Nocc,Nvrt,String,lprint):
    import aces2py as a2
    print(String)

def get_2dVec(String,iList,idx,Nocc,Nvrt):
#   -------------------------------------
#   Get the size of raw and final tensors
#   -------------------------------------
    print('\n * '+String+' : '+str(iList)+' -> '+idx)
    dim=np.zeros(2)
    idx1=idx[:2]
    idx2=idx[2:4]

    if idx1=='oo': dim[0]=Nocc; dim[1]=Nocc; RawDim1=Nocc*Nocc
    if idx1=='vv': dim[0]=Nvrt; dim[1]=Nvrt; RawDim1=Nvrt*Nvrt
    if idx1=='ov': dim[0]=Nocc; dim[1]=Nvrt; RawDim1=Nocc*Nvrt
    if idx1=='vo': dim[0]=Nvrt; dim[1]=Nocc; RawDim1=Nvrt*Nocc

    if idx2=='aa': ispin=1
    if idx2=='bb': ispin=2

    dim=dim.astype(int)
    RawDim=int(RawDim1)
#   print('RawDim = '+str(RawDim))

#   ------------------
#   Read F from AcesPy
#   ------------------
    Irrep=1
    RawMat=np.zeros(RawDim)
    a2.getlstpy(1,1,Irrep,ispin,iList,RawMat,RawDim)
    print('* Read '+String+' from AcesPy(getlstpy). RawDim='+str(RawDim))
    Mat=RawMat.reshape(dim[0],dim[1])
    return Mat,dim


def get_CIS_vector(Nocc,Nvrt):
    import aces2py as a2
    Irrep=1
    Nov=Nocc*Nvrt
    Mat=read_getlstpy(94,1,Nov,Nov)
    util.print_mat2('CIS vector',Mat,Nov,Nov)



