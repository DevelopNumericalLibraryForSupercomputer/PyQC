import os
import numpy as np

def make_header(string):
    nstr=len(string)+4
    print('\n'+(' '*3) + ('='*nstr))
    print((' '*5) + string + '  ')
    print((' '*3) + ('='*nstr) +'\n')
    return


def write_data(String,mat,EnvVal):

    if not os.path.exists("Data"):
       os.makedirs("Data")
    DataDir=EnvVal['DATA_DIR']
    fout=DataDir+'/'+String
    print("\n - "+String+" is written in "+DataDir)
    if os.path.exists(fout):
       os.remove(fout)

    dim=mat.shape
    f=open(fout,'w')
    if (len(dim)==1):
       f.write("%6d\n" % (dim[0]))
       for i0 in range(dim[0]):
           f.write("%6d  %25.15f\n" % (i0,mat[i0]))
    elif (len(dim)==2):
       f.write("%6d %6d\n" % (dim[0],dim[1]))
       for i0 in range(dim[0]):
           for i1 in range(dim[1]):
               f.write("%6d %6d  %25.15f\n" % (i0,i1,mat[i0,i1]))
    elif (len(dim)==4):
       f.write("%6d %6d %6d %6d\n" % (dim[0],dim[1],dim[2],dim[3]))
       for i0 in range(dim[0]):
           for i1 in range(dim[1]):
               for i2 in range(dim[2]):
                   for i3 in range(dim[3]):
                       f.write("%6d %6d %6d %6d  %25.15f\n" % (i0,i1,i2,i3,mat[i0,i1,i2,i3]))
    f.close()

def read_value(String,EnvVal):
    DataDir=EnvVal['DATA_DIR']
    fout=DataDir+'/'+String
    f=open(fout,'r')
    val=f.readline().strip()
    f.close()
    return val

def read_data(String,EnvVal):
    DataDir=EnvVal['DATA_DIR']
    LDebug=EnvVal['HBAR_DEBUG']
    fout=DataDir+'/'+String
    #print(' - '+fout)

    f=open(fout,'r')
    tmp=f.readline().split()
    dim=[]
    for val in tmp:
        dim.append(int(val))

    if LDebug=='TRUE':
       idx=len(dim)
    else:
       idx=0

    if (len(dim)==1):
       if LDebug=='TRUE':
          idx=0
       else:
          idx=1
       mat=np.zeros(dim[0])
       for i0 in range(dim[0]):
           mat[i0]=float(f.readline().split()[idx])

    elif (len(dim)==2):
       mat=np.zeros((dim[0],dim[1]))
       for i0 in range(dim[0]):
           for i1 in range(dim[1]):
               mat[i0,i1]=float(f.readline().split()[idx])

    elif (len(dim)==4):
       mat=np.zeros((dim[0],dim[1],dim[2],dim[3]))
       for i0 in range(dim[0]):
           for i1 in range(dim[1]):
               for i2 in range(dim[2]):
                   for i3 in range(dim[3]):
                       mat[i0,i1,i2,i3]=float(f.readline().split()[idx])
    f.close()
    return mat


def read_data2(String,lread):
    DataDir=EnvVal['DATA_DIR']
    fout=DataDir+'/Data-'+String+'.txt'
    #print(' - '+fout)

    f=open(fout,'r')
    tmp=f.readline().split()
    dim=[]
    for val in tmp:
        dim.append(int(val))
    if (len(dim)==1):
       mat=np.zeros(dim[0])
       for i0 in range(dim[0]):
           mat[i0]=float(f.readline().split()[1])

    elif (len(dim)==2):
       mat=np.zeros((dim[0],dim[1]))
       for i0 in range(dim[0]):
           for i1 in range(dim[1]):
               mat[i0,i1]=float(f.readline().split()[2])

    elif (len(dim)==4):
       mat=np.zeros((dim[0],dim[1],dim[2],dim[3]))
       for i3 in range(dim[3]):
           for i2 in range(dim[2]):
               for i1 in range(dim[1]):
                   for i0 in range(dim[0]):
                       mat[i0,i1,i2,i3]=float(f.readline().split()[4])
    f.close()
    return mat


def mat_swap(Mat,iswap,isign):
    if isign ==1:
       NewMat=np.transpose(Mat,iswap)
    else:
       #print('sign changed')
       NewMat=isign*np.transpose(Mat,iswap)
    #print('NewMat')
    #print(NewMat)
    return NewMat

def check_sum(String,Mat,dim):
    val1=0.0
    val2=0.0
#   print(Mat.shape)
    if (dim==1):
       for i0 in range(len(Mat)):
           val1+=Mat[i0]
           val2+=Mat[i0]*Mat[i0]
    elif (dim==2):
       for i0 in range(len(Mat)):
           for i1 in range(len(Mat[0])):
               val1+=Mat[i0,i1]
               val2+=Mat[i0,i1]*Mat[i0,i1]
    elif (dim==4):
       for i0 in range(len(Mat)):
           for i1 in range(len(Mat[0])):
               for i2 in range(len(Mat[0][0])):
                   for i3 in range(len(Mat[0][0][0])):
                       val1+=Mat[i0,i1,i2,i3]
                       val2+=Mat[i0,i1,i2,i3]*Mat[i0,i1,i2,i3]
    print(' @ Checksum: %s = %20.12f , %20.12f' % (String,val1,val2))
    #print(np.sum(Mat))

    return val2

def update_R2mod(R2,Nvrt,Nocc,string):
    NewR2=np.zeros([Nvrt,Nvrt,Nocc,Nocc])
    if string=='FF':
       for a in range(Nvrt):
           for b in range(Nvrt):
               for i in range(Nocc):
                   for j in range(Nocc):
                       NewR2[a,b,i,j]=R2[a,b,i,j]
    elif string=='LF' or string=='UF':
       for a in range(Nvrt):
           for b in range(a):
               for i in range(Nocc):
                   for j in range(Nocc):
                       NewR2[a,b,i,j]=R2[a,b,i,j]
    elif string=='FL' or string=='FU':
       for a in range(Nvrt):
           for b in range(Nvrt):
               for i in range(Nocc):
                   for j in range(i):
                       NewR2[a,b,i,j]=R2[a,b,i,j]
    elif string=='LL' or string=='UU':
       for a in range(Nvrt):
           for b in range(a):
               for i in range(Nocc):
                   for j in range(i):
                       NewR2[a,b,i,j]=R2[a,b,i,j]


    if string=='UF': NewR2=R2-NewR2
    if string=='FU': NewR2=R2-NewR2
    if string=='UU': NewR2=R2-NewR2

    return NewR2

def print_mat4(String,Mat,dim1,dim2,dim3,dim4):
    print('\n - '+String)
    ind=1
    for j in range(dim4):
        for i in range(dim3):
            for b in range(dim2):
                for a in range(dim1):
                    print(' %7d   %20.12f ' % (ind,Mat[a,b,i,j]))

                    ind=ind+1
    return

def print_mat2(String,Mat,dim1,dim2):
    print('\n - '+String)
    ind=1
    for j in range(dim2):
        for i in range(dim1):
            print(' %7d   %20.12f ' % (ind,Mat[i,j]))
            ind=ind+1
    return

def compare_two_vectors(EnvVal,finp1,finp2):
    DataDir=EnvVal['DATA_DIR']
    finp1=DataDir+'/'+finp1
    finp2=DataDir+'/'+finp2
    print('\n - Check differences of two files:')
    print('   + file1 : '+finp1)
    print('   + file2 : '+finp2)

    f1=open(finp1,'r')
    f2=open(finp2,'r')

    tmp1=f1.readline().split()
    tmp2=f2.readline().split()
    dim0=1
    dim1=[]
    dim2=[]
    for val in tmp1:
        dim1.append(int(val))
        dim0=dim0*int(val)
    for val in tmp2:
        dim2.append(int(val))
    #print(dim0)
    print('   + Dimension of file1 = '+str(dim1))
    print('   + Dimension of file1 = '+str(dim2))
    valmax=0.0
    for i in range(dim0):
        val1=float(f1.readline().split()[-1])
        val2=float(f2.readline().split()[-1])
        #val=abs(val1-val2)
        val=abs(abs(val1)-abs(val2))
        if val>valmax:
           valmax=val
           sval1=val1
           sval2=val2
        #print(val)
    print('   + Maximum difference = '+str(valmax))
    print('   + val1 = '+str(sval1))
    print('   + val2 = '+str(sval2))
    f1.close()
    f2.close()


def get_SubMO(string,EnvVal):
    Nocc=EnvVal['NOCC']
    Nvrt=EnvVal['NVRT']
    Nbas=EnvVal['NBAS']
    MO=read_data('Int-MO.txt',EnvVal)
    ind=[]
    for i in range(4):
        if string[i]=='o': 
           ind.append(slice(0,Nocc,1))
        elif string[i]=='v': 
           ind.append(slice(Nocc,Nbas,1))

    SubMO=MO[ind[0],ind[1],ind[2],ind[3]]

    if (EnvVal['HBAR_SUBMO']=='FILE'):
       write_data('Int-G'+string+'.dat',SubMO,EnvVal)
       return
    else:
       return SubMO
