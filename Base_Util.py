import os
import numpy as np

def make_header(string):
    nstr=len(string)+4
    print('\n'+(' '*3) + ('='*nstr))
    print((' '*5) + string + '  ')
    print((' '*3) + ('='*nstr) +'\n')
    return


def write_data(String,mat,dim,lwrite):
    if (lwrite==False):
       return

    if not os.path.exists("Data"):
       os.makedirs("Data")
    fout='Data/Data-'+String+'.txt'
    if os.path.exists(fout):
       os.remove(fout)
    f=open(fout,'w')
    if (len(dim)==1):
       f.write("%7d\n" % (dim[0]))
       for i0 in range(dim[0]):
           f.write("%7d  %15.7f\n" % (i0,mat[i0]))
    elif (len(dim)==2):
       f.write("%7d  %7d\n" % (dim[0],dim[1]))
       for i0 in range(dim[0]):
           for i1 in range(dim[1]):
               f.write("%7d  %7d  %15.7f\n" % (i0,i1,mat[i0,i1]))
    elif (len(dim)==4):
       f.write("%7d  %7d  %7d  %7d\n" % (dim[0],dim[1],dim[2],dim[3]))
       for i0 in range(dim[0]):
           for i1 in range(dim[1]):
               for i2 in range(dim[2]):
                   for i3 in range(dim[3]):
                       f.write("%7d  %7d  %7d  %7d  %15.7f\n" % (i0,i1,i2,i3,mat[i0,i1,i2,i3]))
    f.close()

def read_data(String,EnvVal):
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
       for i0 in range(dim[0]):
           for i1 in range(dim[1]):
               for i2 in range(dim[2]):
                   for i3 in range(dim[3]):
                       mat[i0,i1,i2,i3]=float(f.readline().split()[4])
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
