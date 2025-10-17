#!/usr/bin/env python3
#-*- coding: utf-8 -*-
from scipy.interpolate import griddata
from multiprocessing import Process,Value,Array,Manager,Pool,TimeoutError
import numpy as np
import cmath as cpx
# import panda as pd
# import csv
import sys
import os
import time
import math 
import operator
from functools import reduce
# import cantera as ct
import matplotlib as mtp

#---------------------------------------------------------------------
######################################             fonction          #
#---------------------------------------------------------------------
def Parseur(Varg,n0,strEXIT) :
	Sysa=sys.argv[1:]
	Narg=len(Varg)
	Arg,Arg0=[],''
	for arg in Varg :
		arg='--'+arg ; Arg0+=arg+' '
		if arg in Sysa : Arg.append(True) ; Sysa.pop(Sysa.index(arg))
		else           : Arg.append(False)
	NSysa=len(Sysa)
	exit=('--help' in Sysa or NSysa<n0)
	# if   exit and strEXIT : sys.exit(strEXIT)
	if   exit : sys.exit('\n\n '+strEXIT+'\nOpt : '+Arg0+'\n\n')
	# if   exit and strEXIT : sys.exit('\n\n '+strEXIT+'  ,  ARG : '+Arg0+'\n\n')
	# elif exit             : sys.exit('\n\n ARG : '+Arg0+'\n\n')
	Section('Arg (Sysa,Varg,Arg) :',1,5,'b')
	print(Sysa)
	print(Varg)
	print(Arg)
	print('\n')
	return(Sysa,NSysa,Arg)
#---------------------------------------------------------------------
def Cut(V0) :
	V=[]
	for i in V0 : V.append(i)
	M=[]
	while V :
		M.append( [V.pop(0)] )
		while V and abs(V[0]-M[-1][-1])==1 :
			M[-1].append(V.pop(0))
	return(M)
#---------------------------------------------------------------------
def Convert(Im,plt,dirout,fmout) :
	from matplotlib import cm
	IM=plt.imread(Im)
	fig,ax=plt.subplots()
	ax.matshow(IM,cmap=cm.binary_r)
	body=Im.split('.')[0].split('/')[-1]
	print(body)
	fig.savefig(dirout+body+'.'+fmout)
#---------------------------------------------------------------------
def arrayToStr(data,sep):
    """
    param: data : input array
    return the string containing all species of the array and
    separate it by sep
    """
    out=''
    for el in data:
            out+=sep+str(el)
    return out
#---------------------------------------------------------------------
def arrayToStr2(data,sep):
    """
    param: data : input array
    return the string containing all species of the array and
    separate it by sep
    """
    out=''
    for el in data:
            out+=sep+str(el)
    if len(out)>0: out=out[len(sep):]
    return out
#---------------------------------------------------------------------
def printround(array,n1,n2,sep):
	fm=str(n1) + '.' +str(n2) +'e'
	Str=arrayToStr( [ format(x,fm) for x in array ] ,sep)
	return(Str)
#---------------------------------------------------------------------
def printround2(array,n1,n2,sep,form):
    fm=str(n1) + '.' +str(n2) +form
    Str=arrayToStr2( [ format(x,fm) for x in array ] ,sep)
    return(Str)
#---------------------------------------------------------------------
def Entete(N0,Nox,Set1,Set2,Spe):

    print('\033[0m')
    print('\n\033[1m' + N0*'#'+'\033[47m')
    print('\033[31m'+Set1+(N0-len(Set1))*' ')
    print('\033[31m'+Set2+(N0-len(Set2))*' ')
    print('\033[0m\033[1m'+N0*'-')
    print('\033[0m Species => ', Spe)
    print('\033[0m\033[1m'+N0*'-')
    if Nox: 
        strNox=arrayToStr(Nox,' --- ')
        print('\033[40m\033[32m Nox => '+ strNox+(N0-len(strNox)-8)*' ')
    print('\033[0m\033[1m' + N0*'#'+'\033[0m')

#---------------------------------------------------------------------
def Entete0(N0,Set1,Set2,title):
    
    N1=(N0-(len(title)))//2
    N2=N0-(N1+len(title))

    print('\n\033[28;1;5;3m' + N0*'#')
        
    print('\033[0m\033[28;1;5;3m'+N0*'-'+'\033[0m')
    print('\033[40m\033[32m'+N1*' '+title+N2*' ' )
    print('\033[0m\033[28;1;5;3m'+N0*'-'+'\033[0m')
    
    print('\033[47m\033[31;1;4;5;3m'+Set1)
    print('\033[47m\033[31;1;4;5;3m'+Set2)

    print('\033[0m\033[28;1;5;3m' + N0*'#'+'\033[0m')
#---------------------------------------------------------------------
def Entete1(N0,Set,title):
    
    N1=(N0-(len(title)))//2
    N2=N0-(N1+len(title))

    print('\n\033[1m' + N0*'#')
        
    print('\033[0m\033[1m'+N0*'-'+'\033[0m')
    print('\033[40m\033[32m'+N0*' ')
    print('\033[40m\033[32;1;4m'+N1*' '+title+N2*' '+'\033[0m' )
    print('\033[40m\033[32m'+N0*' ')
    print('\033[0m\033[1m'+N0*'-'+'\033[0m')
    
    for s in Set: print('\033[47m\033[31m'+s+(N0-len(s))*' ' )

    print('\033[0m\033[1m' + N0*'#'+'\033[0m')
#---------------------------------------------------------------------
def Entete2(N0,Set,title):
    
    N1=(N0-(len(title)))//2
    N2=N0-(N1+len(title))

    print('\n\033[33m' + N0*'#')
        
    print('\033[0m\033[33m'+N0*'-'+'\033[0m')
    print('\033[40m\033[32m'+N0*' ')
    print('\033[40m\033[32;1;4m'+N1*' '+title+N2*' '+'\033[0m' )
    print('\033[40m\033[32m'+N0*' ')
    print('\033[0m\033[33m'+N0*'-'+'\033[0m')
    
    for s in Set: print('\033[47m\033[31m'+s+(N0-len(s))*' ' )

    print('\033[0m\033[33m' + N0*'#'+'\033[0m')
#---------------------------------------------------------------------
ENTETE=Entete1
#---------------------------------------------------------------------
def Col(col,txt) :
	if   (col=='r') : return( '\033[1;31m'+txt+'\033[0m' ) 
	elif (col=='y') : return( '\033[1;33m'+txt+'\033[0m' ) 
	elif (col=='g') : return( '\033[1;32m'+txt+'\033[0m' ) 
	elif (col=='b') : return( '\033[1;34m'+txt+'\033[0m' ) 
	else            : return( '\033[0m '+txt+'\033[0m' )
#---------------------------------------------------------------------
def SectionT(txt,t0,N0,N1,col) :
	print(
		N0*'\n'
		+N1*'='+'> '
		+' Dt : {:.3f} | '.format(time.time()-t0)
		+'\033[1;31m'*(col=='r')
		+'\033[1;33m'*(col=='y')
		+'\033[1;32m'*(col=='g')
		+'\033[1;34m'*(col=='b')
		+txt+'\033[0m'
		+N0*'\n'
	)
#---------------------------------------------------------------------
def Section(txt,N0,N1,col) : 
	print(
		N0*'\n'
		+N1*'='+'> '
		+'\033[1;31m'*(col=='r')
		+'\033[1;33m'*(col=='y')
		+'\033[1;32m'*(col=='g')
		+'\033[1;34m'*(col=='b')
		+txt+'\033[0m'
		+N0*'\n'
	)
#---------------------------------------------------------------------
def Section2(t0,txt,N0,N1,col) : 
	print(
		N0*'\n'
		+N1*'='+'> '
		+'t : {0:.3f} s  ,  '.format(time.time()-t0)
		+'\033[31m'*(col=='r')
		+'\033[33m'*(col=='y')
		+'\033[32m'*(col=='g')
		+'\033[34m'*(col=='b')
		+txt+'\033[0m'
		+N0*'\n'
	)
#---------------------------------------------------------------------
# def ColorVJR(E,V) : return( int(E>V[0])*'\033[42m'+int(E>V[1])*'\033[43m'+int(E>V[2])*'\033[41m'  )
def ColorVJR(E,V) : return( int(E>V[0])*'\033[32m'+int(E>V[1])*'\033[33m'+int(E>V[2])*'\033[31m'  )
#---------------------------------------------------------------------
def EXIT() : sys.exit('\n\n EXIT \n\n')
#---------------------------------------------------------------------
def Error(msg) : Section(msg,0,3,'r')
#---------------------------------------------------------------------
def ERROR(msg) : sys.exit(3*'\n'+'\033[31m '+5*'='+'>  ERROR :'+msg+'\033[0m'+3*'\n')
#---------------------------------------------------------------------
######################################             Plot              #
#---------------------------------------------------------------------
def PlotIm(ax,pic,Ext) :
	im=mtp.image.imread(pic)
	(Ni,Nj,Nk)=im.shape
	ax.imshow( im,extent=Ext )
	ax.set_aspect( (Ext[1]/Ext[3])*(Ni/Nj) )
#---------------------------------------------------------------------
def Plot0():
	import matplotlib.pyplot as plt
	import matplotlib
	# from matplotlib import rc

	# matplotlib.rcParams['text.usetex'] = True
	# matplotlib.rcParams['text.latex.unicode'] = True
	# matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
	# matplotlib.rc('text', usetex=True)
	matplotlib.rc('font', family='sans-serif') 
	# matplotlib.rc('font', serif='Helvetica Neue') 
	matplotlib.rc('font', serif='Times New Romand') 
	matplotlib.rc('font', weight='bold') 
	# matplotlib.rc('text', usetex=True)
	matplotlib.rc('text', usetex=False)
	matplotlib.rc('font', size=15 )
	# matplotlib.rcParams.update({'font.size': 25})
	## for Palatino and other serif fonts use:
	# matplotlib.rc('font',**{'family':'serif','serif':['Palatino']})

	return(plt,matplotlib)
#---------------------------------------------------------------------
def Plot(ax,X,Y,N,Fmt,size_font,size_ticks,xlab,ylab,Title,scix,sciy,xsc,ysc,lab) :
	# N=len(X)
	if size_ticks : ax.tick_params(labelsize=size_ticks)
	if Title      : ax.set_title( Title,fontsize=size_font)
	if xlab       : ax.set_xlabel(xlab ,fontsize=size_font)
	if ylab       : ax.set_ylabel(ylab ,fontsize=size_font)
	if scix       : ax.ticklabel_format(axis='x',style='sci',scilimits=scix)
	if sciy       : ax.ticklabel_format(axis='y',style='sci',scilimits=sciy)
	if xsc        : ax.set_xscale(xsc)
	if ysc        : ax.set_yscale(ysc)
	if not lab : lab=N*['']

	for n in range(N): ax.plot(X[n],Y[n],Fmt[n],label=lab[n])

	return(ax)
#---------------------------------------------------------------------
def NewPos(Ax,ax,ay,dx,dy) :
	Pos=Ax.get_position()
	Posp=[list(Pos.p0),list(Pos.p1)]
	Dx=Posp[1][0]-Posp[0][0]
	Dy=Posp[1][1]-Posp[0][1]
	Pos.p0=[Posp[0][0]+dx,Posp[0][1]+dy]
	Pos.p1=[Pos.p0[0]+ax*Dx,Pos.p0[1]+ay*Dy]
	Ax.set_position(Pos)
#---------------------------------------------------------------------
def SaveFig(fig,name) :
	fig.savefig(name)
	Section('Fig saved : '+name,0,2,'g')
#---------------------------------------------------------------------
######################################             Math              #
#---------------------------------------------------------------------
def Nearest(array, value):
    array = np.asarray(array)
    # print('Nearest : {0}   {1}   {2}'.format(array[0],value,array[-1]) )
    idx = (np.abs( array - value )).argmin()
    return( array[idx] , idx )
#---------------------------------------------------------------------
def Suite(x0,x1,dx) : 
	N=1+round(abs(x1-x0)/dx)
	# print(x0,x1)
	# print(N)
	return( list(np.linspace(x0,x1,N))   )
#---------------------------------------------------------------------
def GradDifFin(Np,od,stri,grid,Ngrid,f):
    trace=0 ; Nt=Np[1]+Np[0]+1
    if Ngrid>Nt and Nt>1:
        Df,Test=[],[]
        for p in range(Np[0],Ngrid-Np[1]) :
            MD=np.zeros((Nt,Nt))
            for i in range(Nt):
                for j in range(-Np[0],Np[1]+1) : 
                	# print(i,j,j+Np[0],p+j,p, Ngrid  )
                	MD[i,j+Np[0]]=(grid[p+j]-grid[p])**i/math.factorial(i)
            detD=np.linalg.det(MD) #; print('Det : ',detD)
            # print(MD)
            bD=np.zeros(Nt) ; bD[od]=1
            VD=np.linalg.solve(MD,bD)
            bD2=np.dot(MD,VD)
            test=[ abs(bD[n]-bD2[n]) for n in range(Nt) if n!=od ] ; Test.append(max(test))
            CD=sum( [ VD[n+Np[0]]*(grid[p+n]-grid[p])**od/math.factorial(od) for n in range(-Np[0],Np[1]+1) ] )
            if Test[-1]<stri:
                Df.append( sum( [ f[p+n]*VD[n+Np[0]] for n in range(-Np[0],Np[1]+1) ])/CD )
                if abs(CD-1)>stri : print('CD not 1',CD-1)
            else :
                Df.append(( f[p+1]-f[p])/(grid[p+1]-grid[p]) )
                Test.append(0) ; CD=0; detD=0; trace=1
                print('=> Wrong inversion')
    else :
        print('Not enough point for the chosen order => Ngrid : {0} , Nt : {1}'.format(Ngrid,Nt) )
        Df=[ (f[p+1]-f[p])/(grid[p+1]-grid[p]) for p in range(Ngrid-1) ]
        Test=[0] ; CD=0 ; detD=0 ; trace=2*(Nt>1)

    return(Df,Test,detD,CD,trace)
#---------------------------------------------------------------------
def InterOd1(x0,xn,f0,fn,x): 
	return( ((fn-f0)/(xn-x0))*x + ((fn*x0-f0*xn)/(x0-xn)) )
#---------------------------------------------------------------------
def InterP(F,X1,N,Nx1):

    M=np.zeros((N,N))
    B=np.zeros(N)

    for i in range(N):
        # print('L : ',i)
        # B[i]=sum([F[k]*X1[k]**i for k in range(Nx1)])
        B[i]=sum( np.array(F)*np.array(X1)**i )
        for j in range(N):
            # M[i,j]=sum( [x**(i+j) for x in X1] )
            M[i,j]=sum( np.array(X1)**(i+j) )
        # print( abs(M[i,i]) > sum([ abs(M[i,j]) for j in range(N) if j!=i ]) )

    # print('=> Solving')
    det=np.linalg.det(M)
    V=np.linalg.solve(M,B)
    B2=np.dot(M,V)
    E=[abs(B2[i]-B[i]) for i in range(N)]

    return(V,E,det)
#---------------------------------------------------------------------
def Interp2D(M0,N,X0,Y0,MX1,MY1) :
	(Ni,Nj)=N
	Field=[]
	for i in range(Ni) :
		for j in range(Nj) :
			Field.append(M0[i][j])
	return( griddata( (X0,Y0) , Field ,(MX1,MY1),method='cubic') )
#---------------------------------------------------------------------
def interp2D(M0,N,X0,Y0,MX1,MY1,n,ntot) :
	OUT.append(Interp2D(M0,N,X0,Y0,MX1,MY1))
	print('Matrice n {0}/{1}'.format(n,ntot))
#------------------------------
def Interp2DPara(MM0,N,X0,Y0,MX1,MY1,Np) :
	Ntot=len(MM0)
	with Manager() as man :
		global OUT ; OUT=man.list( [] ) 
		n=1
		with Pool(processes=Np) as pool :
			for M0 in MM0 : pool.apply_async( interp2D , (M0,N,X0,Y0,MX1,MY1,n,Ntot) ) ; n+=1
			while len(OUT)<Ntot : time.sleep(1e-1) #; print(len(OUT))
		OUT2=[ m for m in OUT ]
	return(OUT2)
#---------------------------------------------------------------------
def EvalP(A,X):
    P=[]
    A=list(A)
    A.reverse()
    for x in X:
        p=x*A[0]
        for a in A[1:-1]:
            p+=a
            p*=x
        P.append(p+A[-1])

    return(P)
#---------------------------------------------------------------------
def Roots3(A,er):
    [a,b,c]=A # P(x)=x3+ax2+bx+c

    p=b-(a**2)/3
    q=( (2*(a**2)-9*b)*(a/27) + c )

    D=4*(p**3)+27*(q**2)

#    print('p   q   D =>',p,q,D)

    D2=( D/(4*27) )**(1/2) #; print(D)
    u0=-0.5*q + D2
    v0=-0.5*q - D2
    
#    print('u0 v0 D2',u0,v0,D2)

    if abs(D)<er or np.isnan(D2) :
        x1=3*q/p
        x2=x3=-3*q/(2*p)
        if abs(p)<er and abs(q)>er: print('x1 undetermined',D,p,q,x1)
        if abs(D)>100*er and np.isnan(D2) : print('NAN D2 + D not 0',D,D2)
    else :
        if D>0 :
            u=np.sign(u0)*( abs(u0) )**(1/3)
            v=np.sign(v0)*( abs(v0) )**(1/3)
        elif D<0 :
            u=(u0)**(1/3)
            v=(v0)**(1/3)

#        print('u v',u,v)
        if np.isnan(D2) : print('NAN D2 - D p q =>',D2,D,p,q)
        if any( [np.isnan(u) , np.isnan(v)] ) : print('NAN uv =>',u,v)
        x1= u + v - a/3
        Re = -0.5*u -0.5*v -a/3
        Im = ( 0.5*(3)**0.5 )*( u-v )
        x2=Re+Im*1j
        x3=Re-Im*1j

#    print('x1',x1)
#    print('x2 x3',x2,x3)

    Roots=[ x1,x2,x3 ]

    return(Roots,D)
#---------------------------------------------------------------------
def found0(fe,Param,I,e,Nx):

    Er=10*e
    while (I[1]-I[0])>e or Er>e:

        X=np.linspace(I[0],I[1],Nx)
        # X=np.logspace(I[0],I[1],Nx)
        E=[ fe(x,Param) for x in X ] ; Er=min(E)
        Id=[ int(i) for i in range(Nx) if (abs(E[i]-Er))==0 ] #<e*1e-8 ]

#        print('Found 0 =====> ',I[0],I[1]-I[0],Er,Id)
        
        Multi=False
        LId=len(Id) #; print(Id[-1],Id[0]+len(Id)-1)
        if LId==1: I=[X[ max(0,Id[0]-1) ], X[ min(Id[0]+1,Nx-1) ]]
#        elif LId>Nx//3 : Er=0 ; I=[0,0] ; xf=X[Nx//2] ; Multi=True
        elif Id[-1]==(Id[0]+len(Id)-1) and X[Id[-1]]-X[Id[0]]<e : Er=0 ; I=[0,0] ;  xf=sum([X[i] for i in Id])/LId ; Multi=True
        elif LId>Nx//2 : Er=0 ; I=[0,0] ; xf=sum([X[i] for i in Id])/LId ; Multi=True
        elif (I[1]-I[0])>e*10 or Er>e :
            print('MULTI =====> ',Id,I,'     ',((I[1]-I[0])>e or Er>e))
            Er=0 ; I=[0,0] ; xf=[] ; Multi=True
        
            for i0 in Id : xf.append( found0(fe,Param, [X[ max(0,i0-1) ], X[ min(i0+1,Nx-1) ]],e,Nx ) )
            Nxf=len(xf)
            MD=max(xf)-min(xf)
            if MD<e : print('Equal Out Found0 Dxmax =>',MD) ; xf=xf[Nxf//2]
        else : 
            xf=[ X[i] for i in Id ]
            Er=0 ; I=[0,0] ; Multi=True

    if not Multi: xf=0.5*sum(I)

    try : print('Found zero out ===> ',xf,[fe(x,Param) for x in xf])
    except : 
        er=fe(xf,Param)
        if er>e : print(ColorVJR(er,[e,e*1e2,e*1e4])+'Found zero out ===> \033[0m',xf,er)

    return(xf)
#---------------------------------------------------------------------
def Found0(fe,Set,X,e,Nx,SHOW):

	def Len(a,b) : return( abs(a-b) )

	Er=10*e
	LrI=Len(X[-1],X[0])
	F=[ fe(x,Set) for x in X ]# ; MF=max([abs(fl) for fl in F])

	X0,Nxl,I=[],len(X),[]

	if SHOW>0 :
		print('\n\n=====> Start Found0\n len(F) , min(F), max(F), X[0], X[-1], Isign, LI, LrI \n')

	while LrI>e :
		
		F=[ fe(x,Set) for x in X ]# ; MF=max([abs(fl) for fl in F])
		Isign=[ i for i in range(Nxl-1) if F[i+1]*F[i]<=0 ] ; LI=len(Isign)
		# Isign=[ i for i in range(Nxl-1) if F[i+1]*F[i]<e**2 ] ; LI=len(Isign)
		if SHOW>1 :
			print(X)
			(plt)=Plot0() ; fig,ax=plt.subplots() 
			ax.plot(X,F,'.k')
			ax.plot([ X[i] for i in Isign],[ F[i] for i in Isign],'*r')
			plt.show()
		if SHOW>0 :
			print('{0} {1:2.4e} {2:2.4e} {3:2.4e} {4:2.4e} {5} {6} {7:2.4e}'.format(len(F),min(F),max(F),X[0],X[-1],Isign,LI,LrI))

		if LI>1 :
			LrI,I=0,[0,0]
			print('=> MULTI')
			X0=[ Found0(fe,Set,np.linspace(X[i],X[i+1],Nx),e,Nx,SHOW) for i in Isign ]
		elif LI==1 :
			I=[X[ Isign[0] ], X[ Isign[0]+1 ]]
			LrI=Len(X[-1],X[0])
		else :
			LrI,I=0,[0,0]
			print('=> No zero')
			X0=[]

		X=np.linspace(I[0],I[1],Nx) ; Nxl=Nx

	if LrI>0 and I : X0=[0.5*sum(I)]
	elif not I : X0=0.5*(X0[0]+X0[1])

	return(X0)
#---------------------------------------------------------------------
def Found0_2(fe,Set,X,e,Nx,SHOW):

	def Len(a,b) : return( abs(2*(a-b)/(a+b)) )

	Er=10*e
	LrI=Len(X[-1],X[0])
	F=[ fe(x,Set) for x in X ] ; MF=max([abs(fl) for fl in F])

	X0,Nxl,I=[],len(X),[]

	if SHOW>1 :
		print('\n\n=====> Start Found0\n len(F) , min(F), max(F), X[0], X[-1], Isign, LI, LrI \n')

	while MF>e/100 and LrI>e :
		
		F=[ fe(x,Set) for x in X ] ; MF=max([abs(fl) for fl in F])
		Isign=[ i for i in range(Nxl-1) if F[i+1]*F[i]<0 ] ; LI=len(Isign)
		if SHOW>2 :
			print(X)
			(plt)=Plot0() ; fig,ax=plt.subplots() 
			ax.plot(X,F,'.k')
			ax.plot([ X[i] for i in Isign],[ F[i] for i in Isign],'*r')
			plt.show()
		if SHOW>1 :
			print('{0} {1:2.4e} {2:2.4e} {3:2.4e} {4:2.4e} {5} {6} {7:2.4e}'.format(len(F),min(F),max(F),X[0],X[-1],Isign,LI,LrI))

		if LI>1 :
			LrI,I=0,[0,0]
			print('=> MULTI')
			X0=[ Found0(fe,Set,np.linspace(X[i],X[i+1],Nx),e,Nx,SHOW) for i in Isign ]
		elif LI==1 :
			I=[X[ Isign[0] ], X[ Isign[0]+1 ]]
			LrI=Len(X[-1],X[0])
		else :
			LrI,I=0,[0,0]
			print('=> No zero')
			X0=[]

		X=np.linspace(I[0],I[1],Nx) ; Nxl=Nx

	if LrI>0 and I : X0=[0.5*sum(I)]
	elif not I : X0=0.5*(X0[0]+X0[1])

	return(X0)
#---------------------------------------------------------------------
def Found0_0(fe,Set,X,e,Nx,SHOW):

	Er=10*e
	LrI=abs( X[-1]-X[0] )
	X0,Nxl,I=[],len(X),[]

	# print(LrI,e)
	while LrI>e :
		
		F=[ fe(x,Set) for x in X ]
		Isign=[ i for i in range(Nxl-1) if F[i+1]*F[i]<=0 ] ; LI=len(Isign)

		if LI>1 :
			LrI,I=0,[0,0]
			X0=[ Found0(fe,Set,np.linspace(X[i],X[i+1],Nx),e,Nx,SHOW) for i in Isign ]
		elif LI==1 :
			I=[X[ Isign[0] ], X[ Isign[0]+1 ]]
			LrI=abs( (I[1]-I[0]) )

		else :
			LrI,I=0,[0,0]
			X0=[]

		X=np.linspace(I[0],I[1],Nx) ; Nxl=Nx

	if LrI>0 and I : X0=[0.5*sum(I)]
	elif not I : X0=0.5*(X0[0]+X0[1])

	return(X0)
#---------------------------------------------------------------------
def Norme2(X1,X2) : return( np.linalg.norm(np.array(X1)-np.array(X2)) )
#---------------------------------------------------------------------
def Normale(X1,X2) : (x,y)=X2-X1 ; r=x/y ; n1=1/np.sqrt(1+r**2) ; n2=-r*n1 ;  return( np.array([n1,n2]) )
#---------------------------------------------------------------------
def Angle(X1,X2) : #return( np.dot(X1,X2)/(np.linalg.norm(X1)*np.linalg.norm(X2)) )
	NX1=np.linalg.norm(X1)
	NX2=np.linalg.norm(X2)
	if NX1>0 and NX2>0 : return( np.dot(X1,X2)/(NX1*NX2) )
	else               :
		Section('Zero norm vector',0,5,'r')
		return(60)
#---------------------------------------------------------------------
def SolveLinSys(A,X,X0) :
		(Vp,P)=np.linalg.eig(A) ; P1=np.linalg.inv(P)
		return( [ np.dot(np.dot(np.dot( P,np.diag([np.exp(e*x) for e in Vp])),P1 ),X0) for x in X ] )
#---------------------------------------------------------------------
def Moy2(X,Y,Z,N) :
	seg=np.hypot( X[1:]-X[:-1],Y[1:]-Y[:-1] )
	return( sum(0.5*(Z[1:]+Z[:-1])*seg)/sum(seg) )
#---------------------------------------------------------------------
def Moy(X,Y,N) :
	L=X[-1]-X[0]
	return( sum([ 0.5*(Y[n+1]+Y[n])*(X[n+1]-X[n]) for n in range(N-1) ])/L )
#---------------------------------------------------------------------
def NoiseFilter(Fin,Nt,Nt2,dt,fc1,fc2) :
	four=np.fft.fftshift(np.fft.fft(Fin,n=Nt2))
	freq=np.fft.fftshift(np.fft.fftfreq(Nt2,dt))

	# filt=np.exp(- (freq/fc)**2)
	# filt=(fc1<abs(freq))*(abs(freq)<fc2)+(freq==0)
	filt=np.array(freq==0,dtype=float)
	for n in range(1,10) : filt+=1.0*(abs(abs(freq)-n*fc1)<8e-3) #6e-3
	four2=four*filt
	
	Sign2=np.abs(np.fft.ifft(four2))[:Nt]

	return( Sign2,four,freq,filt )
#---------------------------------------------------------------------
def prod(V) : return( reduce(operator.mul, V, 1) )
#---------------------------------------------------------------------
def SmoothMeanLoop(F,P,Np,Nl) :
        for n in range(Nl) : [F,P]=SmoothMean(F,P,Np)
        return(F,P)
#---------------------------------------------------------------------
def SmoothMean(F,P,Np) :
        F_out=[F[0]] ; D_out=[P[0]] ; Nf=len(F)
        for p in range(1,Np) :
                N,D=0,0
                # for n in range(p+Np+1) :
                for n in range(2*p+1) :
                        N+=F[n]*P[n]
                        D+=P[n]
                if D : F_out.append( N/D )
                else : F_out.append( 0 )
                D_out.append(D)
        for p in range(Np,Nf-Np) :
                N,D=0,0
                for n in range(p-Np,p+Np+1) :
                        N+=F[n]*P[n]
                        D+=P[n]
                if D : F_out.append( N/D )
                else : F_out.append( 0 )
                D_out.append(D)
        for p in range(Nf-Np,Nf-1) :
                N,D=0,0
                for n in range(2*p-Nf+1,Nf) :
                        N+=F[n]*P[n]
                        D+=P[n]
                if D : F_out.append( N/D )
                else : F_out.append( 0 )
                D_out.append(D)
        F_out.append(F[-1])
        D_out.append(P[-1])
        return( np.array(F_out) , np.array(D_out) )
#---------------------------------------------------------------------
# def NoiseFilter(Fin,Nf,Nf2,ft,fc1,fc2) :
# 	MoyF=np.mean(Fin)
# 	F1=np.array(Fin)-MoyF
# 	SupF=max(F1)
# 	F2=F1/SupF
# 	# F2=Fin
# 	#=====> Transform
# 	CFint=np.fft.rfft(F2,n=Nf2)
# 	#=====> Polar
# 	PCFint=np.array([ cpx.polar(c) for c in CFint ])
# 	rCFint,pCFint=PCFint[:,0],PCFint[:,1]*(ft/(2*np.pi))
# 	#=====> Sort
# 	spCFint=np.sort(pCFint)
# 	srCFint=np.array([ rCFint[pCFint==p][0] for p in spCFint ])
# 	#=====> Filter
# 	CFint[   (0<abs(pCFint))*(abs(pCFint)<fc1) ]=0
# 	CFint[ (fc2<abs(pCFint))                   ]=0
# 	#=====> Transform back
# 	Fint=np.fft.irfft( CFint , Nf2 )[:Nf] ; SupFint=max(Fint)
# 	Fint2=Fint*(SupF/SupFint)+MoyF
# 	# Fint2=Fint
# 	return( Fint2,CFint,srCFint,spCFint )
#---------------------------------------------------------------------
def Stat(V) : return( np.min(V),np.mean(V),np.max(V) )
#---------------------------------------------------------------------
######################################             File              #
#---------------------------------------------------------------------
def Stri2Float(stri) :
	try               : return(float(stri))
	except ValueError : return(False)
#---------------------------------------------------------------------
def Stri2Array(stri,sep) : return( [ x for x in [Stri2Float(s) for s in  stri.split(sep)] if x ] )
#---------------------------------------------------------------------
def Sed(L,S,F):
	cmd='sed'
	for i in range(len(L)) :
		cmd+=" -e '{0:.0f}c {1}'".format(L[i],S[i])
	cmd+=' {0} > TEMP ; mv TEMP {0}'.format(F)
	os.system(cmd)
#---------------------------------------------------------------------
def ReadFile(name,sep,Ntitle) :
	L,M,T=True,[],[]
	with open(name) as file :
		for n in range(Ntitle) : T.append(file.readline())
		L=file.readline()
		while L :
			m=[]
			for strl in L.split( sep ):
				try : m.append(float(strl))
				except : pass
			M.append( m )
			L=file.readline()
	file.closed
	return(M,T)
#---------------------------------------------------------------------
# def Titre(f) : return(os.popen('head -n 1 '+f).read())
def Titre(f,sep) : f1=open(f) ; L0=f1.readline() ; f1.closed ; return(L0[:-1].split(sep))
#---------------------------------------------------------------------
def Id(var,V) : return([ i for i in range(len(V)) if var in V[i] ][0])
#---------------------------------------------------------------------
# def ReadCSV(name,sep,Ntitle) :
	# L,M,T=True,[],[]
	# with open(name,newline='') as csvfile :
		# csv_reader=csv.reader(csvfile,delimiter=',')
		# T=csv_reader[:Ntitle]
	# return(M,T)
#---------------------------------------------------------------------
def Line(I,FILE) : return(os.popen("sed -n '{0:.0f}p' {1}".format(I,FILE)).read()[:-1])
#---------------------------------------------------------------------
def MKDIR(dire):
	if not os.path.exists(dire) :
		print('\n=====> Make dir : ',dire)
		# os.system('mkdir '+dire)
		os.mkdir(dire)
#---------------------------------------------------------------------
def WriteMat(M,Title,name,n1,n2,sep,form) :
	with open(name,'w') as file :
		for t in Title : file.write( t+'\n' )
		for m in M     : file.write( printround2(m,n1,n2,sep,form)+'\n' )
	file.closed
