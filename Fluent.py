#!/usr/bin/env python3
#-*- coding: utf-8 -*-

from numpy import *
import h5py as h5
import Utilities as util
from mpl_toolkits.axes_grid1 import make_axes_locatable

(plt,mtp)=util.Plot0()
#===================================================================
def Tri(X,Y) :
# def Tri(X,Y,Geo) :
    # Lin,Rin,Din,Rou=Geo
    tri=mtp.tri.Triangulation( X,Y )
    Xtri=X[ tri.triangles ].mean(axis=1)
    Ytri=Y[ tri.triangles ].mean(axis=1)
    # Mask_ou=(Xtri<Lin)*(Ytri>Rou)
    # Mask_in=(Xtri<Lin)*(Rin<Ytri)*(Ytri<Rin+Din)
    # tri.set_mask(Mask_ou+Mask_in)
    return(tri,Xtri,Ytri)
#===================================================================
# def TriMask(Xc,Yc,Fm) :
    # tri=mtp.tri.Triangulation( Xc,Yc ) ; tri.set_mask(Mask(tri,Xc,Yc))
#===================================================================
def Read_Faces(f_cas,num,names1) :
    cas=h5.File(f_cas)
    Coor=cas['meshes']['1']['nodes']['coords'][num][:]
    names0=str(cas['meshes']['1']['faces']['zoneTopology']['name'][:][0])[2:-1].split(';')
    minId0=cas['meshes']['1']['faces']['zoneTopology']['minId'][:]
    maxId0=cas['meshes']['1']['faces']['zoneTopology']['maxId'][:]
    Fno0=cas['meshes']['1']['faces']['nodes']['1']['nnodes'][:] #; Nf=len(Fno0)
    Fno1=cas['meshes']['1']['faces']['nodes']['1']['nodes'][:]
    Nn=cas['meshes']['1']['nodes']['zoneTopology']['maxId'][0] ; # Number of nodes
    ID_inF=[] # Id of the points in each faces
    FA_inF=[] # Id of the faces connected to the points belonging to each label
    Fcent =[] # center of faces
    Iface =[] # Id face in name
    for f in names1 :
        idf=names0.index(f)
        id0=minId0[idf]
        id1=maxId0[idf]
        Rn=array(range(Nn))
        n0=sum(Fno0[:id0])
        Idf=array(range(id0,id1)) ; Nf=len(Idf)
        Points_faces=[ [] for n in range(Nn) ] # Faces conected to each points
        # Faces_points=[ [] for n in range(Nf) ] # Points in each faces
        fcent=[]
        for i in Idf :
            np=Fno0[i]
            nodes=Fno1[n0:n0+np] # Nodes in the face i
            for n in nodes : Points_faces[n-1].append(i)
            fcent.append(mean(Coor[nodes,:],axis=0))
            n0+=np
        PinF=[ len(g)>0 for g in Points_faces ]
        ID_inF.append(Rn[PinF])
        FA_inF.append([ g for g in Points_faces if g ])
        Fcent.append(array(fcent))
        Iface.append(Idf)
    # Coor=append(Coor,Fcent,axis=0)
    return( Coor,ID_inF,FA_inF,Fcent,Iface,[minId0,maxId0] )
#===================================================================
def Read_Mesh(cas,num) :
    num=list(cas['meshes']['1']['nodes']['coords'].keys())[0]
    Coor=cas['meshes']['1']['nodes']['coords'][num][:] ; Np=len(Coor)
    #=====> Faces
    Fno0=cas['meshes']['1']['faces']['nodes']['1']['nnodes'][:] ; Nf=len(Fno0)
    Fno1=cas['meshes']['1']['faces']['nodes']['1']['nodes'][:]
    Fc0 =cas['meshes']['1']['faces']['c0']['1'][:] ; Nc0=len(Fc0) #print('=> C0 : ',len(Fc0))
    Fc1 =cas['meshes']['1']['faces']['c1']['1'][:] ; Nc1=len(Fc1) #print('=> C1 : ',len(Fc1))
    P_c,P_s , Sel_s = Projection(Fno0,Fno1,Fc0,Fc1,Np,Nc1,Nf) ; print('Side : ',sum(Sel_s))
    return(Coor,P_c)
#===================================================================================
# def Residual(dat,ax) :
#     D_re=dat['results']['residuals']['phase-1']
#     Var_re=['continuity', 'energy', 'fmean', 'fvar', 'k', 'omega', 'premixc', 'premixc-var', 'x-velocity', 'y-velocity']
#     for v in Var_re :
#         #print( '\033[33m=> var : '+v+'\033[0m' )
#         d=array( D_re[v]['data'] )[:,:2]       ; sel0=d[:,1]>0 ; d=d[sel0,0]/d[sel0,1]
#         i=array( D_re[v]['iterations'] )[sel0] 
#         ax.plot( i,d , label=v )
#     ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#===================================================================================
def Projection(Fno0,Fno1,Fc0,Fc1,Np,Nc1,Nf) :
    Ns=Nf-Nc1
    print('=> Projection  Np : {}  Nc1 : {}  Nf : {}  Ns : {}'.format(Np,Nc1,Nf,Ns))
    P_c=[ [] for n in range(Np) ]
    P_s=[ [] for n in range(Np) ]
    n0=0
    for n in range(Nc1) :
        np=Fno0[n]
        for p in Fno1[n0:n0+np] :
            P_c[p-1].append( Fc0[n]-1 )
            P_c[p-1].append( Fc1[n]-1 )
        n0+=np
    for n in range(Nc1,Nf) :
        np=Fno0[n]
        for p in Fno1[n0:n0+np] :
            # P_s[p-1].append( n-Nc1 )
            P_c[p-1].append( Fc0[n]-1 )
        n0+=np
    Sel_s=array([len(p) for p in P_s])>0
    return(P_c,P_s,Sel_s)
#===================================================================
def NodeDataF(D_f,P_f,Nft,Nfi) :
    # print(D_f)
    # print(P_f)
    # if   len(D_f)==Nft     : return( mean(D_f[P_f],axis=0) )
    LDf=len(D_f)
    if   LDf==Nft     : return( array([ mean(D_f[      f                   ]) for f in P_f ]) )
    elif LDf==Nft-Nfi : return( array([ mean(D_f[array(f)-    Nfi          ]) for f in P_f ]) )
    else                   : 
        print('=> Error : Size data not fitting size faces  =>  L data : {}  ,  N faces : {}  ,  N faces 2 : {}'.format(LDf,Nft,Nft-Nfi))
        #return( array([ mean(D_f[array(f)-int(Nft-len(D_f))]) for f in P_f ]) )
#===================================================================
def NodeData(F_c,P_c) : return(array([ mean(F_c[c]) for c in P_c ]))
#===================================================================
def DataF(dat,var) : return(dat['results']['1']['phase-1']['faces'][var]['1'][:])
def DataC(dat,var) : return(dat['results']['1']['phase-1']['cells'][var]['1'][:])
#===================================================================
def Data_Cond(  D_c,P_c           ) : MaxD=max(D_c) ; D_c[D_c==MaxD]=MaxD*(1-1e-13) ; return(        NodeData( D_c,P_c        ))         
def Data_CondF( D_f,P_f,   Nft,Nfi) : MaxD=max(D_f) ; D_f[D_f==MaxD]=MaxD*(1-1e-13) ; return(        NodeDataF(D_f,P_f,Nft,Nfi))
# def Data_CondFC(D_f,P_f,If,Nft,Nfi) : MaxD=max(D_f) ; D_f[D_f==MaxD]=MaxD*(1-1e-13) ; return( append(NodeDataF(D_f,P_f,Nft,Nfi),D_f[If-Nfi*(len(D_f)==(Nft-Nfi))]) )
def Data_CondFC(D_f,P_f,If,Nft,Nfi) : MaxD=max(D_f) ; D_f[D_f==MaxD]=MaxD*(1-1e-13) ; return(        NodeDataF(D_f,P_f,Nft,Nfi))
#===================================================================
def Field_light(fig,ax,tri,F,v,Log,xlim,ylim,cmap,CMask) :
    lab  =v[1]
    lim  =v[3] #; print(lim)
    ticks=v[4]
    ax.set_aspect('equal')
    # ax.set_xticks([])
    # ax.set_yticks([])
    if len(xlim) : ax.set_xlim(xlim[0],xlim[1])
    if len(ylim) : ax.set_ylim(ylim[0],ylim[1])
    if len(lim)==0 : vmin,vmax=min(F),max(F)
    # else : vmin,vmax=lim[0],lim[1] ; ticks=linspace(vmin,vmax,5) ; F[F<vmin]=vmin*(1+1e-12) ; F[F>vmax]=vmax*(1-1e-12)
    else : vmin,vmax=lim[0],lim[1] ; F[F<vmin]=vmin*(1+1e-12) ; F[F>vmax]=vmax*(1-1e-12)
    #=====> Plot
    if Log : f=ax.tricontourf( tri,F,levels=[ 10**n for n in linspace(0,7,101) ] , cmap=cmap ,vmax=vmax , norm = LogNorm() )
    else   : f=ax.tricontourf( tri,F,levels=int(1e2) , cmap=cmap ,vmax=vmax,vmin=vmin )
    # else   : f=ax.tricontourf( tri,F,levels=int(1e2) , cmap=cmap ,vmin=vmin,vmax=vmax )
    #=====> Mask
    if len(CMask) : ax.fill(CMask[0],CMask[1],facecolor='white',edgecolor='black')
    #=====> Colorbar
    divider = make_axes_locatable(ax) ; cax = divider.append_axes("right", size="2%", pad=0.25)
    if len(ticks) : cb=fig.colorbar(f,cax=cax,ticks=ticks)
    else          : cb=fig.colorbar(f,cax=cax)
    cb.set_label(lab,fontsize=20)
#===================================================================
def Field2(tri,F,lab,Log,xlim,ylim,vmax,ticks,cmap,CMask,SAVE,name,fs) :
    fig,ax=plt.subplots(figsize=fs)
    # fig.suptitle(lab,fontsize=20)
    ax.set_title(lab,fontsize=20)
    ax.set_aspect('equal')
    # ax.set_xticks([])
    # ax.set_yticks([])
    if len(xlim) : ax.set_xlim(xlim[0],xlim[1])
    if len(ylim) : ax.set_ylim(ylim[0],ylim[1])
    # if vmax>0 : F[F>vmax]=vmax
    # else      : vmax=max(F) #; print(vmax)
    if not vmax : vmax=max(F)
    # vmax=11
    #=====> Plot
    # ax.triplot(tri,color='red',linewidth=0.1) ; F[:]=0
    if Log : f=ax.tricontourf( tri,F,levels=[ 10**n for n in linspace(0,7,101) ] , cmap=cmap ,vmax=vmax , norm = LogNorm() )
    else   : f=ax.tricontourf( tri,F,levels=int(100) , cmap=cmap ,vmax=vmax )
    #=====> Mask
    if len(CMask) : ax.fill(CMask[0],CMask[1],facecolor='white',edgecolor='black')
    #=====> Colorbar
    divider = make_axes_locatable(ax) ; cax = divider.append_axes("right", size="2%", pad=0.25)
    if len(ticks) : cb=fig.colorbar(f,cax=cax,ticks=ticks,extend='both') #,extendrect=False)
    else          : cb=fig.colorbar(f,cax=cax            ,extend='both') #,extendrect=False)
    # cb.set_label(lab,fontsize=20)
    if SAVE :
        fig.savefig(name) #,dpi=1e3)
        plt.close(fig)
        print('=> {:.3e}  {:.3e}   ,   {} : Saved'.format(min(F),max(F),lab))
    else : return(fig,ax,cb)
#===================================================================
def Field(tri,F,lab,Log,xlim,ylim,vmax,ticks,cmap,CMask,SAVE,name) :
	fig,ax=plt.subplots(figsize=(30,10))
	ax.set_aspect('equal')
	# ax.set_xticks([])
	# ax.set_yticks([])
	if len(xlim) : ax.set_xlim(xlim[0],xlim[1])
	if len(ylim) : ax.set_ylim(ylim[0],ylim[1])
	if not vmax : vmax=max(F)
	# vmax=11
    #=====> Plot
	if Log : f=ax.tricontourf( tri,F,levels=[ 10**n for n in linspace(0,7,101) ] , cmap=cmap ,vmax=vmax , norm = LogNorm() )
	else   : f=ax.tricontourf( tri,F,levels=int(100) , cmap=cmap ,vmax=vmax )
    #=====> Mask
	if len(CMask) : ax.fill(CMask[0],CMask[1],facecolor='white',edgecolor='black')
    #=====> Colorbar
	divider = make_axes_locatable(ax) ; cax = divider.append_axes("right", size="2%", pad=0.25)
	if len(ticks) : cb=fig.colorbar(f,cax=cax,ticks=ticks,extend='both') #,extendrect=False)
	else          : cb=fig.colorbar(f,cax=cax            ,extend='both') #,extendrect=False)
	# cb.set_label(lab,fontsize=20)
	if SAVE :
		# fig.savefig('Plot/Visu-{}-{}.png'.format(Mesh,var))
		fig.savefig(name) #,dpi=1e3)
		plt.close(fig)
		print('=> {:.3e}  {:.3e}   ,   {} : Saved'.format(min(F),max(F),lab))
	else : return(fig,ax,cb)
#===================================================================
def Field_full(tri,F,lab,name,Log,sci,xlim,ylim,vmax,ticks,cmap,SAVE) :
    fig,ax=plt.subplots(figsize=(30,10))
    # fig,ax=plt.subplots(figsize=(10,10))
    ax.set_aspect('equal')
    # ax.set_xticks([])
    # ax.set_yticks([])
    if len(xlim) : ax.set_xlim(xlim[0],xlim[1])
    if len(ylim) : ax.set_ylim(ylim[0],ylim[1])

    if not vmax : vmax=max(F)
    # if Log : f=ax.tricontourf( tri,F,levels=int(1e2) , cmap=cmap ,vmax=vmax , locator=ticker.LogLocator() )
    # if Log : f=ax.tricontourf( tri,F,levels=int(1e2) , cmap=cmap ,vmax=vmax , norm = LogNorm() )
    if Log : f=ax.tricontourf( tri,F,levels=[ 10**n for n in linspace(0,7,101) ] , cmap=cmap ,vmax=vmax , norm = LogNorm() )
    else   : f=ax.tricontourf( tri,F,levels=int(1e2) , cmap=cmap ,vmax=vmax )
    # ax.triplot(tri,color='white')

    divider = make_axes_locatable(ax) ; cax = divider.append_axes("right", size="2%", pad=0.25)
    if len(ticks) : cb=fig.colorbar(f,cax=cax,ticks=ticks)
    else          : cb=fig.colorbar(f,cax=cax)
    cb.set_label(lab,fontsize=20)
    if sci :
        cb.formatter.set_scientific(True)
        cb.formatter.set_powerlimits((-2,2))
        cb.formatter.set_useMathText(True)
        cb.update_ticks()
    if SAVE :
        # fig.savefig('Plot/Field-{}-{}.png'.format(Mesh,var))
        fig.savefig(name)
        plt.close(fig)
        print('=> {:.3e}  {:.3e}   ,   {} : Saved'.format(min(F),max(F),lab))
    # else : return(fig,ax,cb)
    else : return(fig,ax,cb,f)
#===================================================================
def RelSpeed(Curv,U,V) :
    Np=len(Curv)
    Prod=[]
    Norm=[]
    for n in range(1,Np-1) :
        nm=-util.Normale(Curv[n-1,:],Curv[n+1,:]) ; Norm.append(nm) #; print(nm)
        nv=hypot(U[n],V[n])
        Prod.append( abs(nm[0]*U[n]+nm[1]*V[n]) )
        if Prod[-1]/nv>1 : print(util.Col('r','=====> Scalar prod > 1'))
    return( Prod,array(Norm) )
#===================================================================
def Residual(ax,fdat,I0,I1) :
    F_re=h5.File(fdat,'r') ; D_re=F_re['results']['residuals']['phase-1'] ; Var_re=list( D_re.keys() )
    for v in Var_re :
        print( '\033[33m=> var : '+v+'\033[0m' )
        d=array( D_re[v]['data'] )[:,:2]       ; sel0=d[:,1]>0 ; d=d[sel0,0]/d[sel0,1]
        i=array( D_re[v]['iterations'] )[sel0]  
        if I0 or I1 : sel = (I0<=i)*(i<=I1)
        else        : sel=array(len(i)*[True])
        d[d>100]=100
        ax.plot( i[sel],d[sel] , label=v )
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#===================================================================
def Residual_full(ax,fdat,I0,nc,Nm,It) :
    F_re=h5.File(fdat,'r') ; D_re=F_re['results']['residuals']['phase-1'] ; Var_re=list( D_re.keys() )
    for v in Var_re :
        print( '\033[33m=> var : '+v+'\033[0m' )
        d=array( D_re[v]['data'] )[:,:2]       ; print(d.shape,d[0],d[1],d[-2],d[-1]) ; sel0=d[:,1]>0 ; d=d[sel0,0]/d[sel0,1]
        i=array( D_re[v]['iterations'] )[sel0] ; print(i.shape,i[0],i[1],i[-2],i[-1]) ; 
        sel=i>=It[0] #; n1=list(i[sel]).index(I0)
        #sel=d<1
        #n1=list(i).index(I0)
        n1=nc
        # ax.plot( Med2(i[sel],n1,Nm) , Med2(d[sel],n1,Nm) , label=v )
        ax.plot( i[sel],d[sel] , label=v )
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#===================================================================
def Report_read(i_sim,o_sim) :
    D_in=loadtxt(i_sim,skiprows=3,delimiter=' ') #; M_in=D_in[:,1]+D_in[:,2] ; It=D_in[:,0] ; print('=> report : ',M_in.shape) ; print(It[0],It[-1])
    D_ou=loadtxt(o_sim,skiprows=3,delimiter=' ') #; M_ou=D_ou[:,1] ; M_re=M_in+M_ou ; print('=> Final mflow : ', mean(M_re[-100]) )
    return(D_in,D_ou)
#===================================================================
def tick_k(ax) :
    # print(ax.get_xticks())
    xticks=ax.get_xticks()
    # xticks=arange(13.7e3,14.5001e3,200)
    xlabels = ['{:.1f}'.format(x) + 'K' for x in xticks/1000]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)
#===================================================================
def Report(It,Min,Mou,Mth,I0,I1,Nav,coef) :
    Sel=(I0<=It)*(It<=I1) #; print(I0,I1)
    Nin,Nou=len(Min),len(Mou) ; Nmin=min(Nin,Nou)
    Min=Min[:Nmin]*coef
    Mou=Mou[:Nmin]*coef
    It =It [:Nmin]
    Bal=Min+Mou
    Mth*=1e3
    M0=Min[Sel][-Nav:]
    B0=Bal[Sel][-Nav:]
    Ml=mean(M0)
    Bl=mean(B0)
    Vm=sqrt(mean((M0-Ml)**2))/Ml
    Vb=sqrt(mean((B0-Bl)**2))/Ml
    Er=abs(Ml-Mth)/Mth #; print('Ml : {:.2f}  ,  Mth : {:.2f}  ,  Er : {:.1e}'.format(Ml,Mth,Er*100))
    if abs(Bl)>1e3*Ml : tbal=''
    else              : tbal='Balance : {:.1e} [g/s]  ,  RMS : {:.1e} [%]'.format(Bl,Vb*1e2)
    print( 'Min_av : {:.2f}   ,  Min_av_min : {:.2f}   ,  Min_av_max : {:.2f}'.format(Ml,min(M0),max(M0)) )
    print( 'Bal_av : {:.1e}   ,  Bal_av_min : {:.1e}   ,  Bal_av_max : {:.1e}'.format(Bl,min(B0),max(B0)) )
    print( 'RMS m,b : {:.3f}  {:.3f}'.format(Vm*1e2,Vb*1e2) )
    fig_b,ax_b=plt.subplots(figsize=( 15,5),ncols=2)
    #ax_b[0].set_title('Mf : {:.2f}  ,  Mf_th : {:.2f} [g/s]'.format(Ml,Mth*1e3) )
    fig_b.suptitle('Teoretical Mf : {:.2f} [g/s]  ,  Error : {:.2f} [%]'.format( Mth,Er*1e2 ),y=0.995)
    ax_b[0].set_title('Mf : {:.2f} [g/s]  ,  RMS : {:.2e} [%]'.format(Ml,Vm*1e2) )
    ax_b[1].set_title(tbal)
    # ax_b[1].set_title('Balance : {:.2f} [g/s]  ,  RMS : {:.2f} [%]'.format(Bl,Vb*1e2) )
    ax_b[0].ticklabel_format( axis='y' , scilimits=(-2,2) ) ; 
    ax_b[1].ticklabel_format( axis='y' , scilimits=(-2,2) ) ; 
    # ax_b[1].set_ylim((-20,20))
    # ax_b[0].set_yscale('log')
    # ax_b[1].set_yscale('log')
    ax_b[0].plot( It[Sel],Min[Sel],'k' )
    ax_b[1].plot( It[Sel],Bal[Sel],'k' )
    util.NewPos(ax_b[0],1,1,0,-0.04)
    util.NewPos(ax_b[1],1,1,0,-0.04)
    tick_k(ax_b[0])
    tick_k(ax_b[1])
    return(fig_b,ax_b)
#===================================================================
def Temp_FC(dat,P_f,I_f,Nft,Nfi,Thermo) :
    [R,MH2,MO2,MN2,gam]=Thermo ; Patm=101325
    Lvar=list(dat['results']['1']['phase-1']['faces'].keys()) #; print(Lvar)
    if   'SV_FMEAN' in Lvar :
        Z=Data_CondFC( DataF(dat,'SV_FMEAN')   , P_f,I_f,Nft,Nfi ) #; print('Z',Z)
        P=Data_CondFC( DataF(dat,'SV_P')       , P_f,I_f,Nft,Nfi ) #; print('P',P)
        D=Data_CondFC( DataF(dat,'SV_DENSITY') , P_f,I_f,Nft,Nfi ) #; print('D',D)
        YH2=copy(Z)
        YO2=0.233*(1-Z)
        YN2=0.767*(1-Z)
        M=1/( YH2/MH2+YO2/MO2+YN2/MN2 ) #; print('=> Molar Mass : {:.3f}  {:.3f}'.format(min(M),max(M)))
        return(M*(P+Patm)/(D*R))
    elif 'SV_T' in Lvar :
        return(Data_CondFC( DataF(dat,'SV_T'),P_f,I_f,Nft,Nfi )) #; print('=====> No Combu')
#===================================================================
def Vel_FC(dat,P_f,I_f,Nft,Nfi) :
    Var=list(dat['results']['1']['phase-1']['faces'].keys())
    U=Data_CondFC( DataF(dat,'SV_U'),P_f,I_f,Nft,Nfi )
    V=Data_CondFC( DataF(dat,'SV_V'),P_f,I_f,Nft,Nfi )
    if 'SV_W' in Var :
        W=Data_CondFC( DataF(dat,'SV_W'),P_f,I_f,Nft,Nfi )
        return(    sqrt(U**2+V**2+W**2) )
    else : return(hypot(U   ,V        ) )
#===================================================================
def Vel_Face(dat,P_f,Nft,Nfi) :
    Var=list(dat['results']['1']['phase-1']['faces'].keys())
    U=Data_CondF( DataF(dat,'SV_U'),P_f,Nft,Nfi )
    V=Data_CondF( DataF(dat,'SV_V'),P_f,Nft,Nfi )
    if 'SV_W' in Var :
        W=Data_CondF( DataF(dat,'SV_W'),P_f,Nft,Nfi )
        return(    sqrt(U**2+V**2+W**2) )
    else : return(hypot(U   ,V        ) )
#===================================================================
def Vel_C(dat) :
    Var=list(dat['results']['1']['phase-1']['faces'].keys())
    U=DataC(dat,'SV_U')
    V=DataC(dat,'SV_V')
    if 'SV_W' in Var :
        W=DataC(dat,'SV_W')
        return(    sqrt(U**2+V**2+W**2) )
    else : return(hypot(U   ,V        ) )
#===================================================================
def Scalar_Diss(dat,P_f,I_f,Nft,Nfi) :
    #===================> Param
    Cmu=0.09
    Csd=2
    #===================> Data
    K_c=Data_CondFC( DataF(dat,'SV_K')    , P_f,I_f,Nft,Nfi) #dat['results']['1']['phase-1']['cells']['SV_K'   ]['1'][:]
    O_c=Data_CondFC( DataF(dat,'SV_O')    , P_f,I_f,Nft,Nfi) #dat['results']['1']['phase-1']['cells']['SV_O'   ]['1'][:]
    V_c=Data_CondFC( DataF(dat,'SV_FVAR') , P_f,I_f,Nft,Nfi) #dat['results']['1']['phase-1']['cells']['SV_FVAR']['1'][:]
    Z_c=Data_CondFC( DataF(dat,'SV_FMEAN'), P_f,I_f,Nft,Nfi ) #; print(Z)
    #===================> Computation
    E_c=Cmu*K_c*O_c
    X_c=Csd*E_c*V_c/K_c
    return(X_c,E_c,K_c,O_c,V_c,Z_c)
#===================================================================
def Mach_Face(dat,Thermo,PLOT,tri,P_f,I_f,Nft,Nfi,xlim,ylim,cmap,CMask,Mesh,dirP) :
    [R,MH2,MO2,MN2,gam]=Thermo ; Patm=101325
    # print('FIELD IN : ',PLOT)
    # T=Temp_FC(              dat              ,P_f,I_f,Nft,Nfi,Thermo) #; print(T)
    if 'SV_Y' in list(dat['results']['1']['phase-1']['faces'].keys()) :
        P=Data_CondFC( DataF(dat,'SV_P'),P_f,I_f,Nft,Nfi )
        T=Data_CondFC( DataF(dat,'SV_T'),P_f,I_f,Nft,Nfi )
        # YC=DataC(dat,'SV_Y')
        YF=DataF(dat,'SV_Y')
        YH2=Data_CondFC( YF[:,0]        ,P_f,I_f,Nft,Nfi )
        YO2=Data_CondFC( YF[:,1]        ,P_f,I_f,Nft,Nfi )
        YN2=Data_CondFC( YF[:,2]        ,P_f,I_f,Nft,Nfi )
        M=1/( YH2/MH2+YO2/MO2+YN2/MN2 )
        Cel=sqrt( gam*R*T/M )
    elif 'SV_FMEAN' in list(dat['results']['1']['phase-1']['faces'].keys()) :
        Z=Data_CondFC( DataF(dat,'SV_FMEAN')   , P_f,I_f,Nft,Nfi ) #; print(Z)
        P=Data_CondFC( DataF(dat,'SV_P')       , P_f,I_f,Nft,Nfi ) #; print(Z)
        D=Data_CondFC( DataF(dat,'SV_DENSITY') , P_f,I_f,Nft,Nfi ) #; print(Z)
        YH2=copy(Z)
        YO2=0.233*(1-Z)
        YN2=0.767*(1-Z)
        M=1/( YH2/MH2+YO2/MO2+YN2/MN2 ) #; print('=> Molar Mass : {:.3f}  {:.3f}'.format(min(M),max(M)))
        T=M*(P+Patm)/(D*R) #; print(T)
        Cel=sqrt( gam*(P+Patm)/D )
    else : print('=====> Error : No SV_Y or SV_FMEAN in dat')
    Vel=Vel_FC(dat,P_f,I_f,Nft,Nfi)
    if PLOT :
        # print('Plot Mach Fields => dirP : ',dirP)
        Field( tri, P   ,'Pressure [Pa]'      ,False,xlim,ylim,0,linspace(0,1e7,11),cmap,CMask, True,dirP+'/Visu-{}-{}.png'.format(Mesh,'Pres') )
        Field( tri, T   ,'Temperature [K]'    ,False,xlim,ylim,0,[]                ,cmap,CMask, True,dirP+'/Visu-{}-{}.png'.format(Mesh,'Temp') )
        Field( tri, YH2 ,'Y H2 [-]'           ,False,xlim,ylim,0,linspace(0,1  ,11),cmap,CMask, True,dirP+'/Visu-{}-{}.png'.format(Mesh,'YH2' ) )
        Field( tri, YO2 ,'Y O2 [-]'           ,False,xlim,ylim,0,linspace(0,1  ,51),cmap,CMask, True,dirP+'/Visu-{}-{}.png'.format(Mesh,'YO2' ) )
        Field( tri, YN2 ,'Y N2 [-]'           ,False,xlim,ylim,0,linspace(0,1  ,51),cmap,CMask, True,dirP+'/Visu-{}-{}.png'.format(Mesh,'YN2' ) )
        Field( tri, M   ,'Molar mass [kg/mol]',False,xlim,ylim,0,[]                ,cmap,CMask, True,dirP+'/Visu-{}-{}.png'.format(Mesh,'Mm'  ) )
        Field( tri, Vel ,'Velociy [m/s]'      ,False,xlim,ylim,0,[]                ,cmap,CMask, True,dirP+'/Visu-{}-{}.png'.format(Mesh,'Vel' ) )
        Field( tri, Cel ,'Sound speed [m/s]'  ,False,xlim,ylim,0,[]                ,cmap,CMask, True,dirP+'/Visu-{}-{}.png'.format(Mesh,'Cel' ) )
    # return(Vel/Cel)
    return(P,T,YH2,YO2,YN2,M,Vel,Cel,Vel/Cel)
#===================================================================
def Mach_Number(dat,Thermo,PLOT,tri,P_c,xlim,ylim,cmap,CMask,Mesh,dirP) :
	[R,MH2,MO2,MN2,gam]=Thermo
	T=DataC(dat,'SV_T')
	Y=DataC(dat,'SV_Y') ; YH2,YO2,YN2=Y[:,0],Y[:,1],Y[:,2]
	M=1/( YH2/MH2+YO2/MO2+YN2/MN2 )
	Vel=Vel_C(dat)
	Cel=sqrt( gam*R*T/M )
	if PLOT :
		Field( tri, Data_Cond(T  ,P_c),'Temperature [K]'    ,False,xlim,ylim,0,[]              ,cmap,CMask, True,dirP+'/Visu-{}-{}.png'.format(Mesh,'Temp') )
		Field( tri, Data_Cond(YH2,P_c),'Y H2 [-]'           ,False,xlim,ylim,0,linspace(0,1,11),cmap,CMask, True,dirP+'/Visu-{}-{}.png'.format(Mesh,'YH2' ) )
		Field( tri, Data_Cond(YO2,P_c),'Y O2 [-]'           ,False,xlim,ylim,0,linspace(0,1,11),cmap,CMask, True,dirP+'/Visu-{}-{}.png'.format(Mesh,'YO2' ) )
		Field( tri, Data_Cond(YN2,P_c),'Y N2 [-]'           ,False,xlim,ylim,0,linspace(0,1,11),cmap,CMask, True,dirP+'/Visu-{}-{}.png'.format(Mesh,'YN2' ) )
		Field( tri, Data_Cond(M  ,P_c),'Molar mass [kg/mol]',False,xlim,ylim,0,[]              ,cmap,CMask, True,dirP+'/Visu-{}-{}.png'.format(Mesh,'Mm'  ) )
		Field( tri, Data_Cond(Vel,P_c),'Velociy [m/s]'      ,False,xlim,ylim,0,[]              ,cmap,CMask, True,dirP+'/Visu-{}-{}.png'.format(Mesh,'Vel' ) )
		Field( tri, Data_Cond(Cel,P_c),'Sound speed [m/s]'  ,False,xlim,ylim,0,[]              ,cmap,CMask, True,dirP+'/Visu-{}-{}.png'.format(Mesh,'Cel' ) )
	return(Vel/Cel)
#===================================================================
def Boundary(cas,U_p) :
	Type=cas['meshes']['1']['faces']['zoneTopology']['zoneType'][:] 
	Walls=(Type==3) ; Nw=sum(Walls)
	Fno0=cas['meshes']['1']['faces']['nodes']['1']['nnodes']
	Fno1=cas['meshes']['1']['faces']['nodes']['1']['nodes']
	WI0 =cas['meshes']['1']['faces']['zoneTopology']['minId'][Walls]
	WI1 =cas['meshes']['1']['faces']['zoneTopology']['maxId'][Walls]
	for w in range(Nw) :
		i0,i1=WI0[w],WI1[w]
		n0=sum(Fno0[:i0])
		for i in range(i0,i1) :
			n1=Fno0[i0]
			U_p[ Fno1[n0:n0+n1] ]=0
			n0+=n1
#===================================================================
# def Profile(tri,Dp,Xint,Yint,dirP,Mesh,var) :
def Profile(tri,Dp,Long,Tran,dirP,Mesh,var) :
    [Xint,Yint]     =Long
    [Xrad,Yrad,Nrad]=Tran
    #======> Interpolation
    f_M=mtp.tri.LinearTriInterpolator( tri,Dp ) ; Prof_M=f_M(Xint,Yint)
    # f_M=mtp.tri.CubicTriInterpolator( tri,Dp,kind='geom' ) ; Prof_M=f_M(Xint,Yint)
    #======> Output Longi
    f=open(dirP+'/Profile-{}-{}.dat'.format(Mesh,var),'w')
    for n in range(len(Xint)) : f.write('{:.12e},{:.12e},{:.12e}\n'.format(Xint[n],Yint[n],Prof_M[n]))
    f.closed
    #======> Output Transverse
    for i in range(len(Xrad)) :
        xint=Nrad[i]*[Xrad[i]]
        yint=linspace(0,Yrad[i],Nrad[i])
        Prof_l=f_M(xint,yint)
        f=open(dirP+'/Profile-X{}-{}-{}.dat'.format(Xrad[i],Mesh,var),'w')
        for n in range(Nrad[i]) : f.write('{:.12e},{:.12e},{:.12e}\n'.format(xint[n],yint[n],Prof_l[n]))
        f.closed
#===================================================================
def ChocMass(P0,T0,Sc,Thermo) :
	[R,MH2,MO2,MN2,gam]=Thermo
	return( P0*Sc*sqrt(MH2*gam/(R*T0))*(0.5*(gam+1))**(-0.5*(gam+1)/(gam-1)) )
#===================================================================
def LastMass(frep,Nav) :
	D_in=loadtxt(frep,skiprows=3,delimiter=' ') ; Min=D_in[:,1]
	return(mean(Min[:-Nav]))
#===================================================================