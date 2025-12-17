#!/usr/bin/env python3
#-*- coding: utf-8 -*-

from numpy import *
import h5py as h5
import Utilities as util
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as colors

(plt,mtp)=util.Plot0()
Mol_m={
    'C':12,
    'O':16,
    'H':1,
    'CH4':16,
    'CO2':44,
    'CO':28,
    'H2O':18,
    'O2':32,
    'N2':28,
    'H2':2
    }

Spe_Laera_l0=['CH4','H2','O2','CO2']
Spe_Laera_l1=['O2','H2O','CH4','CO','H2','H','O','OH','HO2','H2O2','CH3','CH2O','CH3O','CH3OH','C2H2','C2H4','C2H6','CH2CO','CH','CH2','CH2(S)','HCO','CH2OH','C2H3','C2H5','HCCO','CH2CHO','CO2']
Spe_Laera=['O2','H2O','CH4','CO','CO2','H2','H','O','OH','HO2','H2O2','CH3','CH2O','CH3O','CH3OH','C2H2','C2H4','C2H6','CH2CO','CH','CH2','CH2(S)','HCO','CH2OH','C2H3','C2H5','HCCO','CH2CHO','N2']
# Spe_Walter=['O2','H2O','CH4','CO','H2','C2H6','C3H8','CO2']
Spe_Walter=['CH4','C2H6','C3H8','H2','O2','CO2','CO','H2O','N2']
Spe_H2Air=['H2','O2','H2O','N2']
Spe_UCSD=['H2','H','O2','OH','O','H2O','HO2','H2O2','N2']

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
def PlotVar(ax,c,var,X,M,T,TXT,BD,PARAM) :
    (Xtxt,txt)=TXT
    if var=='mix' :
        if 'ch4' in T : 
            [BC_f,BC_o]=PARAM
            Ic1=T.index('ch4')
            Ic2=T.index('co2')
            Ic3=T.index('co')
            Zc=Yc( {'CH4':M[:,Ic1],'CO2':M[:,Ic2],'CO':M[:,Ic3]} , Mol_m )
            Yc_f=Yc(BC_f,Mol_m)
            Yc_o=Yc(BC_o,Mol_m)
            Var=(Zc-Yc_o)/(Yc_f-Yc_o)
        elif 'h2' in T :
            [BC_f,BC_o]=PARAM
            Ih1=T.index('h2')
            Ih2=T.index('h2o')
            Zh=Yh( {'CH4':0*M[:,Ih1],'H2':M[:,Ih1],'H2O':M[:,Ih2]} , Mol_m )
            Yh_f=Yh(BC_f,Mol_m)
            Yh_o=Yh(BC_o,Mol_m)
            Var=(Zh-Yh_o)/(Yh_f-Yh_o)
    elif var=='co' :
        I=T.index('co')
        Var=M[:,I]*1e2
    elif var=='no' :
        I=T.index('mf-pollut-pollutant-0')
        Var=M[:,I]*1e6
    else :
        I=T.index(var)
        Var=M[:,I]
    ax.plot( X,Var,c )
    ax.text( Xtxt[0],Xtxt[1],txt )
    ax.set_xlim((BD[0],BD[1]))
    ax.set_ylim((BD[2],BD[3]))
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
def Space(T,M) :
    Id,Ms=[],[]
    if 'x-coordinate' in T : Ix=FindData('x-coordinate',T) ; Mx=M[:,Ix] ; Id.append(Ix) ; Ms.append(Mx) #; print('=> Mx : {:.1f} , {:.1f}'.format(min(Mx),max(Mx)))
    if 'y-coordinate' in T : Iy=FindData('y-coordinate',T) ; My=M[:,Iy] ; Id.append(Iy) ; Ms.append(My) #; print('=> My : {:.1f} , {:.1f}'.format(min(My),max(My)))
    if 'z-coordinate' in T : Iz=FindData('z-coordinate',T) ; Mz=M[:,Iz] ; Id.append(Iz) ; Ms.append(Mz) #; print('=> Mz : {:.1f} , {:.1f}'.format(min(Mz),max(Mz)))
    return(Id,Ms)
#===================================================================
def SpaceD(D) :
    Ms=[]
    if 'x-coordinate' in D.keys() : Mx=D['x-coordinate'] ; Ms.append(Mx)
    if 'y-coordinate' in D.keys() : My=D['y-coordinate'] ; Ms.append(My)
    if 'z-coordinate' in D.keys() : Mz=D['z-coordinate'] ; Ms.append(Mz)
    return(Ms)
#===================================================================
def ReadSurf(f) : 
    op=open(f) ; L0=op.readline() ; op.closed
    return( [ s.strip() for s in L0[:-1].split(',')] , loadtxt(f,skiprows=1,delimiter=',') )
#===================================================================
def ReadSurfD(f) :
    op=open(f) ; L0=op.readline() ; op.closed
    T=[ s.strip() for s in L0[:-1].split(',')]
    M=loadtxt(f,skiprows=1,delimiter=',')
    return({ s:M[:,n] for n,s in enumerate(T) })
#===================================================================
# def FindData(v,T) : return([ i for i in range(len(T)) if v in T[i] ][0])
def FindData(v,T) : return( T.index(v) )
#===================================================================
def Close(tri,M,tol) : return(abs(M[tri.triangles]-min(M))<tol)
#===================================================================
def CleanTri(tri,tol) :
    xy=dstack((tri.x[tri.triangles],tri.y[tri.triangles]))
    S2=cross(xy[:,1,:]-xy[:,0,:],xy[:,2,:]-xy[:,0,:])
    return(S2<tol)
#===================================================================
def Visu(surf,var,lab,xlim,ylim,ticks,TICKS,BD,fs,cmap0,name,OPT) :
    cmap=mtp.colormaps[cmap0]
    tol=1e-5
    if var in ['tt'] : Log=True
    else             : Log=False
    # cmesh=1e3
    # xlim=[ x*cmesh for x in xlim ]
    # ylim=[ y*cmesh for y in ylim ]
    vmax,vmin=0,0
    if len(BD)>0 : [vmin,vmax]=BD
    (T,M)=ReadSurf(surf)
    if var=='tourb' :
        Mk=M[:,FindData('turb-kinetic-energy',T)]
        Me=M[:,FindData('turb-diss-rate',T)]
        Mv=100*0.09*Mk**1.5/Me
    elif var=='mixC' :
        iopt=OPT.index('MIXC')
        [Mol_m,Y_f,Y_o]=OPT[iopt+1]
        Ic1=T.index('ch4')
        Ic2=T.index('co2')
        Ic3=T.index('co')
        Yc_f=Yc(Y_f,Mol_m)
        Yc_o=Yc(Y_o,Mol_m)
        Yc_g=Yc({'CH4':M[:,Ic1],'CO2':M[:,Ic2],'CO':M[:,Ic3]},Mol_m)
        Mv=(Yc_g-Yc_o)/(Yc_f-Yc_o)
    elif var=='mixH' or ('Zst' in OPT and var in ['temperature','tt'] and not 'fmean' in T ) :
        iopt=OPT.index('MIXH')
        [Mol_m,Y_f,Y_o]=OPT[iopt+1]
        Y_h2 =M[:,T.index('h2')]
        Y_h2o=M[:,T.index('h2o')]
        if 'ch4' in T : Y_ch4=M[:,T.index('ch4')]
        else          : Y_ch4=0*Y_h2
        # print('=> MIXH : ',Mol_m['H2'],Y_f['H2'],Y_o['H2'])
        Yh_f=Yh(Y_f,Mol_m)
        Yh_o=Yh(Y_o,Mol_m)
        Yh_g=Yh({'H2':Y_h2,'H2O':Y_h2o,'CH4':Y_ch4},Mol_m)
        Yh0=(Yh_g-Yh_o)/(Yh_f-Yh_o)
        if var=='mixH' : Mv=Yh0
        elif var=='tt' :
            Mk=M[:,FindData('turb-kinetic-energy',T)]
            Me=M[:,FindData('turb-diss-rate',T)]
            # Mv=clip(Mk/Me,ticks[0],ticks[-1])
            Mv=clip(Me/Mk,ticks[0],ticks[-1])
        else : Ivr=T.index(var) ; Mv=M[:,Ivr]
    elif var=='co' :
        Ivr=T.index('co') 
        if 'CO' in OPT : Mv=M[:,Ivr]*OPT[OPT.index('CO')+1]
        else           : Mv=M[:,Ivr]
    # else : Ivr=FindData(var ,T) ; Mv=M[:,Ivr]
    elif var=='no' :
        Ivr=T.index('mf-pollut-pollutant-0') ; Mv=M[:,Ivr]*OPT[OPT.index('NO')+1]
    elif var=='tt' :
        Mk=M[:,FindData('turb-kinetic-energy',T)]
        Me=M[:,FindData('turb-diss-rate',T)]
        # Mv=clip(Mk/Me,ticks[0],ticks[-1])
        Mv=clip(Me/Mk,ticks[0],ticks[-1])
    else : Ivr=T.index(var) ; Mv=M[:,Ivr]
    # Ivr=T.index(var)
    Ibd=FindData('boundary-cell-dist',T)
    Ivl=FindData('velocity-magnitude',T)
    if 'x-coordinate' in T : Ix=FindData('x-coordinate',T) ; Mx=M[:,Ix] ; Mx0,Mx1=min(Mx),max(Mx) ; print( '=> Mx : {:.1f} , {:.1f}'.format(Mx0,Mx1) )
    if 'y-coordinate' in T : Iy=FindData('y-coordinate',T) ; My=M[:,Iy] ; My0,My1=min(My),max(My) ; print( '=> My : {:.1f} , {:.1f}'.format(My0,My1) )
    if 'z-coordinate' in T : Iz=FindData('z-coordinate',T) ; Mz=M[:,Iz] ; Mz0,Mz1=min(Mz),max(Mz) ; print( '=> Mz : {:.1f} , {:.1f}'.format(Mz0,Mz1) )
    Selzx=[]
    XY= ('xy' in surf) or (not 'z-coordinate' in T) 
    ZX='zx' in surf                                
    YZ= 'yz' in surf or 'tga' in surf
    IN='in1' in surf or 'in2' in surf
    OU='ou'  in surf              
    if   XY : tri=mtp.tri.Triangulation(Mx,My)
    elif ZX : tri=mtp.tri.Triangulation(Mx,Mz) ; Selzx=any(Close(tri,Mx,tol)*Close(tri,Mz,tol),axis=1)
    elif YZ : tri=mtp.tri.Triangulation(My,Mz)
    elif IN : tri=mtp.tri.Triangulation(My,Mz)
    elif OU : tri=mtp.tri.Triangulation(My,Mx)
    if max(M[:,Ibd])>2 : Mask0=sum(M[tri.triangles,Ibd]<1.01,axis=1)==3
    else               : Mask0=sum(M[tri.triangles,Ivl]==0  ,axis=1)==3
    if len(Selzx)>0 : Mask0[Selzx]=False
    if vmax!=0 : MaskV=all(Mv[tri.triangles]>vmax,axis=1) ; Mask0[MaskV]=True ; Mv[Mv>vmax]=vmax
    if vmin!=0 : MaskV=all(Mv[tri.triangles]<vmin,axis=1) ; Mask0[MaskV]=True ; Mv[Mv<vmin]=vmin
    MS0=CleanTri(tri,1e-10) ; Mask0[MS0]=True
    tri.set_mask( Mask0 )
    print( '=> ',lab,'   :   %.3e  ,  %.3e'%(min(Mv),max(Mv)) ) ; f=0
    if len(OPT)==0 : Field2(tri,Mv,lab,Log,xlim,ylim,vmax,ticks,TICKS,cmap,[],True,name,fs)
    else :
        (fig,ax,cb)=Field2(tri,Mv,lab,Log,xlim,ylim,vmax,ticks,TICKS,cmap,[],False,name,fs)
        if  'LINES' in OPT : #====================> Lines
            print('=> Lines')
            iopt=OPT.index('LINES')
            Vx=OPT[iopt+1]
            if len(ylim)==0 : ylim=[My0,My1]
            for x in Vx : ax.plot( 2*[x],ylim,':w' )
        if 'MIX' in OPT : #====================> Quiver
            print('=> mixture fraction contours')
            iopt=OPT.index('MIX')
            [Mol_m,BC_f,BC_o]=OPT[iopt+1]
            DY=Dic_Y(M,T)
            Zc=ZC(DY,BC_f,BC_o,Mol_m)
            Zh=ZH(DY,BC_f,BC_o,Mol_m)
            Zh_st=Zst(BC_f,BC_o)
            f=ax.tricontour( tri,Zh,levels=[Zh_st],colors='b',linewidths=1 )
        if 'RECIRC' in OPT : #====================> Quiver
            print('=> Recirculation')
            iopt=OPT.index('RECIRC')
            [col]=OPT[iopt+1]
            if 'radial-velocity' in T :
                IUr=FindData('axial-velocity',T)
                Va=M[:,IUr]
            else :
                IUx=FindData('x-velocity',T)
                Va=M[:,IUx]
            f=ax.tricontour( tri,Va,levels=[-1e-3],colors=col,linewidths=1 )
        if 'QUIV' in OPT : #====================> Quiver
            print('=> Quiver')
            iopt=OPT.index('QUIV')
            [Ni,Nj,s]=OPT[iopt+1]
            Vx0=linspace(xlim[0],xlim[1],Ni)
            Vy0=linspace(ylim[0],ylim[1],Nj)
            IUx=FindData('x-velocity',T)
            IUy=FindData('y-velocity',T)
            IUz=FindData('z-velocity',T)
            if   XY : Xi,Xj=Mx,My ; Vi,Vj=M[:,IUx],M[:,IUy]
            elif ZX : Xi,Xj=Mx,Mz ; Vi,Vj=M[:,IUx],M[:,IUz]
            elif YZ : Xi,Xj=My,Mz ; Vi,Vj=M[:,IUy],M[:,IUz]
            f_Vi=mtp.tri.LinearTriInterpolator( tri,Vi )
            f_Vj=mtp.tri.LinearTriInterpolator( tri,Vj )
            (MXi,MXj)=meshgrid(Vx0,Vy0)
            MVi=f_Vi(MXi,MXj)
            MVj=f_Vj(MXi,MXj)
            ax.quiver(MXi,MXj,MVi,MVj,color='w',scale=s)
        if 'PROF' in OPT : #====================> profiles
            print('=> Profiles')
            iopt=OPT.index('PROF')
            [P0,P1,Np,fplot]=OPT[iopt+1] ; Np=int(Np)
            Vi=linspace(P0[0],P1[0],Np)
            Vj=linspace(P0[1],P1[1],Np)
            f=mtp.tri.LinearTriInterpolator( tri,Mv ) ; Prof=f(Vi,Vj)
            fplot(Vi,Vj,Prof,name)
            ax.plot([P0[0],P1[0]],[P0[1],P1[1]],':w')
        if 'Zst' in OPT : #====================> Stoechiometric line 
            print('=> Stoechiometric line')
            iopt=OPT.index('Zst')
            [Zst,cz]=OPT[iopt+1]
            if 'fmean' in T : Yh0=M[:,T.index('fmean')]
            f=ax.tricontour( tri,Yh0,levels=[Zst],colors=cz,linewidths=1 )
        if 'Tiso' in OPT : #====================> Isolines
            print('=> temperature isolines')
            iopt=OPT.index('Tiso')
            Viso=OPT[iopt+1]
            Itr=T.index('temperature')
            f=ax.tricontour( tri,M[:,Itr],levels=Viso,colors='r',linewidths=1 )
        if 'ISO' in OPT : #====================> Isolines
            print('=> isolines')
            iopt=OPT.index('ISO')
            Viso=OPT[iopt+1]
            f=ax.tricontour( tri,Mv,levels=Viso,colors='w',linewidths=1 )
            # for path in f.collections[0].get_paths() :
                #  points=path.vertices
                #  if len(points)>0 : print('=> Xmin : {:.3f} [mm]  ,  Xmax : {:.3f} [mm]  ,  Dx : {:.3f} [mm]'.format(min(points[:,0]),max(points[:,0]),max(points[:,0])-min(points[:,0])))
        if 'INTERP' in OPT : #====================> Interpolation
            print('=> Interpolation')
            f=mtp.tri.LinearTriInterpolator( tri,Mv )
            util.SaveFig(fig,name)
            return(f)
        util.SaveFig(fig,name)
    # return(f)
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
def Field2(tri,F,lab,Log,xlim,ylim,vmax,ticks,cmesh,cmap,CMask,SAVE,name,fs) :
    fig,ax=plt.subplots(figsize=fs)
    # fig.suptitle(lab,fontsize=20)
    ax.set_title(lab,fontsize=20)
    ax.set_aspect('equal')
    if len(xlim) : ax.set_xlim(xlim[0],xlim[1])
    if len(ylim) : ax.set_ylim(ylim[0],ylim[1])
    # if vmax>0 : F[F>vmax]=vmax
    # else      : vmax=max(F) #; print(vmax)
    if not vmax : vmax=max(F)
    # vmax=11
    #=====> Plot
    # ax.triplot(tri,color='red',linewidth=0.1) ; F[:]=0

    if Log : f=ax.tricontourf( tri,F,levels=[ 10**n for n in linspace(0,7,101) ] , cmap=cmap ,vmax=vmax , norm = colors.LogNorm(vmin=F.min(), vmax=F.max()) )
    else   : f=ax.tricontourf( tri,F,levels=int(100) , cmap=cmap ,vmax=vmax )
    #=====> Mask
    # if len(CMask) : ax.fill(CMask[0],CMask[1],facecolor='white',edgecolor='black')
    if len(CMask) : ax.fill(CMask[0],CMask[1],facecolor='white',edgecolor='none')
    #=====> Colorbar
    divider = make_axes_locatable(ax) ; cax = divider.append_axes("right", size="2%", pad=0.25)
    if len(ticks) : cb=fig.colorbar(f,cax=cax,ticks=ticks,extend='both') #,extendrect=False)
    else          : cb=fig.colorbar(f,cax=cax            ,extend='both') #,extendrect=False)
    #=====> Ticks
    if cmesh==0 :
        ax.set_xticks([])
        ax.set_yticks([])
    else :
        xticks=ax.get_xticks()
        yticks=ax.get_yticks()
        ax.set_xticks( xticks )
        ax.set_yticks( yticks )
        ax.set_xticklabels([ '%.0f'%(x*cmesh) for x in xticks ])
        ax.set_yticklabels([ '%.0f'%(y*cmesh) for y in yticks ])
    # cb.set_label(lab,fontsize=20)
    if SAVE :
        util.SaveFig(fig,name)
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
def Probe_read(fprobe) :
    op=open(fprobe)
    for n in range(2) : L0=op.readline()
    L0=op.readline()
    op.closed
    T=[ s.strip()[1:-1] for s in L0[1:-2].split(' ')]
    M=loadtxt(fprobe,skiprows=3,delimiter=' ')
    return({ s:M[:,n] for n,s in enumerate(T) })
#===================================================================
def Probe_plot(fprobe,fplot) :
    D=Probe_read(fprobe)
    It=D['Iteration']
    Var_plot=[ v for v in D.keys() if v!='Iteration' ]
    Nvar=len(Var_plot)
    fig,ax=plt.subplots(figsize=(15,3*Nvar),nrows=Nvar)
    if Nvar==1 : ax=[ax]
    for nv in range(Nvar) :
        v=Var_plot[nv]
        ax[nv].set_title(v)
        ax[nv].plot( It , D[v] , 'k' )
        tick_k(ax[nv])
    fig.tight_layout()
    fig.savefig(fplot)
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
    Nin,Nou=len(Min),len(Mou) ; Nmin=min([Nin,Nou])
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
    if abs(Bl)>1e3*abs(Ml) : tbal='Mean Balance > 1e3 Mean Min'
    else                   : tbal='Balance : {:.1e} [g/s]  ,  RMS : {:.1e} [%]'.format(Bl,Vb*1e2)
    print( 'Min_av : {:.2f}   ,  Min_av_min : {:.2f}   ,  Min_av_max : {:.2f}'.format(Ml,min(M0),max(M0)) )
    print( 'Bal_av : {:.1e}   ,  Bal_av_min : {:.1e}   ,  Bal_av_max : {:.1e}'.format(Bl,min(B0),max(B0)) )
    print( 'RMS m,b : {:.3f}  {:.3f}'.format(Vm*1e2,Vb*1e2) )
    fig_b,ax_b=plt.subplots(figsize=( 15,5),ncols=2)
    if Mth>0 : 
        Er=abs(Ml-Mth)/Mth
        fig_b.suptitle('Teoretical Mf : {:.2f} [g/s]  ,  Error : {:.2f} [%]'.format( Mth,Er*1e2 ),y=0.995)
    ax_b[0].set_title('Mf : {:.2f} [g/s]  ,  RMS : {:.2e} [%]'.format(Ml,Vm*1e2) )
    ax_b[1].set_title(tbal)
    ax_b[0].ticklabel_format( axis='y' , scilimits=(-2,2) ) ; 
    ax_b[1].ticklabel_format( axis='y' , scilimits=(-2,2) ) ; 
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
def Zst(BC_f,BC_o) :
	alp=BC_f['CH4']/BC_f['H2']
	XH2_s =1/(1.5+3*alp)
	XCH4_s=       alp *XH2_s
	XO2_s =(0.5+2*alp)*XH2_s
	Mav_s=Mol_m['H2']*XH2_s+Mol_m['O2']*XO2_s+Mol_m['CH4']*XCH4_s
	YH_s=Mol_m['H']*(4*XCH4_s+2*XH2_s)/Mav_s
	YH_o=Yh(BC_o,Mol_m)
	YH_f=Yh(BC_f,Mol_m)
	# Zst=(YH_s-YH_o)/(YH_f-YH_o) ; print('=> Stoichiometric mixture fraction Zst = %.3f'%(Zst))
	return((YH_s-YH_o)/(YH_f-YH_o))
#===================================================================
def ZH(DY,BC_f,BC_o,Mol_m) :
    Yh_f=Yh(BC_f,Mol_m)
    Yh_o=Yh(BC_o,Mol_m)
    Yh_g=Yh(DY,Mol_m)
    return((Yh_g-Yh_o)/(Yh_f-Yh_o))
#===================================================================
def ZC(DY,BC_f,BC_o,Mol_m) :
    Yc_f=Yc(BC_f,Mol_m)
    Yc_o=Yc(BC_o,Mol_m)
    Yc_g=Yc(DY,Mol_m)
    return((Yc_g-Yc_o)/(Yc_f-Yc_o))
#===================================================================
def Dic_Y(M,T) :
    DY={}
    if 'ch4' in T : DY['CH4']=M[:,T.index('ch4')]
    if 'co2' in T : DY['CO2']=M[:,T.index('co2')]
    if 'co'  in T : DY['CO' ]=M[:,T.index('co' )]
    if 'h2'  in T : DY['H2' ]=M[:,T.index('h2' )]
    if 'o2'  in T : DY['O2' ]=M[:,T.index('o2' )]
    if 'h2o' in T : DY['H2O']=M[:,T.index('h2o')]
    return(DY)
#===================================================================
def Yc(Y,M) : Spe=['CH4','CO2','CO']       ; Id=[1,1,1]   ; return( M['C']*sum([ Id[i]*Y[s]/M[s] for i,s in enumerate(Spe) if s in list(Y.keys()) ],axis=0) )
def Yo(Y,M) : Spe=['O2' ,'CO2','CO','H2O'] ; Id=[2,2,1,1] ; return( M['O']*sum([ Id[i]*Y[s]/M[s] for i,s in enumerate(Spe) if s in list(Y.keys()) ],axis=0) )
def Yh(Y,M) : Spe=['CH4','H2O','H2']       ; Id=[4,2,2]   ; return( M['H']*sum([ Id[i]*Y[s]/M[s] for i,s in enumerate(Spe) if s in list(Y.keys()) ],axis=0) )
#===================================================================
def Umoy(x,V) : return(2*trapz(V*x,x)/(x[-1]-x[0])**2)
#===================================================================
def ConvBC_XY(BC) :
    Spe_bc=[ sp for sp in BC.keys() if sp in Mol_m.keys() ]
    Mm=sum([ BC[sp]*Mol_m[sp] for sp in Spe_bc ])
    for sp in Spe_bc :
        BC[sp]=BC[sp]*Mol_m[sp]/Mm
    return(BC)