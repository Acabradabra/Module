#!/usr/bin/env python3
#-*- coding: utf-8 -*-

from numpy import *

#===================================================================================
def Trans_Coef(L,h1,hN) :
	if h1==hN :
		ap=1
		N=int(round( (L/h1) , 0)+1)
		# N=int(round( (L/h1) , 0))
	else :
		ap=1+(hN-h1)/(L)
		N=int(round( log(hN/h1)/log(ap) , 0)+1)
	return(N,ap)
#===================================================================================
def Normals0() : return([
	array([ 0, 0]),
	array([ 1, 0]),
	array([-1, 0]),
	array([ 0, 1]),
	array([ 0,-1]),
	array([ 1, 1]),
	array([ 1,-1]),
	array([-1, 1]),
	array([-1,-1])
	])
#===================================================================================
def ht_all(H,r,n) : return([ ht(h,r,n) for h in H ])
#===================================================================================
def ht(h0,r,n) : return(h0*(r**n-1)/(r-1))
#===================================================================================
def Corner(Pm,N0,N1,n0,n1,nm,l,h,Nc) :
	Vp=[]
	P0=array(Pm)+Nc*l*array(N0)
	P1=array(Pm)+Nc*l*array(N1)
	Vp.append(geom.add_point( P0+h*array(n0),l ))
	Vp.append(geom.add_point( Pm+h*array(nm),l ))
	Vp.append(geom.add_point( P1+h*array(n1),l ))
	return( Vp )
#===================================================================================
def P_add(geom,Points,Normals,H0in,L0in,rbd,n) :
	N=len(Points)
	if n==0 :
		Ht=N*[0]
	else :
		Ht=ht_all( H0in,rbd,n )
	return([ geom.add_point(Points[p]+Normals[p]*Ht[p],L0in[p]) for p in range(N) ])
#===================================================================================
def L_add(geom,P0,P,P1) :
	N=len(P)
	L=[  geom.add_line(P0  ,P[0]  ) ]
	L+=[ geom.add_line(P[p],P[p+1]) for p in range(N-1) ]
	if P1 :
		L+=[ geom.add_line(P[-1],P1) ]
	return(L)
#===================================================================================
def S_add(geom,VL) :
	ll=geom.add_curve_loop( VL ) ; return(geom.add_plane_surface(ll))
#===================================================================================
def RotationX(a) : return( array([ [1,0,0],[0,cosi(a),sinu(a)],[0,-sinu(a),cosi(a)] ]) )
#===================================================================================
def Rotate(s) : gm.model.geo.rotate( s.dim_tags ,  0,0,0  , 1,0,0 , -2.5*pi/180 )
#===================================================================================
def Revolve(s) : return( geom.revolve( s , (1,0,0) , (0,0,0) , 5*pi/180 , recombine=True , num_layers=1 ) )
#===================================================================================
def MakePolyCells(mesh,name,Nsv,PLOT) :
	from scipy.spatial import Voronoi, voronoi_plot_2d
	import matplotlib.pyplot as plt
	import matplotlib

	if PLOT : fig_v,ax_v=plt.subplots()

	#==============> Mesh Points
	VP0=mesh.points[:,:2] ; NP0=len(VP0[:,0]) # original points
	Vor=Voronoi(VP0)                          # Voronoi
	Rnp0=array(range(NP0))                    # Point indice array

	#==============> Boundary points
	L0=mesh.get_cells_type('line') # Lines
	MedL0=mean(VP0[L0],axis=1)     # Center line coordinates
	BDP0=list(L0[:,0])             # First points of each lines
	Nbp0=len(BDP0)

	#==============> Voronoi elements
	Nc=Vor.npoints  # number of cells 
	P0=Vor.points   # Original points
	P1=Vor.vertices                ; NP1=len(P1) # New points
	P2=array(list(P1)+list(MedL0)) ; NP2=len(P2) # New points with boundary median
	C1=Vor.regions                  # Cells
	PR=Vor.point_region             # Voronoi region for each input point
	if len(C1[0])==0 : C1=C1[1:]

	#==============> Plotting
	if PLOT :
		# fig_v = voronoi_plot_2d(Vor)
		ax_v.plot(    VP0[:,0],  VP0[:,1] , 'ok' )
		ax_v.plot(  MedL0[:,0],MedL0[:,1] , 'oy' )
		ax_v.plot(     P0[:,0],   P0[:,1] , '.g' )
		ax_v.plot(     P2[:,0],   P2[:,1] , '.k' )
		ax_v.plot(  P0[BDP0,0],P0[BDP0,1] , '.b' )

	#==============> Boundary cells
	BDC=[]
	BDV=[]
	for p in BDP0 :
		if PLOT : ax_v.plot( P0[p,0],P0[p,1],'xr' )
		ic=PR[p]-1
		c=C1[ic]
		c2=list(copy(c))

		#=====> Insert boundary
		n0=c2.index(-1)
		i0=BDP0.index(p) # p is the index of the boundary points in P0
		i1=NP1+i0        # Index of the boundary point in P2 just after p which is in P0
		if i0>0 : i2=i1-1
		else    : i2=NP2-1
		#=====> Choose order
		pb=c2[n0-1]
		N1=linalg.norm( P2[pb]-P2[i1] )
		N2=linalg.norm( P2[pb]-P2[i2] )
		if N1<N2 :
			c2[n0]=i2
			c2.insert(n0,i1)
		else :
			c2[n0]=i1
			c2.insert(n0,i2)

		#=====> Check corner
		V1=P2[i1,:]-P0[p,:] ; NV1=linalg.norm(V1)
		V2=P2[i2,:]-P0[p,:] ; NV2=linalg.norm(V2) 
		scal=sum(V1*V2)/(NV1*NV2) #; print(scal)
		if scal>-1 :
			P2=append(P2,[P0[p,:]],axis=0)
			if PLOT : ax_v.plot( P2[-1,0],P2[-1,1],'+g' )
			i3=len(P2)-1
			c2.insert(n0+1,i3)
			BDV.append([i2,i3])
			BDV.append([i3,i1])
			BDC.append(ic)
			BDC.append(ic)
		else :
			BDV.append([i2,i1])
			BDC.append(ic)

		# C1.append(c2)
		C1[ic]=c2
		c2+=c2[:1]
		if PLOT : ax_v.plot( P2[c2,0],P2[c2,1],'--r' )

	#==============> Smoothing
	#=====> Conections detecion
	RV=array(Vor.ridge_vertices) # new points in vertices
	RP=array(Vor.ridge_points)   # old points on both sides
	SRV1=(RV[:,0]==-1)
	SRV2=(RV[:,1]==-1)
	Cnct=[]
	for p in range(NP1) :
		Sc1=(RV[:,0]==p)
		Sc2=(RV[:,1]==p)
		Vc=list(  RV[Sc1,1])
		Vc.extend(RV[Sc2,0])
		if -1 in Vc :
			Rbd=Sc1*SRV2+Sc2*SRV1 # Select ridge with -1 and p
			Pr1,Pr2=RP[Rbd,:][0]
			ir1=BDP0.index(Pr1)
			ir2=BDP0.index(Pr2)
			i0=min(ir1,ir2) # p is the index of the boundary points in P0
			if i0==0 and max(ir1,ir2)>1 : i1=NP2-1
			else                        : i1=NP1+i0 # Index of the boundary point in P2 just after p which is in P0
			# Vc=[ p for p in Vc if p>-1 ]
			Vc[Vc.index(-1)]=i1    #; print(NP1,i0,i1,ir1,ir2)
			Sm1 = RV[Rbd,:][0]==-1 #; print(Sm1,i1,RV[Rbd,Sm1])
			RV[Rbd,Sm1]=i1
			# print(P2[RV[Rbd,:][0],0])
			# print(P2[RV[Rbd,:][0],1])
			if PLOT : Line=RV[Rbd,:][0] ; ax_v.plot( P2[Line,0],P2[Line,1],'-b' )
		Cnct.append(Vc)
	#=====> Smoothing
	for n in range(Nsv) :
		for p in range(NP1) :
			P2[p,:]=mean(P2[Cnct[p],:],axis=0)
	# Ncon=100
	# for Ncon in range(NP1) :
	# 	ax_v.plot( P2[Ncon,0],  P2[Ncon,1] , '*r' )
	# 	for p in Cnct[Ncon] :
	# 		ax_v.plot( P2[[Ncon,p],0],  P2[[Ncon,p],1] , 'r' )

	#==============> Ploting
	if PLOT :
		for n in range(Nc) : 
			c=C1[n]
			if len(c)==0 : print(c)
			elif -1 in c        : pass
			elif len(c)==6  :
				c2=c+[c[0]]
				ax_v.plot( P2[c2,0],P2[c2,1],'k' )
			elif len(c)==7  :
				c2=c+[c[0]]
				ax_v.plot( P2[c2,0],P2[c2,1],'g' )
			elif len(c)==5  :
				c2=c+[c[0]]
				ax_v.plot( P2[c2,0],P2[c2,1],'r' )
			elif len(c)==4  :
				c2=c+[c[0]]
				ax_v.plot( P2[c2,0],P2[c2,1],'b' )
		fig_v.savefig(name,dpi=1e3) ; plt.close(fig_v)

	#==============> Center
	Cc =array([ mean(P2[c],axis=0) for c in C1 ])
	Rc =array([ mean(P2[c],axis=0) for c in RV ])
	BDc=array([ mean(P2[c],axis=0) for c in BDV ])

	class Poly :
		def __init__(self) :
			self.points=P2
			self.cells=C1
			self.Np_c=NP1
			self.Np_t=len(P2)
			self.Nf_c=len(RV)
			self.Nf_b=len(BDC)
			self.ridge_vert=RV
			self.ridge_cells=PR[RP]
			self.center_cells=Cc
			self.center_ridge=Rc
			self.center_ridge_bd=BDc
			self.boundary_vert =array(BDV)
			self.boundary_cells=array(BDC)
	return(Poly())
#===================================================================================
def PlotPoly(Poly,name) :
	import matplotlib.pyplot as plt
	import matplotlib
	Pvo=Poly.points
	Cvo=Poly.cells
	Nbd=Poly.Np_c
	RV =Poly.ridge_vert
	RC =Poly.ridge_cells
	Cc =Poly.center_cells
	BDV=Poly.boundary_vert
	BDC=Poly.boundary_cells
	fig_v,ax_v=plt.subplots()
	ax_v.plot(Pvo[:Nbd,0],Pvo[:Nbd,1],'.b')
	ax_v.plot(Pvo[Nbd:,0],Pvo[Nbd:,1],'.g')
	for c in Cvo :
		c2=c+[c[0]]
		ax_v.plot(           Pvo[c2,0],Pvo[c2 ,1],  'k' )
	for r in RV  : ax_v.plot( Pvo[r ,0],Pvo[r  ,1], ':b' )
	for r in RC  : ax_v.plot( Cc[r-1,0],Cc [r-1,1],'--b' )
	for c in Cc  : ax_v.plot(      c[0],      c[1], '*r' )
	for r in range(len(BDV)) : 
		rv=BDV[r]
		Crl=mean(Pvo[rv,:],axis=0)
		Ccl=Cc[BDC[r]]
		ax_v.plot( Pvo[rv,0],Pvo[rv,1], ':g' )
		ax_v.plot( [Ccl[0],Crl[0]],[Ccl[1],Crl[1]], '--g' )
	fig_v.savefig(name,dpi=1e3) ; plt.close(fig_v)
#===================================================================================
def WriteAnsys(filename, mesh,Poly, binary=True):
	with open(filename, "wb") as fh:
		RV =Poly.ridge_vert
		RC =Poly.ridge_cells
		Nf_c=Poly.Nf_c
		Nf_b=Poly.Nf_b
		BDV=Poly.boundary_vert
		BDC=Poly.boundary_cells

		#==========> Header
		fh.write(f'(1 "meshio {__version__}")\n'.encode())

		#============================================================> Dimension
		fh.write(f'\n(0 "Dimension")\n'.encode())
		#============================================================
		num_points, dim = mesh.points.shape
		if dim not in [2, 3]:
			raise WriteError(f"Can only write dimension 2, 3, got {dim}.")
		fh.write((f"(2 {dim})\n").encode())

		#============================================================> General description
		fh.write(f'\n(0 "Nodes,Cells,Faces")\n'.encode())
		#============================================================
		#==========> total number of cells
		total_num_cells = sum(len(c) for c in mesh.cells)
		fh.write((f"(12 (0 1 {total_num_cells:x} 0))\n").encode())
		#==========> total number of faces
		total_num_faces = Nf_c+Nf_b
		fh.write((f"(13 (0 1 {total_num_faces:x} 0))\n").encode())
		#==========> total number of Nodes
		first_node_index = 1
		fh.write((f"(10 (0 {first_node_index:x} {num_points:x} 0 2))\n").encode())

		#============================================================> Write faces
		fh.write(f'\n(0 "Cells")\n'.encode())
		#============================================================
		fh.write((f"(12 (4 1 {total_num_cells:x} 1 0))\n").encode())

		#============================================================> Write faces
		fh.write(f'\n(0 "Faces Internal")\n'.encode())
		#============================================================
		key = "3013" if binary else "13"
		fh.write(f"({key} (2 {1:x} {Nf_c:x} 2 {2:x})(\n".encode() )
		for n in range(Nf_c) : 
			fh.write(f"{RV[n,0]:x} {RV[n,1]:x} {RC[n,0]:x} {RC[n,1]:x}".encode() )
			if n<Nf_c-1 : fh.write(f"\n".encode() )
			else        : fh.write(f"))\n".encode() )
		
		#============================================================> Write faces
		fh.write(f'\n(0 "Faces Boundary")\n'.encode())
		#============================================================
		key = "3013" if binary else "13"
		fh.write(f"({key} (3 {Nf_c+1:x} {total_num_faces:x} 3 {2:x})(\n".encode() )
		for n in range(Nf_b) : 
			fh.write(f"{BDV[n,0]:x} {BDV[n,1]:x} {BDC[n]:x} {0:x}".encode() )
			if n<Nf_b-1 : fh.write(f"\n".encode() )
			else        : fh.write(f"))\n".encode() )

		#============================================================> Write nodes
		fh.write(f'\n(0 "Nodes")\n'.encode())
		#============================================================
		key = "3010" if binary else "10"
		fh.write(f"({key} (1 {first_node_index:x} {num_points:x} 1 {dim:x})(\n".encode() )
		if binary:
			mesh.points.tofile(fh)
			fh.write(b"\n)")
			fh.write(b"End of Binary Section 3010)\n")
		else:
			savetxt(fh, mesh.points, fmt="%.16e")
			fh.write(b"))\n")