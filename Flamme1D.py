import cantera as ct
import numpy as np
import Utilities as util
from scipy.interpolate import interp1d

(plt,mtp)=util.Plot0()

#---------------------------------------------------------------------
def Adia0D(phi,X0,Ti,Pi,INox,Tol,gas):
	

	XY='HP'
	solv='auto'
	MaxStep=10000
	MaxIter=1000

	W=gas.molecular_weights
	a=(1-X0)/X0

	gas.X={'H2':2,'O2': 1/phi ,'N2': a/phi}
	mix=ct.Mixture( [(gas,1.0)] )
	mix.T=Ti
	mix.P=Pi
	mix.equilibrate(XY,solver=solv,max_steps=MaxStep,max_iter=MaxIter,
		estimate_equil=1,log_level=0,rtol=Tol)

	Tad=mix.T
	Xmol=mix.species_moles
	Wt=sum( [ W[j]*Xmol[j] for j in range(gas.n_species) ] )
	XNox=[ Xmol[j] for j in INox ]
	YNox=[ Xmol[j]*W[j]/Wt for j in INox ]

	return(Tad,XNox,YNox)

#---------------------------------------------------------------------
def Adia1D(phi,X0,Ti,Pi,grid0,Name,INox,Tol,Tstep,Rafcrit,gas,view):

	a=(1-X0)/X0
	Xi={'H2':2,'O2': 1/phi ,'N2': a/phi}

	########## --------------------> Gas state + Flame
	gas.TPX=Ti,Pi,Xi
	f=ct.FreeFlame(gas,grid0)

	try : f.restore(filename=Name,loglevel=0)
	except : pass

	f.inlet.T=Ti
	f.inlet.X=Xi

	f.flame.set_steady_tolerances( default=Tol )
	f.flame.set_transient_tolerances( default=[ x*10 for x in Tol ])
	f.max_grid_points=5000
	f.set_time_step(Tstep,[2,5,10,20,50,80,130,200,500,1000])
	
	if Rafcrit==4:
	        f.set_refine_criteria(ratio=5.0, slope=0.005, curve=0.001, prune=0.00075)
	elif Rafcrit==3:
	        f.set_refine_criteria(ratio=2.0, slope=0.05, curve=0.005,  prune=0      )
	elif Rafcrit==2:
	        f.set_refine_criteria(ratio=5.0, slope=0.05, curve=0.01,   prune=0.0075 )
	elif Rafcrit==1:
	        f.set_refine_criteria(ratio=5.0, slope=0.1,  curve=0.05,   prune=0.025  )
	elif Rafcrit==0:
	        f.set_refine_criteria(ratio=5.0, slope=0.5,  curve=0.1,    prune=0.075  )

	f.transport_model='Mix'
	# f.soret_enabled=True
	f.energy_enabled=True
	f.radiation_enabled=False

	conv=True
	try : f.solve(loglevel=view,refine_grid=True,auto=False)
	except : conv=False
	
	if conv: f.save(filename=Name,loglevel=0)

	########## --------------------> Extract Data
	XNox,YNox=[],[]
	for i in INox:
		XNox.append(f.X[i])
		YNox.append(f.Y[i])

	hpr=f.heat_production_rates
	hrr=f.heat_release_rate
	Grid=f.grid

	Nreac=len(gas.reaction_equations())
	Ngrid=len(f.grid)
	MaxHeat,IntH=[],[]
	for j in range(Nreac):
		MaxHeat.append( np.max( [ abs(x) for x in hpr[j] ] ) )
		IntH.append( sum( [ 0.5*( hpr[j][n+1]+hpr[j][n] )*( Grid[n+1]-Grid[n] ) for n in range(Ngrid-1) ] ) )
	Q=sum( [ 0.5*( hrr[n+1]+hrr[n] )*( Grid[n+1]-Grid[n] ) for n in range(Ngrid-1) ] )

	return(f,hpr,hrr,Grid,conv,XNox,YNox,MaxHeat,IntH,Q,Ngrid)

#---------------------------------------------------------------------
#def Adia1D_0(phi,X0,Ti,Pi,grid0,NameIn,NameOut,Tol,Tstep,Rafcrit,gas,Ray,view):
def Adia1D_0(Xi,Ti,Pi,grid0,NameIn,NameOut,Tol,Tstep,Rafcrit,gas,Ray,view):

	# a=(1-X0)/X0
	# Xi={'H2':2,'O2': 1/phi ,'N2': a/phi}

	# print('Adia-A')

	########## --------------------> Gas state + Flame
	gas.TPX=Ti,Pi,Xi
	f=ct.FreeFlame(gas,grid0)

	# print('Adia-B')

	try : f.restore(filename=NameIn,loglevel=0)# ; print('=> Adia restore')
	except : pass

	f.inlet.T=Ti
	f.inlet.X=Xi

	f.flame.set_steady_tolerances( default=Tol )
	#f.flame.set_transient_tolerances( default=[ x*10 for x in Tol ])
	f.flame.set_transient_tolerances( default=[ x/10 for x in Tol ])
	#f.flame.set_transient_tolerances( default=Tol )
	f.max_grid_points=50000
	f.set_time_step(Tstep,[2,5,10,20,50,100,200,500,1000,1e5,1e6,1e7])
	
	if Rafcrit==4:
	        f.set_refine_criteria(ratio=5.0, slope=0.005, curve=0.001, prune=0.00075)
	elif Rafcrit==3:
	        f.set_refine_criteria(ratio=2.0, slope=0.05, curve=0.005,  prune=0      )
	elif Rafcrit==2.5:
	        f.set_refine_criteria(ratio=2.0, slope=0.05, curve=0.01,   prune=0.0075 )
	elif Rafcrit==2:
	        f.set_refine_criteria(ratio=5.0, slope=0.05, curve=0.01,   prune=0.0075 )
	elif Rafcrit==1:
	        f.set_refine_criteria(ratio=5.0, slope=0.1,  curve=0.05,   prune=0.025  )
	elif Rafcrit==0:
	        f.set_refine_criteria(ratio=5.0, slope=0.5,  curve=0.1,    prune=0.075  )

	f.transport_model='Mix'
	f.energy_enabled=True
	f.radiation_enabled=Ray

	# print('Adia-C')

	conv=True
	try : f.solve(loglevel=view,refine_grid=Rafcrit,auto=False)
	except : conv=False
	
	# print('Adia-D')

	if conv: f.save(filename=NameOut,loglevel=0) # ; print('=> Adia Save')

	return(f,conv)

#---------------------------------------------------------------------
def ColdFlow(Xi,Ti,Pi,grid0,NameIn,NameOut,Tol,Tstep,Rafcrit,gas,Ray,view):
	gas.TPX=Ti,Pi,Xi
	f=ct.FreeFlame(gas,grid0)
	f.inlet.T=Ti
	f.inlet.X=Xi
	f.transport_model='Mix'
	f.energy_enabled=True
	f.radiation_enabled=Ray
	return(f,True)
#---------------------------------------------------------------------
def Burner1D(phi,X0,mdot,Ti,Pi,grid0,Name,INox,Tol,Tstep,Rafcrit,gas,view):
	
	a=(1-X0)/X0
	Xi={'H2':2,'O2': 1/phi ,'N2': a/phi}

	
	########## --------------------> Gas state + Flame
	gas.TPX=Ti,Pi,Xi
	f=ct.BurnerFlame(gas,grid0)

	try: f.restore(filename=Name,loglevel=0)
	except: pass

	f.burner.T=Ti
	f.burner.X=Xi
	f.burner.mdot=mdot
	
	f.flame.set_steady_tolerances( default=Tol)
	f.flame.set_transient_tolerances( default=[ x*10 for x in Tol ])
	f.max_grid_points=5000
	
	f.set_time_step(Tstep,[2,5,10,20,50,80,130,200,500,1000])
	if Rafcrit==4:
	        f.set_refine_criteria(ratio=5.0, slope=0.005, curve=0.001, prune=0.00075)
	elif Rafcrit==3:
	        f.set_refine_criteria(ratio=2.0, slope=0.05, curve=0.005,  prune=0 )
	elif Rafcrit==2:
	        f.set_refine_criteria(ratio=5.0, slope=0.05, curve=0.01,   prune=0.0075 )
	elif Rafcrit==1:
	        f.set_refine_criteria(ratio=5.0, slope=0.1,  curve=0.05,   prune=0.025  )
	elif Rafcrit==0:
	        f.set_refine_criteria(ratio=5.0, slope=0.5,  curve=0.1,    prune=0.075  )
	
	f.transport_model='Mix'
	f.energy_enabled=True
	f.radiation_enabled=False
	
	conv=True
	try : f.solve(loglevel=view,refine_grid=True,auto=False)
	except : conv=False
	
	if conv: f.save(filename=Name,loglevel=0)

	########## --------------------> Extract Data
	# XNox,YNox=[],[]
	# for i in INox:
	# 	XNox.append(f.X[i])
	# 	YNox.append(f.Y[i])
	
	hpr=f.heat_production_rates
	hrr=f.heat_release_rate
	Grid=f.grid

	Nreac=len(gas.reaction_equations())
	Ngrid=len(f.grid)
	MaxHeat,IntH=[],[]
	for j in range(Nreac):
		MaxHeat.append( np.max( [ abs(x) for x in hpr[j] ] ) )
		IntH.append( sum( [ 0.5*( hpr[j][n+1]+hpr[j][n] )*( f.grid[n+1]-f.grid[n] ) for n in range(Ngrid-1) ] ) )
	Q=sum( [ 0.5*( hrr[n+1]+hrr[n] )*( Grid[n+1]-Grid[n] ) for n in range(Ngrid-1) ] )

	# return(f,hpr,hrr,Grid,conv,XNox,YNox,MaxHeat,IntH,Q,Ngrid)
	return(f,hpr,hrr,Grid,conv,MaxHeat,IntH,Q,Ngrid)

#---------------------------------------------------------------------
def Burner1D_0(Xi,mdot,Ti,Pi,grid0,NameIn,NameOut,Tol,Tstep,Rafcrit,gas,Ray,view):
	
	# a=(1-X0)/X0
	# Xi={'H2':2,'O2': 1/phi ,'N2': a/phi}

	
	########## --------------------> Gas state + Flame
	gas.TPX=Ti,Pi,Xi
	f=ct.BurnerFlame(gas,grid0)

	try: f.restore(filename=NameIn,loglevel=0)
	except: pass

	f.burner.T=Ti
	f.burner.X=Xi
	f.burner.mdot=mdot
	
	f.flame.set_steady_tolerances( default=Tol)
	f.flame.set_transient_tolerances( default=[ x*10 for x in Tol ])
	f.max_grid_points=5000
	
	f.set_time_step(Tstep,[2,5,10,20,50,80,130,200,500,1000,5000,10000])
	if Rafcrit==4:
	        f.set_refine_criteria(ratio=5.0, slope=0.005, curve=0.001, prune=0.00075)
	elif Rafcrit==3:
	        # f.set_refine_criteria(ratio=2.0, slope=0.05, curve=0.005,  prune=0 )
	        f.set_refine_criteria(ratio=2.0, slope=0.025, curve=0.005,  prune=0.0025 )
	elif Rafcrit==2:
	        f.set_refine_criteria(ratio=5.0, slope=0.05, curve=0.01,   prune=0.0075 )
	elif Rafcrit==1:
	        f.set_refine_criteria(ratio=5.0, slope=0.1,  curve=0.05,   prune=0.025  )
	elif Rafcrit==0:
	        f.set_refine_criteria(ratio=5.0, slope=0.5,  curve=0.1,    prune=0.075  )
	
	f.transport_model='Mix'
	f.energy_enabled=True
	f.radiation_enabled=Ray
	
	conv=True
	try : f.solve(loglevel=view,refine_grid=True,auto=False)
	except : conv=False
	
	if conv: f.save(filename=NameOut,loglevel=0)

	return(f,conv)

#---------------------------------------------------------------------
def CounterFlow(Xi,Mdot,Ti,Pi,grid0,NameIn,NameOut,Tol,Tstep,Rafcrit,gas,Ray,view):

	########### --------------------> Gas state + Flame
	gas.TPX=Ti,Pi,Xi
	f=ct.CounterflowPremixedFlame(gas,grid0)
	
	try: f.restore(filename=NameIn,loglevel=0)
	except: pass
	
	f.reactants.mdot = Mdot[0]
	f.products.mdot  = Mdot[1]
	
	f.flame.set_steady_tolerances( default=Tol)
	f.flame.set_transient_tolerances( default=[ x*10 for x in Tol ])
	f.max_grid_points=5000
	
	f.set_time_step(Tstep,[2,5,10,20,50,80,130,200,500,1000,5000,10000])
	if Rafcrit==4:
	        f.set_refine_criteria(ratio=5.0, slope=0.005, curve=0.001, prune=0.00075)
	elif Rafcrit==3:
	        # f.set_refine_criteria(ratio=2.0, slope=0.05, curve=0.005,  prune=0 )
	        f.set_refine_criteria(ratio=2.0, slope=0.025, curve=0.005,  prune=0.0025 )
	elif Rafcrit==2:
	        f.set_refine_criteria(ratio=5.0, slope=0.05, curve=0.01,   prune=0.0075 )
	elif Rafcrit==1:
	        f.set_refine_criteria(ratio=5.0, slope=0.1,  curve=0.05,   prune=0.025  )
	elif Rafcrit==0:
	        f.set_refine_criteria(ratio=5.0, slope=0.5,  curve=0.1,    prune=0.075  )
	
	f.transport_model='Mix'
	f.energy_enabled=True
	f.radiation_enabled=Ray
	
	conv=True
	try : f.solve(loglevel=view,refine_grid=True,auto=False)
	except : conv=False
	
	if conv: f.save(filename=NameOut,loglevel=0)
	
	return(f,conv)
#---------------------------------------------------------------------
def TwinFlames(Xi,mdot,Ti,Pi,grid0,NameIn,NameOut,Tol,Tstep,Rafcrit,gas,Ray,view):

	########### --------------------> Gas state + Flame
	gas.TPX=Ti,Pi,Xi
	f=ct.CounterflowTwinPremixedFlame(gas,grid0)
	
	try: f.restore(filename=NameIn,loglevel=0)
	except: pass
	
	f.reactants.mdot = mdot
	
	f.flame.set_steady_tolerances( default=Tol)
	f.flame.set_transient_tolerances( default=[ x*10 for x in Tol ])
	f.max_grid_points=5000
	
	f.set_time_step(Tstep,[2,5,10,20,50,80,130,200,500,1000,5000,10000])
	if Rafcrit==4:
	        f.set_refine_criteria(ratio=5.0, slope=0.005, curve=0.001, prune=0.00075)
	elif Rafcrit==3:
	        # f.set_refine_criteria(ratio=2.0, slope=0.05, curve=0.005,  prune=0 )
	        f.set_refine_criteria(ratio=2.0, slope=0.025, curve=0.005,  prune=0.0025 )
	elif Rafcrit==2:
	        f.set_refine_criteria(ratio=5.0, slope=0.05, curve=0.01,   prune=0.0075 )
	elif Rafcrit==1:
	        f.set_refine_criteria(ratio=5.0, slope=0.1,  curve=0.05,   prune=0.025  )
	elif Rafcrit==0:
	        f.set_refine_criteria(ratio=5.0, slope=0.5,  curve=0.1,    prune=0.075  )
	
	f.transport_model='Mix'
	f.energy_enabled=True
	f.radiation_enabled=Ray
	
	conv=True
	try : f.solve(loglevel=view,refine_grid=True,auto=False)
	except : conv=False
	
	if conv: f.save(filename=NameOut,loglevel=0)
	
	return(f,conv)

#---------------------------------------------------------------------
def LocateNOx(Spe,Nspe):
	INOx=[]
	NOx=[]
	for i in range(Nspe):
	        spe=Spe[i]
	        spe2=set(spe)
	        spe2=[ x for x in spe2 if not x.isnumeric() ]
	        spe2=sorted(spe2)
	        if len(spe2)==2 and spe2[0]=='N' and spe2[1]=='O':
	                INOx.append(Spe.index(spe))
	                NOx.append(spe)
	NNOx=len(NOx)

	return(NOx,NNOx,INOx)

#---------------------------------------------------------------------
def ExtractNox(f,INox):
	XNox,YNox=[],[]
	for i in INox:
		XNox.append(f.X[i])
		YNox.append(f.Y[i])

	return(XNox,YNox)

#---------------------------------------------------------------------
def ExtractHeat(f,Nreac):
	hpr=f.heat_production_rates  # [W/m^3]
	hrr=f.heat_release_rate      # [W/m^3]
	Grid=f.grid
	Ngrid=len(Grid)

	MH,IH,XH=[],[],[]
	for j in range(Nreac):
		MH.append( np.max( [ abs(x) for x in hpr[j] ] ) )
		XH.append( [ Grid[n] for n in range(Ngrid) if abs(hpr[j][n])==MH[-1] ][0] )
		IH.append( sum( [ 0.5*( hpr[j][n+1]+hpr[j][n] )*( Grid[n+1]-Grid[n] ) for n in range(Ngrid-1) ] ) )

	Q=sum( [ 0.5*( hrr[n+1]+hrr[n] )*( Grid[n+1]-Grid[n] ) for n in range(Ngrid-1) ] )
	MQ=np.max(hrr)
	VXQ=[ Grid[n] for n in range(Ngrid) if hrr[n]==MQ ]
	if len(VXQ)>1 : print('\033[41m ==> Several XQ  \033[0m')
	XQ=VXQ[0]

	od,stri=5,1e-8
	if XQ==0:
		print('\033[1m XQ=0 \033[0m')
		(DQ,Test,detD,CD,trace)=util.GradDifFin(od,stri,Grid,Ngrid,hrr)
		p=[ n for n in range(Ngrid-od-1) if DQ[n]>0 and DQ[n+1]<0 ]
		
		if p: 
			XQ=Grid[p[0]]
			MQ=hrr [p[0]]
		else:
			print('\033[42m no Df=0  \033[0m')
			(D2Q,Test,detD,CD,trace)=util.GradDifFin(od,stri,Grid[:-od],Ngrid-od,DQ)
			# absDQ=[abs(x) for x in DQ]
			# minDQ=min(absDQ)
			# p=[ n for n in range(Ngrid-od-1) if absDQ[n]==minDQ ]
			p=[ n for n in range(Ngrid-2*od-1) if D2Q[n]>0 and D2Q[n+1]<0 ]
			trace=2 

			if p: 
				XQ=Grid[p[0]]
				MQ=hrr [p[0]]
			else: print('\033[41m no D2f=0  \033[0m')

	return(hpr,hrr,MH,XH,IH,MQ,XQ,Q,Grid,Ngrid)

#---------------------------------------------------------------------
def ExtractHeat2(f,Nreac):
	hpr=f.heat_production_rates # [W/m^3]
	hrr=f.heat_release_rate # [W/m^3]
	Grid=f.grid
	Ngrid=len(Grid)

	MH,IH,XH=[],[],[]
	for j in range(Nreac):
		MH.append( np.max( [ abs(x) for x in hpr[j] ] ) )
		XH.append( [ Grid[n] for n in range(Ngrid) if abs(hpr[j][n])==MH[-1] ][0] )
		IH.append( sum( [ 0.5*( hpr[j][n+1]+hpr[j][n] )*( Grid[n+1]-Grid[n] ) for n in range(Ngrid-1) ] ) )
	Q=sum( [ 0.5*( hrr[n+1]+hrr[n] )*( Grid[n+1]-Grid[n] ) for n in range(Ngrid-1) ] )
	# absHrr=[abs(x) for x in hrr]
	MQ=np.max(hrr)

	p=[ n for n in range(Ngrid) if hrr[n]==MQ ][0]
	XQ=Grid[p]
	# XQ=[ Grid[n] for n in range(Ngrid) if hrr[n]==MQ ][0]

	od,stri=5,1e-8

	if XQ==0:
		print('\033[1m XQ=0 \033[0m')
		trace=1
		(DQ,Test,detD,CD,trace)=util.GradDifFin(od,stri,Grid,Ngrid,hrr)
		p=[ n for n in range(Ngrid-od-1) if DQ[n]>0 and DQ[n+1]<0 ]
		
		if p: 
			p=p[0]
			XQ=Grid[p]
			MQ=hrr [p]
		else:
			print('\033[42m no Df=0  \033[0m')
			(D2Q,Test,detD,CD,trace)=util.GradDifFin(od,stri,Grid[:-od],Ngrid-od,DQ)
			# absDQ=[abs(x) for x in DQ]
			# minDQ=min(absDQ)
			# p=[ n for n in range(Ngrid-od-1) if absDQ[n]==minDQ ]
			p=[ n for n in range(Ngrid-2*od-1) if D2Q[n]>0 and D2Q[n+1]<0 ]
			trace=2 

			if p:
				p=p[0]
				XQ=Grid[p]
				MQ=hrr [p]
			else: print('\033[41m no D2f=0 \033[0m')

	return(hpr,hrr,MH,XH,IH,MQ,XQ,Q,Grid,Ngrid,p)

#---------------------------------------------------------------------
def ExtractDif(f):
	cpm     =f.cp_mass                # [J/kg/K]
	cpm_spe =f.partial_molar_cp       # [J/kmol/K]

	lamb    =f.thermal_conductivity   # [W/m/K]
	diffth  =f.thermal_diff_coeffs    # [kg/m/s]

	difm    =f.mix_diff_coeffs_mass   # [m^2/s]

	return(cpm,cpm_spe,lamb,diffth,difm)
#---------------------------------------------------------------------
def OldCSV(f,name,spe,Idou,view):

    if view>0: print('\n=====> Starting OldCsv write \n \n => Reading data')

    Spe=f.gas.species_names ; Spe=[ Spe[i] for i in Idou ]
    NSpe=len(Spe)
    NReac=len(f.gas.reaction_equations())

    grid=f.grid
    U   =f.u
    T   =f.T
    rho =f.density_mass
    P   =f.P
    if   spe=='X': Zspe=f.X
    elif spe=='Y': Zspe=f.Y
    else : print('\n Wrong spe (X,Y) =>',spe)

    ProdRate=f.net_production_rates
   
    FwRate =f.forward_rates_of_progress
    RvRate =f.reverse_rates_of_progress
    NetRate=f.net_rates_of_progress
   
    HeatR  =f.heat_release_rate

    Ngrid=len(grid)
    Sl   =U[0]

    if view>0: print('\n => Data readed \n \n==> Writing started')
    
    Line1=Spe
    Line1.extend(['w_'+x for x in Spe ])
    Line1.extend(['FwRate_'+str(x+1) for x in range(NReac)])
    Line1.extend(['RvRate_'+str(x+1) for x in range(NReac)])
    Line1.extend(['NetRate_'+str(x+1) for x in range(NReac)])
    Line1.append('Heat_release')

    Line=[[] for i in range(Ngrid)]
    for i in range(Ngrid):
        Line[i]=[grid[i],U[i],T[i],rho[i],P]

        Z=[Zspe[j][i] for j in Idou] ; SumZ=sum(Z) ; Z=[z/SumZ for z in Z]
        Line[i].extend(Z)
 
        Line[i].extend([ProdRate[j][i] for j in Idou])
 
        Line[i].extend([FwRate [j][i] for j in range(NReac)])
        Line[i].extend([RvRate [j][i] for j in range(NReac)])
        Line[i].extend([NetRate[j][i] for j in range(NReac)])
 
        Line[i].append(HeatR[i])

    with open(name+'.csv','w') as file:
        file.write('Grid Points: ,{0}, Sl= ,{1:.12f}\n'.format(Ngrid,Sl))
        file.write('x_axis,u,Temperature,rho,Pressure'+util.arrayToStr(Line1,',')+'\n')
        for l in Line:
            file.write(util.printround2(l,'',12,',','g')+'\n')
    file.closed

    if view>0: print('\n==> Writing done')
#---------------------------------------------------------------------
def ExtractSubmech(schem,view):
    all_spe=ct.Species.listFromFile(schem+'.cti')
    Spe=[ s for s in all_spe 
    if  not 'C'  in s.composition 
    and not 'Ar' in s.composition 
    and not 'He' in s.composition  ]
    Spe_name=[ s.name for s in Spe ]
    Nspe=len(Spe_name)

    INOx=[]
    NOx=[]
    for i in range(Nspe):
            spe=Spe_name[i]
            spe2=set(spe)
            spe2=[ x for x in spe2 if not x.isnumeric() ]
            spe2=sorted(spe2)
            if len(spe2)==2 and spe2[0]=='N' and spe2[1]=='O':
                    INOx.append(Spe_name.index(spe))
                    NOx.append(spe)
    NNOx=len(NOx)

    if view>0:
        print('\n=====> Species selected',schem,'\n')
        print(Spe_name)
        print('Nspe',Nspe,'\n')
        print(NOx)
        print(INOx)
        print('NNOx',NNOx,'\n')


    
    all_reac=ct.Reaction.listFromFile(schem+'.cti')
    Reac=[ r for r in all_reac 
    if all(Reactants in Spe_name for Reactants in r.reactants) 
    and all(Products in Spe_name for Products  in r.products) ]
    
    gas=ct.Solution(thermo='IdealGas',kinetics='GasKinetics',
        species=Spe,reactions=Reac)
    Reac=gas.reaction_equations()
    Nreac=len(Reac)

    return(gas,Spe_name,Nspe,NOx,INOx,NNOx,Reac,Nreac)
#---------------------------------------------------------------------
def Thick(X,T,N) :
	od,stri=5,1e-7,
	(DT,Test,detD,CD,trace)=util.GradDifFin( [0,od] , 1 ,stri,X,N,T)
	# (DT,Test,detD,CD,trace)=util.GradDifFin( [3,3] , 1 ,stri,X,N,T)
	# Imax=DT.argmax()
	return( (max(T)-min(T))/max(DT) )
	# return( (max(T)-min(T))/DT[Imax],Imax )
#---------------------------------------------------------------------
def Thick2(X,T,N) :
	(DT,Test,detD,CD,trace)=util.GradDifFin( [3,3] , 1 ,1e-8,X,N,T) ; DT=np.array(DT)
	Imax=DT.argmax()
	return( (max(T)-min(T))/DT[Imax],Imax,DT )
#---------------------------------------------------------------------
def Setting(Pi,schem,Rafcrit,Tol,Tstep,X0,Ti,L):
	Set1='H2_Pi={0:1.3e}_{1}'.format(Pi,schem)
	Set2='_Raf{0:.0f}_Tol={1:.0e}_Tstep={2:.0e}_L{5:1.2e}_X0={3:.2f}_Ti={4:.2f}'.format(Rafcrit,Tol[0]*Tol[1],Tstep,X0,Ti,L)
	return(Set1,Set2,Set1 + Set2)
#---------------------------------------------------------------------
def ConsumptionSpeed(f) :

    Tb    = max(f.T      )
    Tu    = min(f.T      )
    rho_u = max(f.density)

    qx = f.heat_release_rate/f.cp

    Q0 = np.trapz(qx, f.grid)
    Sc = Q0/((Tb - Tu)*rho_u)

    return(Sc)
#---------------------------------------------------------------------
def StrainRates(f):
    grid=f.grid
    U   =f.velocity
    K = util.Deriv(grid,U)

    # fig,ax=plt.subplots()
    # ax.plot( grid,U,'k' )
    # ax.plot( grid,K,'r' )
    # plt.show()
    # plt.close(fig)

    IMaxS = abs(K).argmax()
    IMinU = f.velocity[:IMaxS].argmin()

    Ic = abs(K[:IMinU]).argmax()
    Kc = abs(K[Ic])

    Hr=f.heat_release_rate
    ImaxHr=Hr.argmax()
    cutHr=0.01*max(Hr)
    f_U=interp1d( Hr[:ImaxHr],U   [:ImaxHr] )
    f_K=interp1d( Hr[:ImaxHr],K   [:ImaxHr] )
    f_X=interp1d( Hr[:ImaxHr],grid[:ImaxHr] )

    return( K, Ic,Kc , f_X(cutHr),f_U(cutHr),f_K(cutHr) , IMaxS,IMinU,ImaxHr )
#---------------------------------------------------------------------