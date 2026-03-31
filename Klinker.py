#!/usr/bin/env python3
#-*- coding: utf-8 -*-

#======> Heat Capacity
cp={ # Kcal/K Kg
	'AlOOH' : 0.26 ,
	'Al2O3' : 0.28 ,
	'CACO3' : 0.26 ,
	'CAC'   : 0.26 ,
}
#======> Composition
Compo_feed={'AlOOH' : 0.512 ,'CaCO3' : 0.488 } # Compo raw mix 18-23 May
Compo_AlOOH={'H2O' : 0.112 , 'Fe2O3' : 0.2323 , 'SiO2' : 0.0427 , 'TiO2' : 0.0283 , 'AL2O3' : 0.5668 , 'Vap' : 0.005 } # Compo AlOOH
Compo_CACO3={'CO2' : 0.426 , 'CaO' : 0.539 , 'Vap' : 0.35e-2 } # Compo CaCO3

#======> Enthalpy
DH_dhyd=179 # [kcal/kg(AlOOH)] Enthalpy of dehydatation of AlOOH at 500 °C
DH_decb=396 # [kcal/kg(CaCO3)] Enthalpy of decarbonation of CaCO3 at 900 °C
DH_rcac=90  # [kcal/kg(CAC)]   Enthalpy of reaction of CAC 1300 °C
DH_fcac=145 # [kcal/kg(CAC)]   Enthalpy of fusion of CAC 1300 °C

#=================================================================================
def Compo2( mdot_kk,Compo_feed,Compo_AlOOH,Compo_CACO3 ) : 
	return( Compo( mdot_kk,Compo_feed['CaCO3'],Compo_feed['AlOOH'],Compo_CACO3['Vap'],Compo_AlOOH['Vap'],Compo_AlOOH['H2O'],Compo_CACO3['CO2'] ) )
def Compo(mdot_kk,y_caco3_feed,y_alooh_feed,y_vap_caco3,y_vap_alooh,y_h2o_alooh,y_co2_caco3) :
	y_h2o_feed=(y_vap_caco3*y_caco3_feed+(y_h2o_alooh+y_vap_alooh)*y_alooh_feed) # mass fraction of humidity in raw material
	y_co2_feed=(y_co2_caco3*y_caco3_feed                                       ) # mass fraction of CO2      in raw material
	y_kk=1-(y_h2o_feed+y_co2_feed) # mass fraction of klinker in raw material
	mdot_rm=mdot_kk/y_kk           # mass flow rate of raw material

	mdot_h2o=y_h2o_feed*mdot_rm 
	mdot_co2=y_co2_feed*mdot_rm 

	mdot_h2o*=(1e3/(3600*24)) # [kg/s]
	mdot_co2*=(1e3/(3600*24)) # [kg/s]

	return(mdot_h2o,mdot_co2,y_kk)
#=================================================================================
def Conv_thh(pow) : return(pow*4184e3/3600) # [th/h] -> [W] Thermal power