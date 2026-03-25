#!/usr/bin/env python3
#-*- coding: utf-8 -*-

#=================================================================================
def Compo(mdot_kk,y_caco3_feed,y_alooh_feed,y_h2o_caco3,y_h2o_alooh,y_co2_caco3) :
	y_h2o_feed=(y_h2o_caco3*y_caco3_feed+y_h2o_alooh*y_alooh_feed) # mass fraction of humidity in raw material
	y_co2_feed=(y_co2_caco3*y_caco3_feed                         ) # mass fraction of CO2      in raw material
	y_kk=1-(y_h2o_feed+y_co2_feed) # mass fraction of klinker in raw material
	mdot_rm=mdot_kk/y_kk           # mass flow rate of raw material

	mdot_h2o=y_h2o_feed*mdot_rm 
	mdot_co2=y_co2_feed*mdot_rm 

	mdot_h2o*=(1e3/(3600*24)) # [kg/s]
	mdot_co2*=(1e3/(3600*24)) # [kg/s]