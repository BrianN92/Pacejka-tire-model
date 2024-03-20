"""
Reference textbook: 
@ Vehicle Dynamics and Control, 2nd Edition
1). Normal load on the tires is influenced by:
		a) fore-aft location of the c.g.
		b) longitudinal acceleration of the vehicle
		c) aerodynamic drag forces on the vehicle
		d) grade (inclination) of the road 

2). Aerodynamic drag force:
	The equivalent aerodynamic drag force on a vehicle
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import yaml
import time
sys.path.append('tire_data/config/')
import vehicle


class NormalForces:
	def __init__(self, Vcx, acc, gamma, vehicle_params):
		self.Vcx    = np.ones(1) * Vcx
		self.acc    = np.ones(1) * acc
		self.gamma  = np.ones(1) * gamma
		self.h_cg   = vehicle_params['hcog']
		self.h_aero = vehicle_params['haero']
		self.lf     = vehicle_params['lf']
		self.lr     = vehicle_params['lr']
		self.mass   = vehicle_params['mass']
	
	def calculateFaero(self, Vwind):
		# Theoretical values 
		# Refer to: Wong, J.Y., Theory of Ground Vehicles, Wiley-Interscience, ISBN 0-471-35461-9, Third Edition, 2001
		rho = 1.225  						 	 # Mass density of air [kg/m^3]
		Af  = 1.6 + 0.00056 * (self.mass - 765)  # Frontal area
		Cd  = 0.8 								 # Aerodynamic drag coefficient: 0.7 to 1.1 typical values for Formula Once car
		# Calculate aerodynamic drag force
		F_aero  = 1/2 * rho * Cd * Af * (self.Vcx + Vwind)**2		# [Eq: 4.2 Book: VDC]
		return F_aero

	def calculateRx(self, Vwind):
		# The rolling resistance is modeled as being roughly proportionalto the normal force on each set of tires
		f_roll   = 0.015					# rolling resistance coefficient f varies in the range 0.01 to 0.04. (Wong, 2001)
		Fzf, Fzr = self.calculateFz(Vwind)  
		R_total  = f_roll * (Fzf + Fzr)		# [Eq: 4.15 Book: VDC]
		return R_total

	def calculateFz(self, Vwind):
		h_aero = self.h_aero
		h_cg   = self.h_cg
		l_r    = self.lr
		l_f    = self.lf
		mass   = self.mass
		theta  = self.gamma
		acc    = self.acc
		g      = 9.81			# gravitational acceleration
		F_aero = self.calculateFaero(Vwind)
		# Normal forces on front tire and rear tire
		Fzf = (-F_aero * h_aero - mass * acc * h_cg - mass * g * h_cg * np.sin(theta) + mass * g * l_r * np.cos(theta)) / (l_f + l_r) # [Eq: 4.17 Book: VDC]
		Fzr = ( F_aero * h_aero + mass * acc * h_cg + mass * g * h_cg * np.sin(theta) + mass * g * l_f * np.cos(theta)) / (l_f + l_r) # [Eq: 4.18 Book: VDC]
		return Fzf, Fzr

if __name__ == '__main__':
	# Test case
	nPoints = 500                                     
	gamma	= np.linspace(-0.0873,0.0873,nPoints)                                                      
	Vx   	= np.ones(nPoints) * 37.1    
	Vwind   = Vx - 15
	acc     = np.linspace(0., 1., nPoints)           
	# Python solver
	start    = time.time()
	vehicle_params = vehicle.AV21()
	mf       = NormalForces(Vx, acc, gamma, vehicle_params)
	Fzf, Fzr = mf.calculateFz(Vwind)
	Rx       = mf.calculateRx(Vwind)
	end = time.time()
	print('Calculation time: %s sec' %(end - start))

	# Visualizing the results
	plt.figure()
	plt.plot(Fzf, label='Front tire')
	plt.plot(Fzr, label='Rear tire')
	plt.ylabel('$F_z$ [N]')
	plt.title("Normal Force $F_z$")
	plt.grid()
	plt.legend()

	plt.figure()
	plt.plot(Rx)
	plt.ylabel('$R_x$ [N]')
	plt.title("Rolling Resistance")
	plt.show()  