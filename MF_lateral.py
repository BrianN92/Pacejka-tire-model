'''
Pacejka Tire Model
Calculate parameters and forces of lateral dynamics ONLY
B: Stiffness factor
C: Shape factor
D: Peak factor
E: Curvature factor
Sh: Horizontal shift
Sv: Vertical shift
'''
import numpy as np
import os
import warnings
import yaml
import math

class MF_lateral:
	def __init__(self, filename, data_path):
		self.param_tires = {}
		for ele in filename:
			self.param_tires[ele] = self.load_params(ele, data_path)
		self.epsilon   = 1e-6
		self.zeta = 1
		
	def online_params(self, Fz, p, gamma, kappa, alpha, Vcx, tire_select):
		self.params = self.param_tires[tire_select]
		self.Fz    = Fz
		self.p     = p
		self.gamma = gamma
		self.kappa = kappa
		self.alpha = alpha
		self.Vcx   = Vcx
		self.Vsx   = -self.kappa * np.abs(self.Vcx)
		self.Vsy   = np.tan(self.alpha) * np.abs(self.Vcx) 
		self.Vs = np.sqrt(np.power(self.Vsx,2) + np.power(self.Vsy,2))
		self.LMUY_star = self.params['LMUY'] / (1 + self.params['LMUY'] * self.Vs / self.params['LONGVL']) 	
		self.LMUY_prime = self.LMUY_star / (1 * self.LMUY_star) 
		self.Fz0_prime =  self.params['LFZO'] * self.params['FNOMIN'] 	
		# Unpack Parameters
		self.dfz = (self.Fz - self.Fz0_prime) / self.Fz0_prime	
		self.dpi = (self.p - self.params['NOMPRES']) / self.params['NOMPRES'] 

	def load_params(self, filename, data_path):
		data_ls = [data for data in os.listdir(data_path)]
		idx_ls  = [idx for idx, ltr in enumerate(data_ls) if filename in ltr]
		name    = data_ls[idx_ls[0]]
		with open(data_path + name, 'r') as stream:
			try:
				parsed_yaml = yaml.safe_load(stream)
			except yaml.YAMLError as exc:
				print(exc)
		return parsed_yaml

	def calculateFy0(self):
		Fz    = self.Fz
		# [SCALING_COEFFICIENTS]
		LCY   = self.params['LCY']   # Scale factor of Fy shape factor
		LEY   = self.params['LEY']   # Scale factor of Fy curvature factor
		LKY   = self.params['LKY']   # Scale factor of Fy cornering stiffness
		LHY   = self.params['LHY']   # Scale factor of Fy horizontal shift
		LVY   = self.params['LVY']   # Scale factor of Fy vertical shift
		LKYC  = self.params['LKYC']  # Scale factor of camber force stiffness
		
		# [LATERAL_COEFFICIENTS]
		PCY1  =  self.params['PCY1'] 	#Shape factor Cfy for lateral forces
		PDY1  =  self.params['PDY1'] 	#Lateral friction Muy
		PDY2  =  self.params['PDY2'] 	#Variation of friction Muy with load
		PDY3  =  self.params['PDY3'] 	#Variation of friction Muy with squared camber
		PEY1  =  self.params['PEY1'] 	#Lateral curvature Efy at Fznom
		PEY2  =  self.params['PEY2'] 	#Variation of curvature Efy with load
		PEY3  =  self.params['PEY3'] 	#Zero order camber dependency of curvature Efy
		PEY4  =  self.params['PEY4'] 	#Variation of curvature Efy with camber
		PEY5  =  self.params['PEY5'] 	#Variation of curvature Efy with camber squared
		PKY1  =  self.params['PKY1'] 	#Maximum value of stiffness Kfy./Fznom
		PKY2  =  self.params['PKY2'] 	#Load at which Kfy reaches maximum value
		PKY3  =  self.params['PKY3'] 	#Variation of Kfy./Fznom with camber
		PKY4  =  self.params['PKY4'] 	#Curvature of stiffness Kfy
		PKY5  =  self.params['PKY5'] 	#Peak stiffness variation with camber squared
		PKY6  =  self.params['PKY6'] 	#Fy camber stiffness factor
		PKY7  =  self.params['PKY7'] 	#Vertical load dependency of camber stiffness
		PHY1  =  self.params['PHY1'] 	#Horizontal shift Shy at Fznom
		PHY2  =  self.params['PHY2'] 	#Variation of shift Shy with load
		PVY1  =  self.params['PVY1'] 	#Vertical shift in Svy./Fz at Fznom
		PVY2  =  self.params['PVY2'] 	#Variation of shift Svy./Fz with load
		PVY3  =  self.params['PVY3'] 	#Variation of shift Svy./Fz with camber
		PVY4  =  self.params['PVY4'] 	#Variation of shift Svy./Fz with camber and load
		PPY1  =  self.params['PPY1'] 	#influence of inflation pressure on cornering stiffness
		PPY2  =  self.params['PPY2'] 	#influence of inflation pressure on dependency of nominal tyre load on cornering stiffness
		PPY3  =  self.params['PPY3'] 	#linear influence of inflation pressure on lateral peak friction
		PPY4  =  self.params['PPY4'] 	#quadratic influence of inflation pressure on lateral peak friction
		PPY5  =  self.params['PPY5'] 	#Influence of inflation pressure on camber stiffness
		

		Kya = PKY1 * self.Fz0_prime * (1 + PPY1 * self.dpi) * (1 - PKY3 * abs(self.gamma)) * np.sin(PKY4 * np.arctan((Fz / self.Fz0_prime) / 
				((PKY2+PKY5 * np.power(self.gamma,2)) * (1+PPY2 * self.dpi)))) * self.zeta * LKY 	
		SVyg = Fz * (PVY3 + PVY4 * self.dfz) * self.gamma * LKYC * self.LMUY_prime * self.zeta			
		# MF6.1 and 6.2 equatons
		Kyg0 = Fz * (PKY6 + PKY7 * self.dfz) * (1 + PPY5 * self.dpi) * LKYC 							

		signKya  = np.sign(Kya)
		
		SHy = (PHY1 + PHY2 * self.dfz) * LHY + ((Kyg0 * self.gamma - SVyg) / (Kya + self.epsilon * signKya)) * self.zeta + self.zeta - 1 # (4.E27) 
		SVy = Fz * (PVY1 + PVY2 * self.dfz) * LVY * self.LMUY_prime * self.zeta + SVyg  # (4.E29)

		# low speed mode
		Wvlow 		   = 0.5 * (1 + np.cos(np.pi * (self.Vcx / self.params['VXLOW'])))
		self.rdctSmth  = 1 - Wvlow
		if self.Vcx < self.params['VXLOW']:
			SVy *= self.rdctSmth 
			SHy *= self.rdctSmth 
		alphay 	   = self.alpha + SHy 		# (4.E20)
		Cy 	       = PCY1 * LCY				# (4.E21)
		muy        = (PDY1 + PDY2 * self.dfz) * (1 + PPY3 * self.dpi + PPY4 * np.power(self.dpi,2)) * (1 - PDY3 * np.power(self.gamma,2)) * self.LMUY_star # (4.E23) 
		Dy 	       = muy * Fz * self.zeta		# (4.E22)
		signAlphaY = np.sign(alphay)
		Ey = (PEY1 + PEY2 * self.dfz) * (1 + PEY5 * np.power(self.gamma,2) - (PEY3 + PEY4 * self.gamma) * signAlphaY) * LEY					# (<=1)(4.E24)

		# Limits check
		if Ey>1:
			warnings.warn("Warning CoeffChecks: Ey over limit (>1), Eqn(4.E14)")
			Ey = 1
		signDy   = np.sign(Dy)			# If [Dy = 0] then [sign(0) = 0]. This is done to avoid [Kya / 0 = NaN] in Eqn 4.E26
		By = Kya / (Cy * Dy + self.epsilon * signDy) 																	# (4.E26) [sign(Dy) term explained on page 177]
		Fy0 = Dy * np.sin(Cy * np.arctan(By * alphay - Ey * (By * alphay - np.arctan(By * alphay))))+ SVy   		# (4.E19)
		muy = 0 if Fz == 0 else muy
		params = {
			'By': By,
			'Cy': Cy,
			'Dy': Dy,
			'Ey': Ey,
			'SVy': SVy,
			'SHy': SHy,
			'muy': muy
		}
		return Fy0, params
		
	def calculateFy(self):
		Fz    = self.Fz
		kappa = self.kappa
		alpha_star = self.alpha 
		gamma_star = self.gamma 
		Fy0, params = self.calculateFy0()
		muy = params['muy']
		# [SCALING_COEFFICIENTS]
		LYKA  =  self.params['LYKA']    #  Scale factor of alpha influence on Fx
		LVYKA =  self.params['LVYKA']   #  Scale factor of kappa induced Fy
		
		# [LATERAL_COEFFICIENTS]
		RBY1 =   self.params['RBY1']  # Slope factor for combined Fy reduction
		RBY2 =   self.params['RBY2']  # Variation of slope Fy reduction with alpha
		RBY3 =   self.params['RBY3']  # Shift term for alpha in slope Fy reduction
		RBY4 =   self.params['RBY4']  # Influence of camber on stiffness of Fy combined
		RCY1 =   self.params['RCY1']  # Shape factor for combined Fy reduction
		REY1 =   self.params['REY1']  # Curvature factor of combined Fy
		REY2 =   self.params['REY2']  # Curvature factor of combined Fy with load
		RHY1 =   self.params['RHY1']  # Shift factor for combined Fy reduction
		RHY2 =   self.params['RHY2']  # Shift factor for combined Fy reduction with load
		RVY1 =   self.params['RVY1']  # Kappa induced side force Svyk./Muy.*Fz at Fznom
		RVY2 =   self.params['RVY2']  # Variation of Svyk./Muy.*Fz with load
		RVY3 =   self.params['RVY3']  # Variation of Svyk./Muy.*Fz with camber
		RVY4 =   self.params['RVY4']  # Variation of Svyk./Muy.*Fz with alpha
		RVY5 =   self.params['RVY5']  # Variation of Svyk./Muy.*Fz with kappa
		RVY6 =   self.params['RVY6']  # Variation of Svyk./Muy.*Fz with atan(kappa)
		
		DVyk = muy * Fz * (RVY1 + RVY2 * self.dfz + RVY3 * gamma_star) * np.cos(np.arctan(RVY4 * alpha_star)) * self.zeta  #  (4.E67)
		SVyk = DVyk * np.sin(RVY5 * np.arctan(RVY6 * kappa)) * LVYKA  											  #  (4.E66)
		SHyk = RHY1 + RHY2 * self.dfz  #  (4.E65)
		Eyk  = REY1 + REY2 * self.dfz  #  (<=1) (4.E64)
		Cyk = RCY1				# (4.E63)
		Byk = (RBY1 + RBY4 * np.power(gamma_star,2)) * np.cos(np.arctan(RBY2 * (alpha_star - RBY3))) * LYKA       # (> 0) (4.E62)
		kappas = kappa + SHyk   # (4.E61)
		
		if Eyk>1:
			warnings.warn("Warning CoeffChecks: Ey over limit (>1), Eqn(4.E14)")
			Eyk = 1

		Gyk0 = np.cos(Cyk * np.arctan(Byk * SHyk - Eyk * (Byk * SHyk - np.arctan(Byk * SHyk))))                   # (4.E60)
		Gyk = np.cos(Cyk *np.arctan(Byk * kappas - Eyk * (Byk * kappas - np.arctan(Byk * kappas)))) / Gyk0        # (> 0)(4.E59)
		# Low speed model
		Vx    = self.Vcx
		if Vx <= self.params['VXLOW']:
			SVyk *= self.rdctSmth
		Fy = Gyk * Fy0 + SVyk  # (4.E58)
		return Fy, Gyk
	


		
