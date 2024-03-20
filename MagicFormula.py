"""
MF62-Tyre: Pacejka Tire Model
The equations coded in this script are published in the book:
Title: Tire and Vehicle Dynamics
Author: Hans Pacejka
Edition: 3, revised
Publisher: Elsevier, 2012
ISBN: 0080970176, 9780080970172
Link: https://www.elsevier.com/books/tire-and-vehicle-dynamics/pacejka/978-0-08-097016-5
"""

import numpy as np
import matplotlib.pyplot as plt
import yaml
import os
import warnings
import time


class MagicFormula:
	def __init__(self, filename, data_path, mode):
		self.param_tires = {}
		for ele in filename:
			self.param_tires[ele] = self.load_params(ele, data_path)
		self.mode      = mode
		tenth          = self.tenth_digit()
		hundredth	   = self.hundredth_digit()
		self.useStar   = True if tenth % 2 == 0 else False      # Only valid for the book implementation. TNO MF-Tyre does not use this.
		self.limtCheck = False if hundredth % 2 == 0 else True  # Check the sign of the coefficient of friction. The calculation of Fy is not affected by the sign of muy
		self.turnslip  = False 									# No turn slip and small camber angles
		self.epsilon   = 1e-6 									# [Eqn (4.E6a) Page 178 - Book]
		self.zeta = 1											# Create a vector with numbers between 0 and 1 to apply a reduction factor with smooth transitions.

	def online_params(self, Fz, p, gamma, kappa, alpha, Vcx, tire_select, phit=None, omega=None):
		self.params    = self.param_tires[tire_select]
		self.Fz        = Fz          # vertical load        (N)
		self.omega     = omega		# rotational speed     [rad/s]
		self.kappa     = kappa       # long. slip ratio     (-) (-1 = locked wheel)
		self.alpha     = alpha       # side slip angle      (rad)
		self.gamma     = gamma       # inclination angle    (rad)
		self.phit      = phit         # turn slip            (1/m)
		self.Vcx       = Vcx         # forward velocity     (m/s) 
		self.p         = p
		self.Fz0_prime =  self.params['LFZO'] * self.params['FNOMIN'] 	
		Wvlow 		   = 0.5 * (1 + np.cos(np.pi * (self.Vcx / self.params['VXLOW'])))
		self.rdctSmth  = 1 - Wvlow if self.limtCheck else 1
		self.dfz = (self.Fz - self.Fz0_prime) / self.Fz0_prime	
		self.dpi = (self.p - self.params['NOMPRES']) / self.params['NOMPRES'] 

	def tenth_digit(self):
		digits = []
		number = self.mode
		
		while number != 0:
			digits.append(number % 10)
			number = number // 10
		return digits[1]

	def hundredth_digit(self):
		number = self.mode
		digit = int(number // 100) % 10
		return digit

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

	def calculateBasics(self):
		# Unpack Parameters
		self.V0   = self.params['LONGVL']  									# Nominal speed
		self.pi0  = self.params['NOMPRES'] 			   						# Nominal tyre inflation pressure
		self.Fz0  = self.params['FNOMIN']  									# Nominal wheel load
		self.LFZO = self.params['LFZO']  									# Scale factor of nominal (rated) load
		self.LMUX = self.params['LMUX']  									# Scale factor of Fx peak friction coefficient
		self.LMUY = self.params['LMUY']  									# Scale factor of Fy peak friction coefficient
		self.LMUV = 0														# Scale factor with slip speed Vs decaying friction
		# Velocities in point S (slip point)
		Vsx = -self.kappa * np.abs(self.Vcx) 								# [Eqn (4.E5) Page 181 - Book]
		Vsy = np.tan(self.alpha) * np.abs(self.Vcx) 						# [Eqn (2.12) Page 67 - Book] and [(4.E3) Page 177 - Book]
		Vs = np.sqrt(np.power(Vsx,2) + np.power(Vsy,2)) 					# [Eqn (3.39) Page 102 - Book] -> Slip velocity of the slip point S
		# Velocities in point C (contact)
		self.Vcy = Vsy 														# Assumption from page 67 of the book, paragraph above Eqn (2.11)
		self.Vc = np.sqrt(np.power(self.Vcx,2) + np.power(self.Vcy,2)) 		# Velocity of the wheel contact centre C, Not described in the book 
																			# but is the same as [Eqn (3.39) Page 102 - Book]
		# Effect of having a tire with a different nominal load
		Fz0_prime =  self.params['LFZO'] * self.params['FNOMIN'] 				# [Eqn (4.E1) Page 177 - Book]																
		# Normalized change in vertical load
		dfz = (self.Fz - Fz0_prime) / Fz0_prime									# [Eqn (4.E2a) Page 177 - Book]
		# Normalized change in inflation pressure
		dpi = (self.p - self.params['NOMPRES']) / self.params['NOMPRES'] 		# [Eqn (4.E2b) Page 177 - Book]
		# Use of star (*) definition. Only valid for the book implementation. TNO MF-Tyre does not use this.
		if self.useStar:
			alpha_star = np.tan(self.alpha) * np.sign(self.Vcx) 			# [Eqn (4.E3) Page 177 - Book]
			gamma_star = np.sin(self.gamma) 								# [Eqn (4.E4) Page 177 - Book]
		else:
			alpha_star = self.alpha
			gamma_star = self.gamma
		# For the aligning torque at high slip angles
		self.signVc = np.sign(self.Vc)
		self.signVc = 1 if np.sign(self.Vc) == 0 else self.signVc
		epsilonv = self.epsilon
		Vc_prime = self.Vc + epsilonv * self.signVc 						# [Eqn (4.E6a) Page 178 - Book] [sign(Vc) term explained on page 177]
		alpha_prime = np.arccos(self.Vcx / Vc_prime)						# [Eqn (4.E6) Page 177 - Book]
		# Slippery surface with friction decaying with increasing (slip) speed
		LMUX_star = self.LMUX / (1 + self.LMUV * Vs / self.V0) 				# [Eqn (4.E7) Page 179 - Book]
		LMUY_star = self.LMUY / (1 + self.LMUV * Vs / self.V0) 				# [Eqn (4.E7) Page 179 - Book]
		# Digressive friction factor
		# On Page 179 of the book is suggested Amu = 10, but after
		# comparing the use of the scaling factors against TNO, Amu = 1
		# was giving perfect match
		Amu = 1
		LMUX_prime = Amu * LMUX_star / (1+(Amu-1) * LMUX_star) # [Eqn (4.E8) Page 179 - Book]
		LMUY_prime = Amu * LMUY_star / (1+(Amu-1) * LMUY_star) # [Eqn (4.E8) Page 179 - Book]
		# Pack output
		starVar = {
			'alpha_star': alpha_star,
			'gamma_star': gamma_star,
			'LMUX_star' : LMUX_star,
			'LMUY_star' : LMUY_star,
		}
		primeVar = {
			'Fz0_prime'  : Fz0_prime,
			'alpha_prime': alpha_prime,
			'LMUX_prime' : LMUX_prime,
			'LMUY_prime' : LMUY_prime,
		}
		incrVar = {
			'dfz': dfz,
			'dpi': dpi,
		}
		slipVar = {
			'Vsx': Vsx,
			'Vsy': Vsy,
			'Vs' : Vs
		}
		return starVar, primeVar, incrVar, slipVar

	def calculateFx0(self):
		# [SCALING_COEFFICIENTS]
		LCX  =  self.params['LCX']  #Scale factor of Fx shape factor
		LEX  =  self.params['LEX']  #Scale factor of Fx curvature factor
		LKX  =  self.params['LKX']  #Scale factor of Fx slip stiffness
		LHX  =  self.params['LHX']  #Scale factor of Fx horizontal shift
		LVX  =  self.params['LVX']  #Scale factor of Fx vertical shift
			
		# [LONGITUDINAL_COEFFICIENTS]
		PCX1  =   self.params['PCX1']  #Shape factor Cfx for longitudinal force
		PDX1  =   self.params['PDX1']  #Longitudinal friction Mux at Fznom
		PDX2  =   self.params['PDX2']  #Variation of friction Mux with load
		PDX3  =   self.params['PDX3']  #Variation of friction Mux with camber squared
		PEX1  =   self.params['PEX1']  #Longitudinal curvature Efx at Fznom
		PEX2  =   self.params['PEX2']  #Variation of curvature Efx with load
		PEX3  =   self.params['PEX3']  #Variation of curvature Efx with load squared
		PEX4  =   self.params['PEX4']  #Factor in curvature Efx while driving
		PKX1  =   self.params['PKX1']  #Longitudinal slip stiffness Kfx./Fz at Fznom
		PKX2  =   self.params['PKX2']  #Variation of slip stiffness Kfx./Fz with load
		PKX3  =   self.params['PKX3']  #Exponent in slip stiffness Kfx./Fz with load
		PHX1  =   self.params['PHX1']  #Horizontal shift Shx at Fznom
		PHX2  =   self.params['PHX2']  #Variation of shift Shx with load
		PVX1  =   self.params['PVX1']  #Vertical shift Svx./Fz at Fznom
		PVX2  =   self.params['PVX2']  #Variation of shift Svx./Fz with load
		PPX1  =   self.params['PPX1']  #linear influence of inflation pressure on longitudinal slip stiffness
		PPX2  =   self.params['PPX2']  #quadratic influence of inflation pressure on longitudinal slip stiffness
		PPX3  =   self.params['PPX3']  #linear influence of inflation pressure on peak longitudinal friction
		PPX4  =   self.params['PPX4']  #quadratic influence of inflation pressure on peak longitudinal friction

		# Unpack parameters
		epsilonX = self.epsilon
		starVar, primeVar, incrVar, slipVar = self.calculateBasics()
		
		dfz = incrVar['dfz']
		dpi = incrVar['dpi']

		LMUX_star  = starVar['LMUX_star']
		LMUX_prime = primeVar['LMUX_prime']

		Fz    = self.Fz
		kappa = self.kappa
		gamma = self.gamma
		Vx    = self.Vcx
		zeta1 = 1

		useLimitsCheck = self.limtCheck 	# The limit checks verify that the inputs are inside the stable range of the model

		Cx 	    = PCX1 * LCX 					# (> 0) (4.E11)
		mux     = (PDX1 + PDX2 * dfz) * (1 + PPX3 * dpi + PPX4 * np.power(dpi,2)) * (1 - PDX3 * np.power(gamma,2)) * LMUX_star
		if Fz == 0:
			mux = 0						# Zero Fz correction
		Dx 	    = mux * Fz * zeta1 				# (> 0) (4.E12)
		Kxk     = Fz * (PKX1 + PKX2 * dfz) * np.exp(PKX3 * dfz) * (1 + PPX1 * dpi + PPX2 * np.power(dpi,2)) * LKX  # (= BxCxDx = dFxo./dkx at kappax = 0) (= Cfk) (4.E15)
		signDx = np.sign(Dx)					# If [Dx = 0] then [sign(0) = 0]. This is done to avoid [Kxk / 0 = NaN] in Eqn 4.E16
		signDx = 1 if np.sign(Dx) == 0 else signDx
		Bx  = Kxk / (Cx * Dx + epsilonX * signDx) 						# (4.E16) [sign(Dx) term explained on page 177]
		SHx = (PHX1 + PHX2 * dfz) * LHX			  						# (4.E17)
		SVx = Fz * (PVX1 + PVX2 * dfz) * LVX * LMUX_prime * zeta1 		# (4.E18)

		# Low speed mode
		if Vx <= self.params['VXLOW']: # len(lowspeed_index) > 0 and lowspeed_index[0].size>0:
			SVx *= self.rdctSmth
			SHx *= self.rdctSmth
		
		kappax = kappa + SHx 																	   # (4.E10)
		Ex 	   = (PEX1 + PEX2 * dfz + PEX3 * np.power(dfz,2)) * (1 - PEX4 * np.sign(kappax)) * LEX # (<=1) (4.E14)

		# Limit check
		if useLimitsCheck and Ex > 1:
			warnings.warn("Warning CoeffChecks: Ex over limit (>1), Eqn(4.E14)")
			Ex = 1
		
		# Pure longitudinal force
		Fx0 = Dx * np.sin(Cx * np.arctan(Bx * kappax - Ex * (Bx * kappax - np.arctan(Bx * kappax)))) + SVx # (4.E9)
		params = {
			'Bx': Bx,
			'Cx': Cx,
			'Dx': Dx,
			'Ex': Ex,
			'SVx': SVx,
			'SHx': SHx,
			'mux': mux
		}
		return Fx0, params


	def calculateFy0(self):
		# Unpack Parameters
		epsilonk = self.epsilon
		epsilony = self.epsilon
		
		starVar, primeVar, incrVar, slipVar = self.calculateBasics()
		
		dfz = incrVar['dfz']
		dpi = incrVar['dpi']

		LMUY_star  = starVar['LMUY_star']
		alpha_star = starVar['alpha_star']
		gamma_star = starVar['gamma_star']
		
		LMUY_prime = primeVar['LMUY_prime']
		Fz0_prime  = primeVar['Fz0_prime']
		
		Fz    = self.Fz
		Vcx   = self.Vcx
		
		useLimitsCheck = self.limtCheck 	# The limit checks verify that the inputs are inside the stable range of the model
		useAlphaStar = self.useStar
		useTurnSlip = self.turnslip			

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
		
		if not useTurnSlip:
			zeta2 = 1
			zeta3 = 1
		Kya = PKY1 * Fz0_prime * (1 + PPY1 * dpi) * (1 - PKY3 * abs(gamma_star)) * np.sin(PKY4 * np.arctan((Fz / Fz0_prime) / 
				((PKY2+PKY5 * np.power(gamma_star,2)) * (1+PPY2 * dpi)))) * zeta3 * LKY 	# (4.E25) 
		SVyg = Fz * (PVY3 + PVY4 * dfz) * gamma_star * LKYC * LMUY_prime * zeta2 			# (4.E28)
		# MF6.1 and 6.2 equatons
		Kyg0 = Fz * (PKY6 + PKY7 * dfz) * (1 + PPY5 * dpi) * LKYC 							# (=dFyo./dgamma at alpha = gamma = 0) (= CFgamma) (4.E30)
		# First paragraph on page 178 of the book
		zeta0 = 1
		zeta4 = 1

		# [sign(Kya) term explained on page 177]
		signKya = np.sign(Kya) # This is done to avoid [num / 0 = NaN] in Eqn 4.E27
		signKya = 1 if np.sign(Kya) == 0 else signKya
		
		SHy = (PHY1 + PHY2 * dfz) * LHY + ((Kyg0 * gamma_star - SVyg) / (Kya + epsilonk * signKya)) * zeta0 + zeta4 - 1 # (4.E27) 
		SVy = Fz * (PVY1 + PVY2 * dfz) * LVY * LMUY_prime * zeta2 + SVyg  # (4.E29)

		# low speed mode
		Vx = self.Vcx
		
		if Vx <= self.params['VXLOW']:
			SVy *= self.rdctSmth
			SHy *= self.rdctSmth
		alphay 	   = alpha_star + SHy 		# (4.E20)
		Cy 	       = PCY1 * LCY				# (4.E21)
		muy        = (PDY1 + PDY2 * dfz) * (1 + PPY3 * dpi + PPY4 * np.power(dpi,2)) * (1 - PDY3 * np.power(gamma_star,2)) * LMUY_star # (4.E23) 
		Dy 	       = muy * Fz * zeta2		# (4.E22)
		
		signAlphaY = np.sign(alphay)
		signAlphaY = 1 if np.sign(alphay) == 0 else signAlphaY
		Ey = (PEY1 + PEY2 * dfz) * (1 + PEY5 * np.power(gamma_star,2) - (PEY3 + PEY4 * gamma_star) * signAlphaY) * LEY					# (<=1)(4.E24)

		# Limits check
		if useLimitsCheck and Ey>1:
			warnings.warn("Warning CoeffChecks: Ey over limit (>1), Eqn(4.E14)")
			Ey = 1

		signDy = np.sign(Dy)			
		signDy = 1 if signDy==0 else signDy   # If [Dy = 0] then [sign(0) = 0]. This is done to avoid [Kya / 0 = NaN] in Eqn 4.E26
		By = Kya / (Cy * Dy + epsilony * signDy) 																	# (4.E26) [sign(Dy) term explained on page 177]
		Fy0 = Dy * np.sin(Cy * np.arctan(By * alphay - Ey * (By * alphay - np.arctan(By * alphay))))+ SVy   		# (4.E19)
		# Zero Fz correction
		if Fz == 0:
			muy = 0
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

	def calculateFx(self):
		# Unpack Parameters
		useTurnSlip    = self.turnslip
		useLimitsCheck = self.limtCheck 
		kappa = self.kappa
		starVar, primeVar, incrVar, slipVar = self.calculateBasics()
		alpha_star = starVar['alpha_star'] 
		gamma_star = starVar['gamma_star'] 

		dfz = incrVar['dfz'] 
		
		# [SCALING_COEFFICIENTS]
		LXAL =  self.params['LXAL']  # Scale factor of alpha influence on Fx
		
		# [LONGITUDINAL_COEFFICIENTS]
		RBX1 = self.params['RBX1']   # Slope factor for combined slip Fx reduction
		RBX2 = self.params['RBX2']   # Variation of slope Fx reduction with kappa
		RBX3 = self.params['RBX3']   # Influence of camber on stiffness for Fx combined
		RCX1 = self.params['RCX1']   # Shape factor for combined slip Fx reduction
		REX1 = self.params['REX1']   # Curvature factor of combined Fx
		REX2 = self.params['REX2']   # Curvature factor of combined Fx with load
		RHX1 = self.params['RHX1']   # Shift factor for combined slip Fx reduction
		
		Cxa = RCX1  			 # (4.E55)
		Exa = REX1 + REX2 * dfz  # (<= 1) (4.E56)
		# Limits check
		if useLimitsCheck and Exa>1:
			warnings.warn("Warning CoeffChecks: Ey over limit (>1), Eqn(4.E14)")
			Exa = 1
		
		SHxa    = RHX1 			     # (4.E53)
		Bxa     = (RBX1 + RBX3 * np.power(gamma_star,2)) *np.cos(np.arctan(RBX2 * kappa)) * LXAL	  						  # (> 0) (4.E54)
		alphas  = alpha_star + SHxa   # (4.E53)
		Gxa0    = np.cos(Cxa * np.arctan(Bxa * SHxa - Exa * (Bxa * SHxa - np.arctan(Bxa * SHxa))))   						  # (4.E52)
		Gxa     = np.cos(Cxa * np.arctan(Bxa * alphas - Exa * (Bxa * alphas - np.arctan(Bxa * alphas)))) / Gxa0 			  # (> 0)(4.E51)
		Fx0,_   = self.calculateFx0()
		Fx      = Gxa * Fx0    
		return Fx, Gxa

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
		

	def calculateRe(self):
		# Unpack Parameters
		Vcx = self.Vcx
		Fz_unlimited = self.Fz
		kappa_unlimited = self.kappa
		
		# Rename the TIR file variables in the Pacejka style
		Fz0	= self.params['FNOMIN'] 		   	  # Nominal (rated) wheel load
		R0	= self.params['UNLOADED_RADIUS']      # Free tyre radius
		V0	= self.params['LONGVL']               # 'Nominal speed (LONGVL)
		Cz0	= self.params['VERTICAL_STIFFNESS']   #  Vertical stiffness
		
		PFZ1	= self.params['PFZ1']    # Pressure effect on vertical stiffness
		BREFF	= self.params['BREFF']   #  %Low load stiffness effective rolling radius
		DREFF	= self.params['DREFF']   #  %Peak value of effective rolling radius
		FREFF	= self.params['FREFF']   #  %High load stiffness effective rolling radius
		Q_RE0	= self.params['Q_RE0']   #  %Ratio of free tyre radius with nominal tyre radius
		Q_V1	= self.params['Q_V1']    # Tyre radius increase with speed
		
		starVar, primeVar, incrVar, slipVar = self.calculateBasics()
		
		dfz = incrVar['dfz']
		dpi = incrVar['dpi']

		# C code declaration
		omega = self.omega 				 # rotational speed (rad/s)
		Romega = R0 * (Q_RE0 + Q_V1 * np.power(((omega * R0) / V0),2)) # [Eqn (1) Page 2 - Paper] - Centrifugal growth of the free tyre radius
		Re = np.ones(Vcx.shape[0]) * (R0 * 0.965)
		# Nominal stiffness (pressure corrected)
		Cz = Cz0 * (1 + PFZ1 * dpi)      # [Eqn (5) Page 2 - Paper] - Vertical stiffness adapted for tyre inflation pressure
		# Check if omega is one of the inputs
		# If it is, use it to calculate Re, otherwise it can be approximated with a short iteration
		if self.omega.any(): 
			eps   = np.nextafter(0,1)
			omega = self.omega + eps # rotational speed (rad/s) [eps is added to avoid Romega = 0]
			Romega = R0  * (Q_RE0 + Q_V1 * np.power(((omega * R0) / V0),2))   # [Eqn (1) Page 2 - Paper] - Centrifugal growth of the free tyre radius
			# Eff. Roll. Radius
			Re = Romega - (Fz0 / Cz) * ( DREFF * np.arctan(BREFF * (Fz_unlimited / Fz0)) + FREFF * (Fz_unlimited / Fz0)) # [Eqn (7) Page 2 - Paper]
		else:
			# Omega is not specified and is going to be approximated
			# Initial guess of Re based on something slightly less than R0
			
			Re_old = Vcx * R0
			
			while (np.abs(Re_old - Re).any()>1e-9):
				Re_old = Re
				
				# Use the most up to date Re to calculate an omega
				# omega = Vcx ./ Re; % Old definition of Henning without kappa, not valid for brake and drive
				omega = np.real((kappa_unlimited * Vcx + Vcx) / Re)  # [Eqn (2.5) Page 65 - Book]
				
				# Then we calculate free-spinning radius
				Romega = R0 * (Q_RE0 + Q_V1 * np.power(((omega * R0) / V0),2))  # [Eqn (1) Page 2 - Paper] - Centrifugal growth of the free tyre radius
				
				# Effective Rolling Radius
				Re = Romega - (Fz0 / Cz) * (DREFF * np.arctan(BREFF * (Fz_unlimited / Fz0)) + FREFF * (Fz_unlimited / Fz0))   # [Eqn (7) Page 2 - Paper]
		return Re
