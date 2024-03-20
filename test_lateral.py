import numpy as np
import time 
import matplotlib.pyplot as plt
from MagicFormula import MagicFormula
from MF_lateral import MF_lateral
import argparse
import os 

# Initialize the argparse parser
parser = argparse.ArgumentParser()

# Add arguments
parser.add_argument("-m", "--model", default='lateral', type=str, help="tire model")

# Parse the arguments
args = parser.parse_args()
script_dir = os.path.dirname(os.path.abspath(__file__))
data_path  = os.path.join(script_dir, 'tire_data', 'yamls/')
# Test case
nPoints = 500
Fz      = np.ones(nPoints) * 3112                      
kappa	= np.linspace(-1., 1., nPoints)                                      
alpha	= np.linspace(-0.523599, 0.523599, nPoints)                                       
gamma	= np.ones(nPoints) * 0.0173                          
phit 	= np.ones(nPoints) * 0            
omega   = np.ones(nPoints) * 0         
Vx   	= np.linspace(-0.523599, 66.523599, nPoints)  
acc     = np.linspace(0., 1., nPoints)            
P       = np.ones(nPoints) * 200568.49 
Vwind   = 0   


###################################
# Python solver
###################################
tire_F = 'FRONT'
tire_R = 'REAR'
usemode = 111
# Python solver
if args.model == 'lateral':
    model = MF_lateral([tire_F, tire_R], data_path)
elif args.model == 'MF':
    model = MagicFormula(usemode)
Fyf    = np.zeros(nPoints)
Fyf0   = np.zeros(nPoints)
Fyr    = np.zeros(nPoints)
Fyr0   = np.zeros(nPoints)
for idx in range(nPoints):
    # start = time.time()
    model.online_params(Fz[idx], P[idx], gamma[idx], kappa[idx], alpha[idx], Vx[idx], tire_F)
    Fyf0[idx], params = model.calculateFy0()
    Fyf[idx], Gy      = model.calculateFy()
    model.online_params(Fz[idx], P[idx], gamma[idx], kappa[idx], alpha[idx], Vx[idx], tire_R)
    Fyr0[idx], params = model.calculateFy0()
    Fyr[idx], Gy      = model.calculateFy()
    # end = time.time()
    # print(end-start)
# Visualizing results
plt.figure()
plt.plot(Fyf0, label='Pure lateral force, [front]')
plt.plot(Fyf, label='Lateral force, [front]')
plt.legend()
plt.grid()

plt.figure()
plt.plot(Fyr0, label='Pure lateral force, [rear]')
plt.plot(Fyr, label='Lateral force, [rear]')
plt.legend()
plt.grid()
plt.show()

print(params)