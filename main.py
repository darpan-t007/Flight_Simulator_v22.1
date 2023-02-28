import numpy as np
import pandas as pd
import math
from density_n_temp import density_n_temp
from Predictor import predictor
from Thrust_Misc import thrust_n_other_things
from side_forces import side_forces
"""
- Main function to initialize, calculate, integrate everything and plot desired results
- Axis System is right hand with X upwards, Y out of the screen and Z going right
- Body axis b1 aligned with X(roll), b2(pitch) with Y and b3(yaw) with z
- Check PDF Titled EOM Derived for all equations [RONAN, MAKE A NEW PDF IN LATEX]
"""

#Initializing all values
#Note: All units are in SI units and the angles are in rads.

acceleration = np.zeros((500000, 1)) # Acceleration of the Rocket(in body direction aligned to roll)
prediction = np.zeros((500000, 1)) # stores predicted apogee at each time step
u = np.zeros((500000, 1)) # Body frame axis aligned in direction of thrust & roll axis
v = np.zeros((500000, 1)) # Body frame axis aligned to pitch axis
w = np.zeros((500000, 1)) # Body frame axis aligned to yaw axis
p = np.zeros((500000, 1)) # Yaw Rate
q = np.zeros((500000, 1)) # Pitch Rate
r = np.zeros((500000, 1)) # Roll Rate
psi = np.zeros((500000, 1)) # Yaw Euler Axis Angle
theta = np.zeros((500000, 1)) # Pitch Euler Axis Angle
phi = np.zeros((500000, 1)) # Roll Euler Axis Angle
Xe = np.zeros((500000, 1)) # X Distance in earth frame co ordinates
Ye = np.zeros((500000, 1)) # Y Distance in earth frame co ordinates
Ze = np.zeros((500000, 1)) # Z Distance in earth frame co ordinates
timer = np.zeros((500000, 1)) # Keeps track of time variable
velx = np.zeros((500000, 1)) # Keeps track of velocity in X(upwards) direction
dt = 0.001 #increment in time steps
counter1 = 1 #A counter for flagging and running the loop
temp = 300 #Ground Temperature for density calculations
CP = np.zeros((500000, 1)) #stores centre of pressure position
C_roll = np.zeros((500000, 1)) #stores roll moment coefficient
Cn_yaw = np.zeros((500000, 1)) #stores normal force coefficient for yaw
Cn_pitch = np.zeros((500000, 1)) #stores normal force coefficient for pitch
Cn_alpha = np.zeros((500000, 1)) #stores dCn/d alpha
Cd = np.zeros((500000, 1)) #stores drag coefficient
Mass = np.zeros((500000, 1)) #stores rocket mass
CG = np.zeros((500000, 1)) #Centre of Gravity
Ixx = np.zeros((500000, 1)) #Axial Moment of Inertia
Iyy = np.zeros((500000, 1)) #Transeverse Moment of Inertia
Stab_Cal = np.zeros((500000, 1)) #stability calibers of the rocket
flag = 0 #Flag for signifying state of the air brakes..0 implies airbrakes retracted, 1 implies airbrakes deployed
phase = 0 #Variable to keep track of flight phase..0 implies ascent phase 1 implies parachute phase
mphase = 1 #flag to keep track of motor burn out... value of 1 means motor is firing 0 means motor has been burnt out
state = np.zeros((500000, 1)) #flag to keep track of airbrake status...value of 2000 means airbrakes are deployed, 500 means airbrakes are retracted..these values and this variable is used to plot airbrake state
delay = 1 #Variable to account for delay in opening the air brakes in s
delaytracker = 0 #Variable to keep track of mechanical delay in the loop
turb_gen = np.zeros((4,4)) # helps in turbulence generation
turb_gen[0,0] = 1
turb_gen[0,3] = 1.5 # initial values for the pink noise generator for modelling turbulence in atmosphere
drift = np.zeros((500000,1)) # keeps track of the drifting of rocket from the launch point
vtrajectory = np.zeros((500000,3)) # variable for plotting rocket trajectory during ascent phase

"""
Assembly of rocket part and calculating their mass, cg and mi

To calculate the transverse moment of inertia, we need to know the position of the center of gravity (CG), which can only be determined in the main loop because the mass changes as the motor burns. During the motor burn, the transverse moment of inertia, mass, and distance of each component's CG from the nose tip are stored in the variable TMI. Once the motor burn is complete, the transverse moment of inertia can be calculated using the values stored in TMI.

"""

data = pd.read_excel('rocket_data.xlsx', header=None, usecols="A", nrows=11)

rmas = data.loc[0][0] #rocket mass without motor
rcg = data.loc[1][0] #rocket CG without motor
AMI = data.loc[2][0] #rocket axial moment of inertia without motor
TMI = data.loc[3][0] #rocket transverse moment of inertia without motor
mrad = data.loc[4][0] #rocket's max radius also used as refrence value for calculating aerodynamic forces and moments
rlen = data.loc[5][0] #rocket's total length
fin_details = data.loc[6:10, 0].values #5x1 matrix for rocket fin dimension [root chord,tip chord,fin semi span,sweep length,distance of fin from nose tip]
mrad = mrad/2
Aref = math.pi*mrad**2
mpos = rlen

#Asking User Input For Initial Pitch and Yaw Angle and atmospheric conditions
psi[0] = float(input('Enter Initial Yaw angle(degrees): '))
psi[0] = math.radians(psi[0])
theta[0] = float(input('Enter Initial Pitch angle(degrees): '))
theta[0] = math.radians(theta[0])

#Asking user for airbrakes
airbrake = input('Enter y to add airbrakes, Enter n to launch without airbrakes: ') #for simulation with or without airbrakes
if airbrake == 'y':
	desired = float(input('Enter Desired Apogee(m): '))

#User input for basic atmospheric conditions on launch site
awiny = -float(input('Enter mean wind speed in north direction: '))
awinz = float(input('Enter mean wind speed in east direction: '))
tur_inten = float(input('Enter turbulent intensity in perecentage in the atmosphere: '))

#Importing rocket's aerodynamic data
data1 = pd.read_excel('aero_data.xlsx', sheet_name=0) #Aerodynamic data without airbrakes
data2 = pd.read_excel('aero_data.xlsx', sheet_name=1) #Aerodynamic data with airbrakes

#Importing data of motor and its thrust-time curve is imported from an excel file
motor_data = pd.read_excel('motor_data.xlsx')
length = motor_data.iloc[0, 1]  # motor length
odia = motor_data.iloc[1, 1]  # motor outer diameter
wetmass = motor_data.iloc[2, 1]  # motor wet mass 
drymass = motor_data.iloc[3, 1]  # motor dry mass 
nozmass = motor_data.iloc[4, 1]  # nozzle mass
nozlen = motor_data.iloc[5, 1]  # nozzle length

thrust = np.zeros((16, 1))
time = np.zeros((16, 1))

for n in range(8, 23):
    thrust[n-8, 0] = motor_data.iloc[n, 1]
    time[n-8, 0] = motor_data.iloc[n, 0]

#Main Loop
while Xe[counter1] >= 0 or timer < time[0]:
	sphi = np.sin(phi[counter1])
	stheta = np.sin(theta[counter1])
	spsi = np.sin(psi[counter1])
	cphi = np.cos(phi[counter1])
	ctheta = np.cos(theta[counter1])
	cpsi = np.cos(psi[counter1])
	
	# mach no 
	mydensity, mytemp = density_n_temp(Xe[counter1], temp) # Calculating density and temperature at current altitude
	v_sound = np.sqrt(1.4 * 287 * mytemp)  # speed of sound at the current altitude
	mach_no = u[counter1] / v_sound  # mach number of rocket

	# calculating mass, centre of gravity, axial moment of inertia and transverse moment of inertia for the current time
	mmass, mcg, mIx, mIy, mythrust = thrust_n_other_things(timer[counter1], length, odia, wetmass, drymass, nozmass, nozlen, time, thrust)  # function call to calculate current thrust, instantaneous motor mass, cg and mi
	# mmas (= motor mass) is in kg, mcg (= position of motor cg wrt motor head) in cm
	# mIx and mIy are mi of motor in kgcm2
	# mythrust is current hrust in Newtons

	Mass[counter1] = rmas + mmass  # rocket's total mass
	CG[counter1] = ((rmas * rcg) + (mmass * (mpos - (length + nozlen) + mcg))) / Mass[counter1]  # net centre of gravity
	Ixx[counter1] = (AMI) + (mIx)  # net axial moment of inertia
	Iyy[counter1] = TMI + (rmas * (CG[counter1] - (rcg))**2) + (mIy) + (mmass * (CG[counter1] - (mpos - length + mcg))**2)  # net transeverse moment of inertia..rocket is assumed to be symmetric about xy and xz plane mi about y and z are same
"""
# condition for motor burn out
	if timer[counter1] > time[-1, 1]:
	    mphase = 0
	    
	if phase == 0:  # condition for rocket to be in ascent phase
	    # function call to interpolate aerodynamic data(cd,cn,cp etc) based on motor and airbrakes status(variables phase1 and flag)
	    CP[counter1], Cn_pitch[counter1], Cn_yaw[counter1], Cn_alpha[counter1], Cd[counter1] = aerodynamic_data(data1, data2, mach_no, psi[counter1], theta[counter1], mphase, flag)

	momentarm = CP[counter1] - CG[counter1]
	Stab_Cal[counter1] = momentarm / (2 * mrad)

	if timer[counter1] > time[-1] and phase == 0 and airbrake == 'y' and u[counter1] < 180:
	    # condition to ensure motor has burnt out, rocket is in ascent phase and the speed is below the speed allowed by structural limits (needs to be found from structural sims)
	    
	    if flag == 0:
	    	# prediction function call to predict apogee based on current altitude, vertical speed
	    	prediction[counter1] = predictor(data1, acceleration[counter1], u[counter1], Xe[counter1], Mass[counter1], Aref, mytemp)
	    
	    elif flag == 1:
	    	prediction[counter1] = predictor(data2, acceleration[counter1], u[counter1], Xe[counter1], Mass[counter1], Aref, mytemp)
		

	if prediction[counter1] > desired:  # Desired Apogee
	    if timer[counter1] - delaytracker >= delay: # check for mechanical delay
	    	flag = 1
	    	delaytracker = timer[counter1] 
		 
		
	elif prediction[counter1] < desired:  # Desired Apogee
	    if timer[counter1] - delaytracker >= delay:  # check for mechanical delay
	    	flag = 0
	    	delaytracker = timer[counter1]
		
		
	else:
		prediction[counter1] = 0


	if flag == 1:
	    state[counter1] = 2000  # Registering state of Air Brakes as open
	else:
	    state[counter1] = 500  # Registering state of Air Brakes as closed


	# Powered Flight & Coasting Phase
	if timer[counter1] < 10 or u[counter1] >= 0:
	    # generation of side wind gusts and turbulence
	    sidey, turb_gen = turbulence_generator(counter1, turb_gen, awiny, tur_inten, dt)
	    sidez, turb_gen = turbulence_generator(counter1, turb_gen, awinz, tur_inten, dt)
	    
	    # calculation of disturbing side forces and moments
	    Fy, Mz = side_forces(sidey, CG[counter1], mydensity, mrad, rlen, fin_details)
	    Fz, My = side_forces(sidez, CG[counter1], mydensity, mrad, rlen, fin_details)
	    
	    # My and Mz are disturbing moment about pitch and yaw respectively
	    Fy = Fy * np.sign(sidey) * np.cos(psi[counter1])
	    Fz = Fz * np.sign(sidez) * np.cos(theta[counter1])
	    My = My * np.sign(sidez) * np.cos(theta[counter1])
	    Mz = -Mz * np.sign(sidey) * np.cos(psi[counter1])
	    
	    if Xe[counter1] > 5.18:  # condition for launch rod clearance
	    	CNY = -(Cn_yaw[counter1] * np.sign(psi[counter1])) - (Cd[counter1] * np.sin(psi[counter1])) - (Cn_alpha[counter1] * np.arctan((q[counter1] * momentarm) / u[counter1]))  # net Cn yaw
	    	CNP = -(Cn_pitch[counter1] * np.sign(theta[counter1])) - (Cd[counter1] * np.sin(theta[counter1])) - (Cn_alpha[counter1] * np.arctan((p[counter1] * momentarm) / u[counter1]))  # net Cn pitch
	    	C_roll[counter1] = -2 * np.pi * np.arctan((1.5 * mrad * r[counter1]) / u[counter1]) * np.sign(r[counter1])  # net roll moment coefficient
		
	    else:
	    	CNY = 0
	    	CNP = 0
	    	C_roll[counter1] = 0
	    	My = 0
	    	Mz = 0
	    	Fy = 0
	    	Fz = 0
		
	drag = 0.5 * Cd[counter1] * u[counter1]**2 * mydensity * Aref # Calculating drag force
	Yaw = (0.5 * CNY * mydensity * u[counter1]**2 * Aref * momentarm) + Mz # Calculating net yaw moment
	Pitch = (0.5 * CNP * mydensity * u[counter1]**2 * Aref * momentarm) + My # Calculating net pitch moment
	Roll = C_roll[counter1] * mydensity * u[counter1]**2 * (fin_details[2][0] * (fin_details[0][0] + fin_details[1][0])) * (1.5 * mrad) # calculating net rolling moment

	# Store x,y,z values for trajectory for ascent
	vtrajectory[counter1][0] = Xe[counter1]
	vtrajectory[counter1][1] = Ye[counter1]
	vtrajectory[counter1][2] = Ze[counter1]

	# Parachute Phase
	else:
	    pdia = 0.95 # Parachute diameter in reefed condition
	    if Xe[counter1] < 460: # Condition for dis reef
		pdia = 3 # Un-reefed parachute diameter

	drag = (0.5 * 1.2 * math.pi * (pdia)**2 / 4) * mydensity * u[counter1]**2 * np.sign(u[counter1]) # 1.2 is parachute drag coefficient found from CFD
	Fy = (0.5 * 0.1 * mydensity * sidey**2 * math.pi * (pdia)**2 / 4 * np.sign(sidey)) - (0.5 * 0.1 * mydensity * v[counter1]**2 * math.pi * (pdia)**2 / 4 * np.sign(v[counter1])) # 0.1 is side force coefficient value from CFD
	Fz = (0.5 * 0.1 * mydensity * sidez**2 * math.pi * (pdia)**2 / 4 * np.sign(sidez)) - (0.5 * 0.1 * mydensity * w[counter1]**2 * math.pi * (pdia)**2 / 4 * np.sign(w[counter1]))
	Yaw = 0
	Pitch = 0
	Roll = 0
	phi[counter1] = 0
	psi[counter1] = 0
	theta[counter1] = 0
	phase = 1
	flag = 0
"""
#Add MATPLOTLIB here
"""
This is just to test the code
""" 
print(rmas)
print(rcg)
print(AMI)
print(TMI)
print(mrad)
print(rlen)
print(fin_details)
print(thrust)
print(time)
