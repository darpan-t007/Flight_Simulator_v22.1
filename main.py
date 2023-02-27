import numpy as np
import pandas as pd
import math
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

data = pd.read_excel('rayquaza data.xlsx', header=None, usecols="A", nrows=11)

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
data1 = pd.read_excel('rayquaza aero data.xlsx', sheet_name=0) #Aerodynamic data without airbrakes
data2 = pd.read_excel('rayquaza aero data.xlsx', sheet_name=1) #Aerodynamic data with airbrakes

#Importing data of motor and its thrust-time curve is imported from an excel file
motor_data = pd.read_excel('M3400 WT.xlsx')
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
