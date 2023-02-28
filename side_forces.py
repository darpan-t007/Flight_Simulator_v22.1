import numpy as np

def side_forces(vel, cm, den, rad, length, fin):
    # Magnitude of forces are calculated in this function. Direction is taken care of inside the main function.
    vel = np.sqrt(vel*vel)
    # Cd of cylinder at different Re
    cdcyl = np.array([[-5, 5], [10, 3.33], [100, 1.334], [1000, 0.738], [10000, 0.766], [100000, 0.56], [1000000, 0.55], [10000000, 0.54], [100000000, 0.54]])
    # Cd of flat plate perpendicular to flow
    cdplt = 1.28
    ryno = (1.225*vel*2*rad)/(1.73*10**-5)
    cn = np.interp(ryno, cdcyl[:,0], cdcyl[:,1])
    F = 0.5*den*vel*vel*((cn*rad*2*length)*(cdplt*(fin[0]+fin[1])*fin[2]))
    # Moment due to fins is only considered. Fuselage is assumed to be a cylindrical body so it does not produce any moments.
    mom = 0.5*den*vel*vel*cdplt*(fin[0]+fin[1])*fin[2]*(fin[4]-cm)
    return F, mom