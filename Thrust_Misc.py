import numpy as np
def thrust_n_other_things(timer, length, odia, wetmass, drymass, nozmass, nozlen, time, thrust):
    # Ix and Iy are axial and transverse mi respectively
    if timer > time[-1]:  # motor has burnt out
        currentthrust = 0
        mass = drymass + nozmass  # total mass
        Ix1 = 0  # mi of motor grain
        Ix2 = drymass * (odia / 2) ** 2  # mi of casing assuming it a thin cylinder
        Ix3 = 0.5 * nozmass * 0.25 * (odia ** 2)  # mi of nozzle assuming it a thick cylinder
        cg = (((mass - nozmass) * (length / 2)) + (nozmass * (length + (nozlen / 2)))) / mass  # cg position wrt to motor head
        Iy1 = 0  # mi of grain
        Iy2 = (0.5 * drymass * (odia / 2) ** 2) + ((1 / 12) * drymass * length ** 2) + (drymass * (cg - (length / 2)) ** 2)  # mi of casing
        Iy3 = ((1 / 12) * nozmass * (((3 / 4) * (odia ** 2)) + nozlen ** 2)) + (nozmass * (cg - (length + (nozlen / 2))) ** 2)  # mi of nozzle
    else:
        currentthrust = np.interp(timer, time, thrust, 'spline')  # thrust for the given time
        burnrate = (wetmass - drymass) / time[-1]
        dm = burnrate * timer
        mass = (wetmass - dm) + nozmass  # total mass at the given time
        Ix1 = 0.5 * (mass - drymass - nozmass) * 0.25 * (odia ** 2)  # mi of motor grain modelled as a thick cylinder
        Ix2 = drymass * (odia / 2) ** 2  # mi of casing modelled as a thin cylinder
        Ix3 = 0.5 * nozmass * 0.25 * (odia ** 2)  # mi nozzle modelled as a thick cylinder
        cg = (((drymass) * (length / 2)) + ((mass - drymass) * ((length) / 2)) + (nozmass * (length + (nozlen / 2)))) / mass  # cg position wrt to motor head
        Iy1 = ((1 / 12) * (mass - drymass) * (((3 / 4) * (odia ** 2)) + (length) ** 2)) + ((mass - drymass) * (cg - ((length) / 2)) ** 2)
        Iy2 = (0.5 * drymass * (odia / 2) ** 2) + ((1 / 12) * drymass * length ** 2) + (drymass * (cg - (length / 2)) ** 2)
        Iy3=((1/12)*nozmass*(((3/4)*(odia**2))+nozlen**2))+(nozmass*(cg-(length+(nozlen/2)))**2)
        Ix=Ix1+Ix2+Ix3
        Iy=Iy1+Iy2+Iy3
        return mass, cg, Ix, Iy, currentthrust
