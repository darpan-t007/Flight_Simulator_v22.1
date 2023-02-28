import numpy as np

def Rk4Derivative(thrust, drag, sidey, sidez, mass, u, v, w, p, q, r, cphi, cpsi, ctheta, sphi, spsi, stheta, ymoment, pmoment, rmoment, Ixx, Iyy, Izz):
    # Eqns from PDF Titled EOM Derived
    udot = (thrust/mass) - (drag/mass) - ctheta*cpsi*9.81 + p*v - q*w
    vdot = (sidey/mass) - 9.81*(sphi*stheta*cpsi - cphi*spsi) - p*u + r*w
    wdot = (sidez/mass) - 9.81*(cphi*stheta*cpsi + sphi*spsi) + q*u - r*v
    rdot = (1/Ixx)*(rmoment + q*p*(Iyy-Izz))
    qdot = (1/Iyy)*(pmoment + r*p*(Izz-Izz))
    pdot = (1/Izz)*(ymoment + r*q*(Ixx-Iyy))
    phidot = r + (q*spsi + p*cphi)*(stheta/ctheta) 
    thetadot = q*cphi - p*sphi
    psidot = (q*sphi + p*cphi)*(1/ctheta)
    xedot = ctheta*cphi*u + (-cphi*spsi + spsi*stheta*cpsi)*v + (sphi*spsi + cphi*stheta*cpsi)*w
    yedot = ctheta*spsi*u + (cphi*cpsi + sphi*stheta*spsi)*v + (-sphi*cpsi + cphi*stheta*spsi)*w 
    zedot = -stheta*u + sphi*ctheta*v + cphi*ctheta*w
    return udot, vdot, wdot, rdot, qdot, pdot, phidot, thetadot, psidot, xedot, yedot, zedot