import numpy as np
from scipy.interpolate import interp1d

def aerodynamic_data(data1, data2, mach, pitch, yaw, phase, flag):
    m = np.arange(0, 2, 0.01)
    cd_poweron = np.zeros((200, 3, 2))
    cd_poweroff = np.zeros((200, 3, 2))
    cn = np.zeros((200, 3))
    cn_alpha = np.zeros((200, 1))
    cp = np.zeros((200, 1))
    angle = np.array([0, 2, 4])
    mach = np.sqrt(mach ** 2)
    pitch = np.sqrt(np.rad2deg(pitch) ** 2)
    yaw = np.sqrt(np.rad2deg(yaw) ** 2)

    # without airbrakes
    cd_poweron[:, 0, 0] = data1[:200, 4]
    cd_poweron[:, 1, 0] = data1[2501:2700, 4]
    cd_poweron[:, 2, 0] = data1[5001:5200, 4]
    cd_poweroff[:, 0, 0] = data1[:200, 3]
    cd_poweroff[:, 1, 0] = data1[2501:2700, 3]
    cd_poweroff[:, 2, 0] = data1[5001:5200, 3]

    cn[:, 0] = data1[:200, 8]
    cn[:, 1] = data1[2500:2700, 8]
    cn[:, 2] = data1[5000:5200, 8]

    cn_alpha[:, 0] = data1[:200, 11]

    cp[:, 0] = data1[:200, 12]

    # with airbrakes
    cd_poweron[:, 0, 1] = data2[:200, 4]
    cd_poweron[:, 1, 1] = data2[2501:2700, 4]
    cd_poweron[:, 2, 1] = data2[5001:5200, 4]
    cd_poweroff[:, 0, 1] = data2[:200, 3]
    cd_poweroff[:, 1, 1] = data2[2501:2700, 3]
    cd_poweroff[:, 2, 1] = data2[5001:5200, 3]

    # interpolating data
    if phase == 1:
        ccd0 = interp1d(m, cd_poweron[:,0,0])(mach)
        cd2 = interp1d(m, cd_poweron[:,1,0])(mach)
        cd4 = interp1d(m, cd_poweron[:,2,0])(mach)
        p1 = np.polyfit(angle, [cd0, cd2, cd4], 2)
        CD = (p1[0]*(pitch**2+yaw**2)) + (p1[1]*np.sqrt(pitch**2+yaw**2)) + p1[2]
        cn0 = 0
        cn2 = interp1d(m, cn[:,1])(mach)
        cn4 = interp1d(m, cn[:,2])(mach)
        p2 = np.polyfit(angle, [cn0, cn2, cn4], 2)
        CNP = (p2[0]*pitch**2) + (p2[1]*pitch) + p1[2]
        CNY = (p2[0]*yaw**2) + (p2[1]*yaw) + p1[2]
        CP = interp1d(m, cp)(mach)

    else:
        if flag == 0:
            cd0 = interp1d(m, cd_poweroff[:, 1, 1])(mach)
            cd2 = interp1d(m, cd_poweroff[:, 2, 1])(mach)
            cd4 = interp1d(m, cd_poweroff[:, 3, 1])(mach)
            p1 = np.polyfit(angle, [cd0, cd2, cd4], 2)
            CD = (p1[0] * (pitch ** 2 + yaw ** 2)) + (p1[1] * np.sqrt(pitch ** 2 + yaw ** 2)) + p1[2]
            cn0 = 0
            cn2 = interp1d(m, cn[:, 2])(mach)
            cn4 = interp1d(m, cn[:, 3])(mach)
            p2 = np.polyfit(angle, [cn0, cn2, cn4], 2)
            CNP = (p2[0] * pitch ** 2) + (p2[1] * pitch) + p2[2]
            CNY = (p2[0] * yaw ** 2) + (p2[1] * yaw) + p2[2]
            CP = interp1d(m, cp.T[0])(mach)
        
        else:
            cd0 = interp1d(m, cd_poweroff[:, 0, 1], kind='linear')(mach)
            cd2 = interp1d(m, cd_poweroff[:, 1, 1], kind='linear')(mach)
            cd4 = interp1d(m, cd_poweroff[:, 2, 1], kind='linear')(mach)
            p1 = np.polyfit(angle, [cd0, cd2, cd4], 2)
            CD = (p1[0]*(pitch**2+yaw**2)) + (p1[1]*np.sqrt(pitch**2+yaw**2)) + p1[2]
            cn0 = 0
            cn2 = interp1d(m, cn[:, 1], kind='linear')(mach)
            cn4 = interp1d(m, cn[:, 2], kind='linear')(mach)
            p2 = np.polyfit(angle, [cn0, cn2, cn4], 2)
            CNP = (p2[0]*(pitch**2)) + (p2[1]*pitch) + p1[2]
            CNY = (p2[0]*(yaw**2)) + (p2[1]*yaw) + p1[2]
            CP = interp1d(m, cp, kind='linear')(mach)

    CN_alpha = np.interp(mach, m, cn_alpha)
    CP = CP/39.37 #Converting inches to meteres. Current datasheet is taken from RASAero II.

    return CN_alpha, CP, CD, CNP, CNY
