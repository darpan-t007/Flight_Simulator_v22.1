import numpy as np

def turbulence_generator(n, a, U, I, dt):
    """
    This function generates random gusts of wind whose values are centered around the mean to
    generate turbulence about the mean velocity using a pink noise generator.
    """
    
    SD = U * I * 0.01
    
    if n == 1:
        a[0,1] = np.random.randn()
    elif n == 2:
        a[1,0] = (2 - 5/6) * (a[0,0] / 2)
        a[1,1] = np.random.randn() - (a[1,0] * a[0,1])
        a[1,2] = a[0,2] + dt
    elif n == 3:
        a[2,0] = (3 - 5/6) * (a[1,0] / 3)
        a[2,1] = np.random.randn() - (a[2,0] * a[1,1]) - (a[1,0] * a[0,1])
        a[2,2] = a[1,2] + dt
        avg = (a[1,1] + a[0,1]) / 2
        a[0,3] = np.sqrt(((a[2,1] - avg)**2 + (a[2,1] - avg)**2) / 2)
    else:
        a[3,0] = (n - 5/6) * (a[2,0] / n)
        a[3,1] = np.random.randn() - (a[3,0] * a[2,1]) - (a[2,0] * a[1,1]) - (a[1,0] * a[0,1])
        a[3,2] = a[2,2] + dt
        avg = (a[2,1] + a[1,1] + a[0,1]) / 2
        a[0,3] = np.sqrt(((a[3,1] - avg)**2 + (a[3,1] - avg)**2 + (a[3,1] - avg)**2) / 3)
    
    u = U + ((a[3,1] / a[0,3]) * SD)
    
    if np.sqrt((u - U)**2) > SD:
        u = U + (np.sign(u) * SD)
    
    return u, a
