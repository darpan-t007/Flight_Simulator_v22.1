def density_n_temp(altitude, t1):
    a = -0.0065  # Constant
    # Eqn for Temperature
    t = t1 + (a * altitude)
    # Density at given Height
    d = 1.23 * ((t / t1) ** ((-9.81 / (a * 287)) - 1))
    return d, t
