import numpy as np


def bilinear_interpolation(X, Y, f, x, y):
    """
    Interpolation methods for estimating surface current values in between grid points. Edge cases are outlined in the
    top of the function and may need to be refactored. NaNs is returned for cases where values can not be interpolated.

    Arguments
    ---------
    X, Y: Coordinate mesh grid
    f: Grid of velocity values that can be accessed as f(j,i) Remeber row, column
    x, y: coordinates to compute interpolation to f(y,x)

    Returns
    ---------
    interp_value: interpolated value of f(y,x)

    """
    # Grid index shape
    M = np.shape(X[:, 0])[0]
    N = np.shape(X[0, :])[0]

    dx, dy = X[0, 1] - X[0, 0], Y[1, 0] - Y[0, 0]
    x_start, y_start = X[0, 0], Y[0, 0]

    # Find the index of each value
    i1, i2 = int((x - x_start) / dx), int((x - x_start) / dx) + 1
    j1, j2 = int((y - y_start) / dy), int((y - y_start) / dy) + 1

    # Boundary Conditions when interpolating near the edge.
    # 1. Eastern boundary
    if (i1 - N) > 1:
        return np.nan
    if i1 >= N - 1 and j1 <= N - 1 and j1 >= 0: # If on the Eastern edge of the boundary
        return f[j1, N - 1]
    if i1 >= N - 1 and j1 <= 0:
        return f[0, N - 1]
    if i1 >= N - 1 and j1 >= N - 1:
        return f[N - 1, N - 1]

    # 2. Western boundary
    if i1 <= 0 and j1 <= N - 1 and j1 >= 0:
        return f[j1, 0]
    if i1 <= 0 and j1 <= 0:
        return f[0, 0]
    if i1 <= 0 and j1 >= N - 1:
        return f[N - 1, 0]

    # 3. Northern boundary
    if j1 >= M - 1 and i1 <= M - 1 and i1 >= 0:
        return f[M - 1, i1]
    if j1 >= N - 1 and i1 <= 0:
        return f[M - 1, 0]

    # 3. Bottom boundary
    if j1 <= 0 and i1 <= M - 1 and i1 >= 0:
        return f[0, i1]
    if j1 <= 0 and i1 >= M - 1:
        return f[M - 1, 0]

    x1, x2 = X[j1, i1], X[j2, i2]
    y1, y2 = Y[j1, i1], Y[j2, i2]

    interp_value = (1 / (x2 - x1) * 1 / (y2 - y1) *
                      (f[j1, i1] * (x2 - x) * (y2 - y) + f[j1, i2] * (x - x1) * (y2 - y)
                       + f[j2, i1] * (x2 - x) * (y - y1) + f[j2, i2] * (x - x1) * (y - y1)))

    return interp_value

def rk4(X, Y, x, y, f, h, dim):
    """
    Solves for the next position of a particle in after time, h, in either the x
    or y using a runge-kutta 4th order scheme.

    TODO: Update function to get next half timestep if goes into next hour

    Arguments
    ---------
    X, Y: mesh grid.
    x, y: coordinates where to begin the evolution.
    f: the current vector f that will be evolved.
    h: the time step
    dim: 0 for x and 1 for y.

    Returns
    ---------
    interp_value: interpolated value of f(y,x)
    """

    k1 = h * bilinear_interpolation(X, Y, f, x, y)
    k2 = h * bilinear_interpolation(X, Y, f, x + 0.5 * h, y + 0.5 * k1)
    k3 = h * bilinear_interpolation(X, Y, f, x + 0.5 * h, y + 0.5 * k2)
    k4 = h * bilinear_interpolation(X, Y, f, x + h, y + k3)
    try:
        if dim == 0:
            return x + 1. / 6 * k1 + 1. / 3 * k2 + 1. / 3 * k3 + 1. / 6 * k4
        elif dim == 1:
            return y + 1. / 6 * k1 + 1. / 3 * k2 + 1. / 3 * k3 + 1. / 6 * k4
    except Exception as e:
        print(e.with_traceback())
        sys.exit()
