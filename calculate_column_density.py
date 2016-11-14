#!/usr/bin/env python
import math
import numpy as np
import plotly.graph_objs as go

from scipy import integrate
from plotly.offline import plot
def _get_weight(d, dim=1):
    """Get the value of the kernel in units of smoothing length."""
    ## Springel 2005 Gadget2 paper Eqn 4 weighting function
    if d < 0.5:
        f = 1. - 6. * d ** 2 + 6. * d ** 3
    elif d < 1:
        f = 2. * (1. - d) ** 3
    else:
        f = 0.

    ## To get agreement with Price 2005 and Springel 2005 which has half the smoothing length definition
    f *= 2. ** dim

    ## From Price 2005 this is the normalisation constants given in Eqn 5
    if dim == 1:
        sigma = 2./3.
    elif dim == 2:
        sigma = 10./(7.*np.pi)
    elif dim == 3:
        sigma = 1./np.pi
    else:
        print("What dimensionality is your problem?")

    return f * sigma

def _Kernel(b, h=1., kernelspline=1., dim=1):
    """ Integrate the Kernel along z at distance b from origin """
    ## Distances are all in smoothing length 'h' and hence z integral runs from 0 to zmax_2D
    ## Integral is symmetric about zero, so (-zmax_2D to 0) + (0 to zmax_2D) is just double (0 to zmax_2D)
    ## Price's eqn 30 but modified to give
    zmax_2D = np.sqrt(kernelspline ** 2 - b ** 2)
    integrated_los = 2. * integrate.quad(lambda z: _get_weight(np.sqrt(z ** 2 + b ** 2), dim=dim), 0, zmax_2D)[0]

    return integrated_los / h ** dim

def testsph(x,y,h,xorigin=0.,yorigin=0.,dim=1):
    """ Given x,y coordinates and smoothing length get SPH weight along los """
    ## Work out impact parameter
    b = np.sqrt( (x-xorigin) ** 2 + (y - yorigin) ** 2)
    ## Put in units of smoothing length, zero contribution when further then 2h
    b /= h

    los_value = _Kernel(b, h=h, dim=dim)

    return los_value

def run_test(dim=1):
    """
        Create a bunch of different smoothing lengths

        For each smoothing length, create a set of pencilbeams
        at various distances from the center of the particle
        (ie different impact points)
        Calculate the weighted integral for each pencilbeam.

        Together these pencilbeams cover half the particle.
        Integrating their values should give us the full mass
        of the particle(?) which should be 1(?)
        This mass should be identical regardless of the
        smoothing length as increasing the smoothing only
        decreases the density (all other things being equal)
        """

    traces = []

    for smoothing in range(10, 101, 10):
        pencilbeams = []
        num_sight_lines = 1000

        # Construct our pencilbeams
        for i in range(0, num_sight_lines+1):
            # Make impact parameters covering the full
            # particle
            x = i / (1. * num_sight_lines) * smoothing

            pencilbeams.append(
                dict(x=x, y=0),
            )

        results = []
        for pencilbeam in pencilbeams:
            result = testsph(h=smoothing, dim=dim, **pencilbeam)
            results.append(result)

        # Integrate the pencilbeam weightings to find the full SPH weighting
        # This is the plane x-z from origin along +ve x-axis (sitting at y=0)
        particle_integral = integrate.trapz([x for x in results], [x['x'] for x in pencilbeams])
        # "All smoothing lengths should integrate to the same value "
        print(1 / (particle_integral ** 2 * math.pi ** 3))
        traces.append(go.Scatter(y=[x for x in results], x=[y['x'] for y in pencilbeams]))

    # The mass of a particle should be the area under each of these curves(?)
    plot(traces)


if __name__ == '__main__':
    run_test()