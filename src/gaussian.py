
import numpy

from scipy.special import gamma


def gaussian(r, p, a):
    # not normalized
    return sum(r**p)*numpy.exp(-a*numpy.dot(r, r))

def normalization(p, a):
    return (
        gamma(0.5*(p[0]+1))*
        gamma(0.5*(p[1]+1))*
        gamma(0.5*(p[2]+1))*
        a**(-1.5-0.5*sum(p))*(
            1 + (-1)**p[0] + (-1)**p[1] + (-1)**p[2] +
            (-1)**(p[0]+p[1]) + (-1)**(p[1]+p[2]) + (-1)**(p[2]+p[0]) +
            (-1)**(p[0]+p[1]+p[2])
    ))/8.0

