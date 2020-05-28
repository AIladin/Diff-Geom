from scipy.misc import derivative
from scipy.linalg import norm
import numpy as np


class NaturalPCurve:
    def __init__(self, x: callable, y: callable, z: callable):
        self._is_planar = z is None

        self._x = x
        self._y = y
        self._z = z

    def x(self, s):
        return self._x(s)

    def y(self, s):
        return self._y(s)

    def z(self, s):
        return self._z(s)

    def __call__(self, s: float):
        return np.array((self.x(s), self.y(s), self.z(s)))

    def deriv(self, s, n=1):
        d_x = derivative(self.x, s, n=n, order=7)
        d_y = derivative(self.y, s, n=n, order=7)
        d_z = derivative(self.z, s, n=n, order=7)
        return np.array((d_x, d_y, d_z))


def tangent(curve: NaturalPCurve, s):
    return curve.deriv(s)


def normal(curve: NaturalPCurve, s):
    n = curve.deriv(s, 2)
    return n/norm(n)


def binormal(curve: NaturalPCurve, s):
    return np.cross(tangent(curve, s), normal(curve, s))


def curvature(curve: NaturalPCurve, s):
    return norm(curve.deriv(s, 2))


def torsion(curve: NaturalPCurve, s):
    d = []
    for i in range(1, 4):
        d.append(curve.deriv(s, i))
    return np.dot(d[0], np.cross(d[1], d[2])).sum()/norm(d[1]) ** 2


def osculating_circle(curve: NaturalPCurve, s):
    radius = 1/curvature(curve, s)
    center = curve(s) + radius*normal(curve, s)
    return center, radius
