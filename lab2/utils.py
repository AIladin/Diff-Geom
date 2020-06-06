import numpy as np
import sympy as sym


def tangent(d1):
    return d1/np.linalg.norm(d1.astype('f'))


def normal(d1, d2):
    return np.cross(binormal(d1, d2), tangent(d1))


def binormal(d1, d2):
    prod = np.cross(d1, d2)
    return prod / np.linalg.norm(prod)


def curvature(d1, d2):
    return np.linalg.norm(np.cross(d1, d2)) / np.linalg.norm(d1) ** 3


def torsion(d1, d2, d3):
    return np.dot(d1, np.cross(d2, d3)).sum() / np.linalg.norm(np.cross(d1, d2)) ** 2


def osculating_circle(p, d1, d2):
    radius = 1/curvature(d1, d2)
    center = p + radius*normal(d1, d2)
    return center, radius


class EqBuilder:
    def __init__(self, p,  d1, d2, d3):
        self._p = p
        self._d1 = d1
        self._d2 = d2
        self._d3 = d3

    @staticmethod
    def _rib_eq_builder(s_0, vector: np.array, digits=5):
        x_0, y_0, z_0 = s_0
        return r"$\frac{x " + ("- " if x_0 > 0 else "+ ") + str(abs(x_0))[:digits] + " }{" +\
               str(vector[0])[:digits] + "}" + "=" +\
               r"\frac{y " + ("- " if y_0 > 0 else "+ ") + str(abs(y_0))[:digits] + " }{" +\
               str(vector[1])[:digits] + "}" + "=" + \
               r"\frac{y " + ("- " if z_0 > 0 else "+ ") + str(abs(z_0))[:digits] + " }{" +\
               str(vector[2])[:digits] + "}$"

    def tangent_eq(self):
        return self._rib_eq_builder(self._p, tangent(self._d1))

    def normal_eq(self):
        return self._rib_eq_builder(self._p, normal(self._d1, self._d2))

    def binormal_eq(self):
        return self._rib_eq_builder(self._p, binormal(self._d1, self._d2))

    @staticmethod
    def _faces_eq_builder(s_0, vector: np.array):
        r = np.array((sym.Symbol('x'), sym.Symbol('y'), sym.Symbol('z')))
        f = np.sum((r-s_0)*vector)
        return str(sym.simplify(sym.N(f, 3))) + " = 0"

    def normal_plane(self):
        return self._faces_eq_builder(self._p, tangent(self._d1))

    def osculating_plane(self):  # стична
        return self._faces_eq_builder(self._p, binormal(self._d1, self._d2))

    def rectifying_plane(self):  # спрямна
        return self._faces_eq_builder(self._p, normal(self._d1, self._d2))

    def osculating_circle(self, digits=5):
        center, radius = osculating_circle(self.p, self.d1, self.d2)
        return f"(x - {str(center[0])[:digits]})**2 + " \
               f"(y - {str(center[1])[:digits]})**2 + " \
               f"(z - {str(center[2])[:digits]})**2 = {radius**2}\n"\
               + self.osculating_plane()

    @property
    def p(self):
        return np.copy(self._p)

    @property
    def d1(self):
        return np.copy(self._d1)

    @property
    def d2(self):
        return np.copy(self._d2)

    @property
    def d3(self):
        return np.copy(self._d3)
