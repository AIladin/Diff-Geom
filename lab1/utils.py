from lab1.Curve import *
import sympy as sym
sym.init_printing(use_unicode=False, wrap_line=True)


class EqBuilder:
    def __init__(self, curve: NaturalPCurve):
        self._curve = curve

    @staticmethod
    def _rib_eq_builder(s_0, vector: np.array, digits=5):
        x_0, y_0, z_0 = s_0
        return r"$\frac{x " + ("- " if x_0 > 0 else "+ ") + str(abs(x_0))[:digits] + " }{" +\
               str(vector[0])[:digits] + "}" + "=" +\
               r"\frac{y " + ("- " if y_0 > 0 else "+ ") + str(abs(y_0))[:digits] + " }{" +\
               str(vector[1])[:digits] + "}" + "=" + \
               r"\frac{y " + ("- " if z_0 > 0 else "+ ") + str(abs(z_0))[:digits] + " }{" +\
               str(vector[2])[:digits] + "}$"

    def tangent_eq(self, s):
        return self._rib_eq_builder(self._curve(s), tangent(self._curve, s))

    def normal_eq(self, s):
        return self._rib_eq_builder(self._curve(s), normal(self._curve, s))

    def binormal_eq(self, s):
        return self._rib_eq_builder(self._curve(s), binormal(self._curve, s))

    @staticmethod
    def _faces_eq_builder(s_0, vector: np.array):
        r = np.array((sym.Symbol('x'), sym.Symbol('y'), sym.Symbol('z')))
        f = np.sum((r-s_0)*vector)
        return str(sym.simplify(sym.N(f, 3))) + " = 0"

    def normal_plane(self, s):
        return self._faces_eq_builder(self._curve(s), tangent(self._curve, s))

    def osculating_plane(self, s):  # стична
        return self._faces_eq_builder(self._curve(s), binormal(self._curve, s))

    def rectifying_plane(self, s):  # спрямна
        return self._faces_eq_builder(self._curve(s), normal(self._curve, s))
