import sympy
from itertools import combinations
from sympy import pprint
sympy.init_printing(use_unicode=True)


class ImplicitFunctionSolver:

    def __init__(self, f, g, point: tuple, var=('x', 'y', 'z')):
        """
        :param f: sympy formula of a plane
        :param g: --//--
        :param point: P(x_0, y_0, z_0)
        """
        self.var = var
        assert f.subs(zip(var, point)) == 0 and g.subs(zip(var, point)) == 0
        self.F = f
        self.G = g
        self.point = point
        print(f"F{var} = {self.F}")
        print(f"G{var} = {self.G}")
        print()
        self.variable = None

    @classmethod
    def from_string(cls, f: str, g: str, point: tuple, var=('x', 'y', 'z')):
        """ Constructor for string arguments

        :param f: F(var) = 0
        :param g: G(var) = 0
        :param point: P(x_0, y_0, z_0)
        :param var: variables in formulas
        :return: ImplicitFunctionSolver class
        """
        sympy.var(var)
        return cls(sympy.simplify(f.split("=")[0]), sympy.simplify(g.split("=")[0]), point)

    def find_variable(self):
        print("Finding dependent variables.")
        for a, b in combinations(self.var, 2):

            matrix = sympy.Matrix([[sympy.diff(self.F, a), sympy.diff(self.F, b)],
                                   [sympy.diff(self.G, a), sympy.diff(self.G, b)]])
            det = matrix.det()
            val = det.subs(zip(self.var, self.point))
            pprint(matrix)
            print(f"det = {det} = {val}")
            if val != 0:
                var = (set(self.var) - {a, b}).pop()
                print(f"Accepted {var}.\n")
                return var
            print()

    def solve(self, steps=3):
        self.variable = self.find_variable()
        assert self.variable is not None

        print("Calculating parametrisation diff.")
        f = sympy.Function('f')(self.variable)
        g = sympy.Function('g')(self.variable)
        tmp = [f, g]
        tmp.insert(self.var.index(self.variable), self.variable)
        r = sympy.Matrix(tmp).transpose()

        print(f"r = ", end='')
        pprint(r)
        print()

        f_main = self.F.subs(zip(self.var, tmp))
        print(f_main, "= 0")
        g_main = self.G.subs(zip(self.var, tmp))
        print(g_main, "= 0")
        print()

        diffs = {f: self.point[tmp.index(f)],
                 g: self.point[tmp.index(g)]}

        var_sub = {sympy.symbols(self.variable): self.point[tmp.index(self.variable)]}
        rez = []
        for step in range(steps):  # TODO
            print("Step ", step+1)
            f_main = f_main.diff(self.variable)
            g_main = g_main.diff(self.variable)

            f_st = sympy.Function('f' + "'"*(step+1))(self.variable)
            g_st = sympy.Function('g' + "'"*(step+1))(self.variable)

            diff_sub = {f.diff(self.variable): f_st,
                        g.diff(self.variable): g_st}
            diffs.update(diff_sub)
            r = r.diff(self.variable).subs(diff_sub)

            print("r" + "'"*(step + 1) + " =", end='')
            pprint(r)

            f_exp = f_main.subs(diffs).subs(var_sub)
            g_exp = g_main.subs(diffs).subs(var_sub)
            print(f_main, "= 0")
            print(g_main, "= 0")
            print()

            print("Solving system.")
            # f_exp = f_exp.subs(diffs.items()).xreplace(var_sub)
            # g_exp = g_exp.subs(diffs.items()).xreplace(var_sub)

            solve = {key.subs(var_sub[sympy.symbols(self.variable)], sympy.symbols(self.variable)): value
                     for key, value in sympy.solve((f_exp, g_exp))[0].items()}

            # f_exp = f_exp.subs(var_sub[sympy.symbols(self.variable)], sympy.symbols(self.variable))
            # g_exp = g_exp.subs(var_sub[sympy.symbols(self.variable)], sympy.symbols(self.variable))

            print(f_exp, "= 0")
            print(g_exp, "= 0")
            print()

            diffs.update(solve)
            print("Solution ", solve)

            f = f_st
            g = g_st

            rez.append(r.subs(diffs.items()))
            print()
        return rez


if __name__ == '__main__':
    solver = ImplicitFunctionSolver.from_string("x + y**2 + 4*z**2  = 0",
                                                "x - 2*y - 8*z - 1 = 0",
                                                (-5, 1, -1))
    for d in solver.solve(1):
        print(d)

