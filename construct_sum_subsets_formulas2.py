# Construct formulas for counting sum of subsets modulo div
# Problem related to youtube video "Olympiad level counting - How many subsets of {1,…,2000} have a sum divisible by 5"
#  by 3Blue1Brown (url: https://www.youtube.com/watch?v=bOXCLR3Wric)

import sympy as sp
#from sympy import Symbol
from sympy import symbols

from sympy.interactive.printing import init_printing

from calculate_sum_subsets import *
from calculate_sum_subsets_cs import *
from construct_sum_subsets_formulas import KeySymbol


class SumSubsetFormulaResult:
    def __init__(self, d, v, zeta):
        self.d = d
        self.divisors = sp.divisors(d)
        self.root = sp.root(1, d, 1)
        self.zeta = zeta

        self.sym_c = symbols('c:' + str(d))
        self.lookup_v_s = {}  # value -> symbol
        self.key_symbols = {}  # key symbol -> KeySymbol objects
        self.sym_c2 = list(self.sym_c)
        self.formulas = {}  # key symbol -> formula string

        # list of dictionary (index 0,...,d-1) with mapping "power of two" to coefficients
        self.list_dict = []
        for i in range(d):
            self.list_dict.append({})
            for j in range(d):
                if (v[j, 0] > 0):
                    self.list_dict[i][v[j, 0]] = 0

    def __shrink_poly_evalf(self, expr, sym, div, list_index, key):
        temp_expr = expr.subs(self.zeta, self.root)
        temp_expr = temp_expr.evalf()
        re = sp.re(temp_expr)
        im = sp.im(temp_expr)
        tolerance = 1e-8
        assert (abs(im) < tolerance)
        re_rounded = round(re)

        if (abs(re - re_rounded) >= tolerance):
            return None
        return sp.Poly([re_rounded], sym)

    def __simplify(self, list_index):
        for key, value in self.list_dict[list_index].items():
            # Eliminate roots in polynomial myself since Sympy seems to be unable to eliminate roots in sum.
            value = shrink_poly(value, self.zeta, self.d,
                                self.divisors, frantic=True)
            # shrink_poly should have already eliminated all roots. Unfortunately, it doesn't. Even when frantic=True :-(
            # It might happen if roots of unity are distributed for different power of twos! E. g.:
            # d=9, list_index=1, key=2, value=Poly(zeta**8 + zeta**7 + zeta**5 + zeta**4 + zeta**2 + zeta, zeta, domain='ZZ')
            # {512: 1, 2: zeta**8 + zeta**7 + zeta**5 + zeta**4 + zeta**2 + zeta, 8: zeta**6 + zeta**3}
            # In spite of frantic=True for shrink_poly, assertion occurs for d=75!
            # d=75, list_index=1, key=32, value=Poly(zeta**70 + zeta**65 + zeta**55 + zeta**40 + zeta**35 + zeta**20 + zeta**10 + zeta**5, zeta, domain='ZZ')
            # {37778931862957161709568: 1, 2: 0, 8: 0, 32: zeta**70 + zeta**65 + zeta**55 + zeta**40 + zeta**35 + zeta**20 + zeta**10 + zeta**5, 32768: zeta**60 + zeta**45 + zeta**30 + zeta**15, 33554432: zeta**50 + zeta**25}
            if (len(value.free_symbols) != 0):
                # Evaluate polynom as floating point and and convert it back to SymPy Poly with integer coefficients.
                new_value = self.__shrink_poly_evalf(
                    value, self.zeta, self.d, list_index, key)
                if new_value is not None:
                    value = new_value
                else:
                    print(
                        f"# d={self.d}, list_index={list_index}, key={key}, value={value}")
                    print(self.list_dict[list_index])
                    assert False, "Cannot eliminate free symbols!"
            value = value.subs(self.zeta, self.root)
            self.list_dict[list_index][key] = value

    def simplify_all(self):
        #two = sp.symbols("two")
        #m = sp.symbols("m")
        # for i in range(self.d):
        #    sum = 0
        #    for key, value in self.list_dict[i].items():
        #        e, log_result = sp.integer_log(key, 2)
        #        assert(log_result)
        #        sum += two**(e*m) * value
        #    print(f"sum for i={i}:", sum)

        for i in range(self.d):
            self.__simplify(i)
        for i in range(self.d):
            key = frozenset(self.list_dict[i].items())
            if key in self.lookup_v_s:
                self.sym_c2[i] = self.lookup_v_s[key]
                self.key_symbols[self.sym_c2[i]].add_symbol(self.sym_c[i])
            else:
                self.lookup_v_s[key] = self.sym_c[i]
                self.key_symbols[self.sym_c[i]] = KeySymbol(self.sym_c[i], i)
        self.__construct_formulas()

    def __construct_formulas(self):
        sorted_dict = {}
        for key, value in self.lookup_v_s.items():
            for entry in key:
                sorted_dict[entry[0]] = entry[1]
            sorted_dict = dict(sorted(sorted_dict.items()))
            formula_numerator = ""
            for key2, value2 in sorted_dict.items():
                if value2 == 0:
                    continue
                if len(formula_numerator) > 0:
                    formula_numerator += (" + " if value2 > 0 else " - ")
                    formula_numerator += (f"{value2}" if value2 >
                                          0 else f"{-value2}")
                else:
                    formula_numerator += f"{value2}"
                formula_numerator += " * "
                e, log_result = sp.integer_log(key2, 2)
                assert (log_result)
                if e == 1:
                    formula_numerator += f"2**(n//{self.d})"
                else:
                    formula_numerator += f"(2**{e})**(n//{self.d})"
            formula_string = f"( {formula_numerator} ) // {self.d}"
            self.formulas[value] = formula_string

    def __str__(self):
        str = self.key_symbols.__str__() + "\nFormulas:\n"
        for key_symbol, formula in self.formulas.items():
            str += f"{key_symbol}(n) = {formula}\n"
        return str

    def print_code(self):
        #indentation = 2
        # self.formula_objects[10] = SumSubsetsFormula(10, [0, 1, 1, 1, 1, 0, 1, 1, 1, 1],
        #    [0, lambda n: (( 4 * (2**2)**(n//10) + 1 * (2**10)**(n//10) ) // 10)],  # c0_10
        #    [1, lambda n: (( -1 * (2**2)**(n//10) + 1 * (2**10)**(n//10) ) // 10)],  # c1_10
        #    )
        str = ""
        for i, symbol in enumerate(self.sym_c2):
            if i == 0:
                str += f"        self.formula_objects[{self.d}] = SumSubsetsFormula({self.d},  ["
            else:
                str += ", "
            str += f"{self.key_symbols[symbol].index}"
        str += "],"
        print(str)
        for key_symbol, formula in self.formulas.items():
            print(
                f"            [{self.key_symbols[key_symbol].index}, lambda n: ({formula})],  # {key_symbol}_{self.d}")
        print("            )")


def construct_sum_subsets_formula(div):
    if div < 2:
        raise ValueError("n must be >= 0 and div must be >=2!")

    d = div
    v = sp.Matrix(d, 1, lambda row, col: 1 if row == 0 else 0)

    v = calculate_product_matrix(d)
    for i in range(d):
        if v[i, 0] > 0:
            e, log_result = sp.integer_log(v[i, 0], 2)
            assert log_result, f"Entry {i} in product matrix is not power of 2!"

    # factors to generating functions
    # row: index for functions, col: index of coefficient
    zeta = sp.symbols("zeta")
    gf_factors = sp.Matrix(d, d, lambda row, col: zeta ** (((d-row)*col) % d))

    # v2 = v.T * gf_factors
    result = SumSubsetFormulaResult(d, v, zeta)
    for i in range(d):
        for j in range(d):
            if (v[j, 0] > 0):
                result.list_dict[i][v[j, 0]] += gf_factors[j, i]
    result.simplify_all()
    # result.print_code()  # output code
    print(f"##### d: {d}\n{result}")  # OR output normally


def test_expressions():
    z = sp.symbols("z")
    expr = sp.sympify("z**14 + z**13 + z**11 + z**8 + z**7 + z**4 + z**2 + z")
    print(f"div=15,", expr, "=>", expr.subs(z, sp.root(1, 15, 1)).evalf())
    expr = sp.sympify(
        "z**20 + z**19 + z**17 + z**16 + z**13 + z**11 + z**10 + z**8 + z**5 + z**4 + z**2 + z")
    print(f"div=21,", expr, "=>", expr.subs(z, sp.root(1, 21, 1)).evalf())
    expr = sp.sympify(
        "z**28 + z**26 + z**22 + z**16 + z**14 + z**8 + z**4 + z**2")
    print(f"div=30,", expr, "=>", expr.subs(z, sp.root(1, 30, 1)).evalf())
    expr = sp.sympify(
        "z**32 + z**31 + z**29 + z**28 + z**26 + z**25 + z**23 + z**20 + z**19 + z**17 + z**16 + z**14 + z**13 + z**10 + z**8 + z**7 + z**5 + z**4 + z**2 + z")
    print(f"div=33,", expr, "=>", expr.subs(z, sp.root(1, 33, 1)).evalf())
    expr = sp.sympify("z**34 + z**33 + z**32 + z**31 + z**29 + z**27 + z**26 + z**24 + z**23 + z**22 + z**19 + z**18 + z**17 + z**16 + z**13 + z**12 + z**11 + z**9 + z**8 + z**6 + z**4 + z**3 + z**2 + z")
    print(f"div=35,", expr, "=>", expr.subs(z, sp.root(1, 35, 1)).evalf())
    expr = sp.sympify("z**38 + z**37 + z**35 + z**34 + z**32 + z**31 + z**29 + z**28 + z**25 + z**23 + z**22 + z**20 + z**19 + z**17 + z**16 + z**14 + z**11 + z**10 + z**8 + z**7 + z**5 + z**4 + z**2 + z")
    print(f"div=39,", expr, "=>", expr.subs(z, sp.root(1, 39, 1)).evalf())
    expr = sp.sympify(
        "z**40 + z**38 + z**34 + z**32 + z**26 + z**22 + z**20 + z**16 + z**10 + z**8 + z**4 + z**2")
    print(f"div=42,", expr, "=>", expr.subs(z, sp.root(1, 42, 1)).evalf())
    expr = sp.sympify(
        "z**42 + z**39 + z**33 + z**24 + z**21 + z**12 + z**6 + z**3")
    print(f"div=45,", expr, "=>", expr.subs(z, sp.root(1, 45, 1)).evalf())
    expr = sp.sympify(
        "z**70 + z**65 + z**55 + z**40 + z**35 + z**20 + z**10 + z**5")
    print(f"div=75,", expr, "=>", expr.subs(z, sp.root(1, 75, 1)).evalf())
    print("")


def test_formulas():
    c9_n = []

    # Formulas:
    # c0(n) = ( 6 * 2**(n//9) + 2 * (2**3)**(n//9) + 1 * (2**9)**(n//9) ) // 9
    # c1(n) = ( -1 * (2**3)**(n//9) + 1 * (2**9)**(n//9) ) // 9
    # c3(n) = ( -3 * 2**(n//9) + 2 * (2**3)**(n//9) + 1 * (2**9)**(n//9) ) // 9

    c9_n.append(sp.sympify(
        "( 6 * 2**(n//9) + 2 * (2**3)**(n//9) + 1 * (2**9)**(n//9) ) // 9"))
    c9_n.append(sp.sympify("( -1 * (2**3)**(n//9) + 1 * (2**9)**(n//9) ) // 9"))
    c9_n.append(sp.sympify(
        "( -3 * 2**(n//9) + 2 * (2**3)**(n//9) + 1 * (2**9)**(n//9) ) // 9"))
    n = sp.symbols("n")
    for i in range(1, 6):
        v = i * 9
        print("n", v, "c0(n)", c9_n[0].subs(n, v))
        print("n", v, "c1(n)", c9_n[1].subs(n, v))
        print("n", v, "c3(n)", c9_n[2].subs(n, v))
        print("csslog: ", calcluate_sum_subsets_logarithmic(v, 9))
        print("cssgcd: ", calcluate_sum_subsets_gcd(v, 9))


def main():
    init_printing(use_unicode=False, wrap_line=False, num_columns=300)

    # test_expressions()
    for i in range(2, 101):
        try:
            construct_sum_subsets_formula(i)
        except AssertionError as ae:
            print(f"##### d: {i}. Assertion occured! ({ae.__str__()})")

    # test_formulas()


if __name__ == '__main__':
    main()
