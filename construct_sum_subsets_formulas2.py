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
        self.lookup_v_s = {} # value -> symbol
        self.key_symbols = {} # key symbol -> KeySymbol objects
        self.sym_c2 = list(self.sym_c)

        self.list_dict =  [] # list of dictionary (index 0,...,d-1) with mapping "power of two" to coefficients
        for i in range (d):
            self.list_dict.append({})
            for j in range (d):
                if(v[j,0] > 0):
                    self.list_dict[i][v[j,0]] = 0

    def shrink_poly_cheat(self, expr, sym, div, list_index, key):
        all_coeffs = sp.poly(expr, sym).all_coeffs()
        set_length = None
        significant_indices_correct = None
        possible_values = None
        value_to_check = None
        if div == 15:
            # "zeta**14 + zeta**13 + zeta**11 + zeta**8 + zeta**7 + zeta**4 + zeta**2 + zeta" => 1
            set_length = 15
            significant_indices_correct = [1, 2, 4, 7, 8, 11, 13, 14]
            possible_values = [1]
        elif div == 21:
            # "zeta**20 + zeta**19 + zeta**17 + zeta**16 + zeta**13 + zeta**11 + zeta**10 + zeta**8 + zeta**5 + zeta**4 + zeta**2 + zeta" => 1
            set_length = 21
            significant_indices_correct = [1, 2, 4, 5, 8, 10, 11, 13, 16, 17, 19, 20]
            possible_values = [1]
        elif div == 30:
            # "zeta**28 + zeta**26 + zeta**22 + zeta**16 + zeta**14 + zeta**8 + zeta**4 + zeta**2" => 1
            set_length = 29
            significant_indices_correct = [2, 4, 8, 14, 16, 22, 26, 28]
            possible_values = [1]
        elif div == 33:
            # "zeta**32 + zeta**31 + zeta**29 + zeta**28 + zeta**26 + zeta**25 + zeta**23 + zeta**20 + zeta**19 + zeta**17 + zeta**16 + zeta**14 "
            # "+ zeta**13 + zeta**10 + zeta**8 + zeta**7 + zeta**5 + zeta**4 + zeta**2 + zeta" => 1
            set_length = 33
            significant_indices_correct = [1, 2, 4, 5, 7, 8, 10, 13, 14, 16, 17, 19, 20, 23, 25, 26, 28, 29, 31, 32]
            possible_values = [1]
        elif div == 35:
            # "zeta**34 + zeta**33 + zeta**32 + zeta**31 + zeta**29 + zeta**27 + zeta**26 + zeta**24 + zeta**23 + zeta**22 + zeta**19 + zeta**18 "
            # "+ zeta**17 + zeta**16 + zeta**13 + zeta**12 + zeta**11 + zeta**9 + zeta**8 + zeta**6 + zeta**4 + zeta**3 + zeta**2 + zeta" => 1
            set_length = 35
            significant_indices_correct = [1, 2, 3, 4, 6, 8, 9, 11, 12, 13, 16, 17, 18, 19, 22, 23, 24, 26, 27, 29, 31, 32, 33, 34]
            possible_values = [1]
        elif div == 39:
            # "zeta**38 + zeta**37 + zeta**35 + zeta**34 + zeta**32 + zeta**31 + zeta**29 + zeta**28 + zeta**25 + zeta**23 + zeta**22 + zeta**20 "
            # "+ zeta**19 + zeta**17 + zeta**16 + zeta**14 + zeta**11 + zeta**10 + zeta**8 + zeta**7 + zeta**5 + zeta**4 + zeta**2 + zeta" => 1
            set_length = 39
            significant_indices_correct = [1, 2, 4, 5, 7, 8, 10, 11, 14, 16, 17, 19, 20, 22, 23, 25, 28, 29, 31, 32, 34, 35, 37, 38]
            possible_values = [1]
        elif div == 42:
            # "zeta**40 + zeta**38 + zeta**34 + zeta**32 + zeta**26 + zeta**22 + zeta**20 + zeta**16 + zeta**10 + zeta**8 + zeta**4 + zeta**2" => 1
            set_length = 41
            significant_indices_correct = [2, 4, 8, 10, 16, 20, 22, 26, 32, 34, 38, 40]
            possible_values = [1]
        elif div == 45:
            # "zeta**42 + zeta**39 + zeta**33 + zeta**24 + zeta**21 + zeta**12 + zeta**6 + zeta**3" => 1
            # "3*zeta**42 + 3*zeta**39 + 3*zeta**33 + 3*zeta**24 + 3*zeta**21 + 3*zeta**12 + 3*zeta**6 + 3*zeta**3" => 3
            set_length = 43
            significant_indices_correct = [3, 6, 12, 21, 24, 33, 39, 42]
            possible_values = [1, 3]
        elif div == 75:
            # "zeta**70 + zeta**65 + zeta**55 + zeta**40 + zeta**35 + zeta**20 + zeta**10 + zeta**5" => 1
            # "5*zeta**70 + 5*zeta**65 + 5*zeta**55 + 5*zeta**40 + 5*zeta**35 + 5*zeta**20 + 5*zeta**10 + 5*zeta**5" => 5
            set_length = 71
            significant_indices_correct = [5, 10, 20, 35, 40, 55, 65, 70]    
            possible_values = [1, 5]
        else:
            return None
        
        assert (len(all_coeffs) == set_length), f"Mismatch! len(all_coeffs): {len(all_coeffs)}. Shall be {set_length}!"
        
        significant_indices = [(set_length - x - 1) for x in significant_indices_correct]
        
        if all_coeffs[significant_indices[0]] not in possible_values:
            return None
        
        value_to_check = all_coeffs[significant_indices[0]]

        for i, coeff in enumerate(all_coeffs):
            if i in significant_indices:
                if coeff != value_to_check:
                    return None
            else:
                if coeff != 0:
                    return None

        #print(f"div: {div}, list_index: {list_index}, key: {key}, new value: {value_to_check}")
        return sp.Poly([value_to_check], sym)

    def simplify(self, list_index):
        for key, value in self.list_dict[list_index].items():
            # Eliminate roots in polynomial myself since Sympy seems to be unable to eliminate roots in sum.
            value = shrink_poly(value, self.zeta, self.d, self.divisors, frantic=True)
            # shrink_poly should have already eliminated all roots. Unfortunately, it doesn't. :-(
            if(len(value.free_symbols)!=0):
                # It might happen if roots of unity are distributed for different power of twos! E. g.:
                # d=9, list_index=1, key=2, value=Poly(zeta**8 + zeta**7 + zeta**5 + zeta**4 + zeta**2 + zeta, zeta, domain='ZZ')
                # {512: 1, 2: zeta**8 + zeta**7 + zeta**5 + zeta**4 + zeta**2 + zeta, 8: zeta**6 + zeta**3}
                # In spite of frantic=True for shrink_poly, assertion occurs for d=75!
                # d=75, list_index=1, key=32, value=Poly(zeta**70 + zeta**65 + zeta**55 + zeta**40 + zeta**35 + zeta**20 + zeta**10 + zeta**5, zeta, domain='ZZ')
                # {37778931862957161709568: 1, 2: 0, 8: 0, 32: zeta**70 + zeta**65 + zeta**55 + zeta**40 + zeta**35 + zeta**20 + zeta**10 + zeta**5, 32768: zeta**60 + zeta**45 + zeta**30 + zeta**15, 33554432: zeta**50 + zeta**25}
                
                # For some values of d, we know the polynomials to be returned and assign exact values for them.
                new_value = self.shrink_poly_cheat(value, self.zeta, self.d, list_index, key)
                if new_value is not None:
                    value = new_value
                else:
                    print(f"# d={self.d}, list_index={list_index}, key={key}, value={value}")
                    print(self.list_dict[list_index])
                    assert False, "Cannot eliminate free symbols!"
            value = value.subs(self.zeta, self.root)
            self.list_dict[list_index][key] = value

    def simplify_all(self):
        #two = sp.symbols("two")
        #m = sp.symbols("m")
        #for i in range(self.d):
        #    sum = 0
        #    for key, value in self.list_dict[i].items():
        #        e, log_result = sp.integer_log(key, 2)
        #        assert(log_result)
        #        sum += two**(e*m) * value
        #    print(f"sum for i={i}:", sum)

        for i in range(self.d):
            self.simplify(i)
        for i in range(self.d):
            key = frozenset(self.list_dict[i].items())
            if key in self.lookup_v_s:
                self.sym_c2[i] = self.lookup_v_s[key]
                self.key_symbols[self.sym_c2[i]].add_symbol(self.sym_c[i])
            else:
                self.lookup_v_s[key] = self.sym_c[i]
                self.key_symbols[self.sym_c[i]] = KeySymbol(self.sym_c[i], i)


    def __str__(self):
        str = self.key_symbols.__str__() + "\nFormulas:\n"
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
                    formula_numerator += (f"{value2}" if value2 > 0 else f"{-value2}")
                else:
                    formula_numerator += f"{value2}"
                formula_numerator += " * "
                e, log_result = sp.integer_log(key2, 2)
                assert(log_result)
                if e == 1:
                    formula_numerator += f"2**(n//{self.d})"
                else:
                    formula_numerator += f"(2**{e})**(n//{self.d})"
            formula_string = f"( {formula_numerator} ) // {self.d}"
            str += f"{value}(n) = {formula_string}\n"

        return str


def construct_sum_subsets_formula(div):
    if div < 2:
        raise ValueError("n must be >= 0 and div must be >=2!")
    
    d = div
    v = sp.Matrix(d, 1, lambda row, col: 1 if row == 0 else 0)

    v = calculate_product_matrix(d)
    for i in range(d):
        if v[i, 0] > 0:
            e, log_result = sp.integer_log(v[i,0], 2)
            assert log_result, f"Entry {i} in product matrix is not power of 2!"

    # factors to generating functions
    # row: index for functions, col: index of coefficient
    zeta = sp.symbols("zeta")
    gf_factors = sp.Matrix(d, d, lambda row, col: zeta ** (((d-row)*col) % d))

    # v2 = v.T * gf_factors
    result = SumSubsetFormulaResult(d, v, zeta)
    for i in range (d):
        for j in range (d):
            if(v[j,0] > 0):
                result.list_dict[i][v[j,0]] += gf_factors[j,i]
    result.simplify_all()

    print(f"##### d: {d}\n{result}")

def test_expressions():
    z = sp.symbols("z")
    expr = sp.sympify("z**14 + z**13 + z**11 + z**8 + z**7 + z**4 + z**2 + z")
    print(f"div=15,", expr, "=>", expr.subs(z, sp.root(1, 15, 1)).evalf())
    expr = sp.sympify("z**20 + z**19 + z**17 + z**16 + z**13 + z**11 + z**10 + z**8 + z**5 + z**4 + z**2 + z")
    print(f"div=21,", expr, "=>", expr.subs(z, sp.root(1, 21, 1)).evalf())
    expr = sp.sympify("z**28 + z**26 + z**22 + z**16 + z**14 + z**8 + z**4 + z**2")
    print(f"div=30,", expr, "=>", expr.subs(z, sp.root(1, 30, 1)).evalf())
    expr = sp.sympify("z**32 + z**31 + z**29 + z**28 + z**26 + z**25 + z**23 + z**20 + z**19 + z**17 + z**16 + z**14 + z**13 + z**10 + z**8 + z**7 + z**5 + z**4 + z**2 + z")
    print(f"div=33,", expr, "=>", expr.subs(z, sp.root(1, 33, 1)).evalf())
    expr = sp.sympify("z**34 + z**33 + z**32 + z**31 + z**29 + z**27 + z**26 + z**24 + z**23 + z**22 + z**19 + z**18 + z**17 + z**16 + z**13 + z**12 + z**11 + z**9 + z**8 + z**6 + z**4 + z**3 + z**2 + z")
    print(f"div=35,", expr, "=>", expr.subs(z, sp.root(1, 35, 1)).evalf())
    expr = sp.sympify("z**38 + z**37 + z**35 + z**34 + z**32 + z**31 + z**29 + z**28 + z**25 + z**23 + z**22 + z**20 + z**19 + z**17 + z**16 + z**14 + z**11 + z**10 + z**8 + z**7 + z**5 + z**4 + z**2 + z")
    print(f"div=39,", expr, "=>", expr.subs(z, sp.root(1, 39, 1)).evalf())
    expr = sp.sympify("z**40 + z**38 + z**34 + z**32 + z**26 + z**22 + z**20 + z**16 + z**10 + z**8 + z**4 + z**2")
    print(f"div=42,", expr, "=>", expr.subs(z, sp.root(1, 42, 1)).evalf())
    expr = sp.sympify("z**42 + z**39 + z**33 + z**24 + z**21 + z**12 + z**6 + z**3")
    print(f"div=45,", expr, "=>", expr.subs(z, sp.root(1, 45, 1)).evalf())
    expr = sp.sympify("z**70 + z**65 + z**55 + z**40 + z**35 + z**20 + z**10 + z**5")
    print(f"div=75,", expr, "=>", expr.subs(z, sp.root(1, 75, 1)).evalf())
    print("")

def test_formulas():
    c9_n = []

    #Formulas:
    # c0(n) = ( 6 * 2**(n//9) + 2 * (2**3)**(n//9) + 1 * (2**9)**(n//9) ) // 9
    # c1(n) = ( -1 * (2**3)**(n//9) + 1 * (2**9)**(n//9) ) // 9
    # c3(n) = ( -3 * 2**(n//9) + 2 * (2**3)**(n//9) + 1 * (2**9)**(n//9) ) // 9

    c9_n.append(sp.sympify("( 6 * 2**(n//9) + 2 * (2**3)**(n//9) + 1 * (2**9)**(n//9) ) // 9"))
    c9_n.append(sp.sympify("( -1 * (2**3)**(n//9) + 1 * (2**9)**(n//9) ) // 9"))
    c9_n.append(sp.sympify("( -3 * 2**(n//9) + 2 * (2**3)**(n//9) + 1 * (2**9)**(n//9) ) // 9"))
    n = sp.symbols("n")
    for i in range(1, 6):
        v = i* 9
        print("n", v, "c0(n)", c9_n[0].subs(n, v))
        print("n", v, "c1(n)", c9_n[1].subs(n, v))
        print("n", v, "c3(n)", c9_n[2].subs(n, v))
        print("csslog: ", calcluate_sum_subsets_logarithmic(v, 9))
        print("cssgcd: ", calcluate_sum_subsets_gcd(v, 9))

def main():
    init_printing(use_unicode=False, wrap_line=False, num_columns=300)
    
    #test_expressions()
    for i in range(2, 51):
        try:
            construct_sum_subsets_formula(i)
        except AssertionError as ae:
            print(f"##### d: {i}. Assertion occured! ({ae.__str__()})")    

    # test_formulas()


if __name__ == '__main__':
    main()

