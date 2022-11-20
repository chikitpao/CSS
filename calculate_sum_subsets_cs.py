# Construct formulars for counting sum of subsets modulo div via cyclic sieving 
# Problem related to youtube video "Olympiad level counting - How many subsets of {1,…,2000} have a sum divisible by 5"
#  by 3Blue1Brown (url: https://www.youtube.com/watch?v=bOXCLR3Wric)

import sympy as sp

from sympy.interactive.printing import init_printing
from sympy.abc import x
import time

from calculate_sum_subsets import *


def shrink_poly(expr, sym, div, divisors):
    new_expr = sp.poly(expr, sym)
    
    # Replace power with (power % div)
    old_coeffs = new_expr.all_coeffs()
    new_coeffs = [0] * div
    old_coeffs_len = len(old_coeffs)
    for i in range(old_coeffs_len):
        new_coeffs[div - 1 - (i % div)] += old_coeffs[old_coeffs_len - 1 - i]
    
    # Cancel out roots using divsors of div.
    for p in divisors:
        if p == div:
            continue
        common_coefficient = None
        for i in range(0, div, p):
            index = div - i - 1
            if common_coefficient is None:
                common_coefficient = new_coeffs[index]
            elif common_coefficient > new_coeffs[index]:
                common_coefficient = new_coeffs[index]
        #print(p, common_coefficient)
        if common_coefficient is not None and common_coefficient > 0:
                for i in range(0, div, p):
                    index = div - i - 1
                    new_coeffs[index] -= common_coefficient

    new_poly = sp.Poly(new_coeffs, sym)

    return new_poly


class CyclicSievingMatrix:
    def __init__(self, d, debug=False):
        self.d = d
        self.debug = debug
        self.zeta = sp.symbols("zeta")
        # M: Matrix of zetas in as functions of sum 
        # row: index for functions, col: index of coefficient
        self.M = sp.Matrix(d, d, lambda row, col: self.zeta**((row * (col + 1)) % d))
        #print("# Matrix:", self.M)
        self.list_of_roots = [sp.root(1, d, k) for k in range(d)]
        #root = self.list_of_roots[1]
        #self.prime_factors = sp.ntheory.factorint(d)
        #print(f"prime factors of {d}: {self.prime_factors}")
        self.divisors = sp.divisors(d)
        if debug:
            print(f"divisors of {d}: {self.divisors}")

 
        # factors to functions
        # row: index for functions, col: index of coefficient
        self.factors = sp.Matrix(d, d, lambda row, col: self.zeta ** (((d-row)*col) % d))

    def print_sum_of_columns(self):
        # row: index for functions, col: index of coefficient
        v = self.M.T * self.factors
        for i in range(0, self.d):
            for j in range(0, self.d):
                v[i, j] = sp.expand(v[i, j])
                v[i, j] = shrink_poly(v[i, j], self.zeta, self.d, self.divisors)
                 # shrink_poly shall have already eliminated all roots.
                v[i, j] = v[i, j].subs(self.zeta,  self.list_of_roots[1])
        for i in range(0, self.d):
            print(f"sum of columns for i={i}: {v[:, i]}")
    
    def return_product_matrix(self):
        # row: index for functions
        v = sp.zeros(self.d, 1)        
        for i in range(0, self.d):
            product = 1
            for j in range(0, self.d):
                product = sp.expand(product * (self.M[i, j] + 1))
            #print(f"i={i}: product={product}, free_symbols: {product.free_symbols}")
            # Eliminate roots in polynomial myself since Sympy seems to be unable to eliminate roots in sum.
            new_poly = shrink_poly(product, self.zeta, self.d, self.divisors)
            # shrink_poly shall have already eliminated all roots.
            assert(len(new_poly.free_symbols)==0)
            v[i] = new_poly.subs(self.zeta,  self.list_of_roots[1])
            # print(f"product of row for i={i}: {v[i]}")
            #v[i] = v[i].evalf()
        if self.debug:
            print(f"product of rows: {v}")
        return v


def calcluate_sum_subsets_cs(n, div, debug=False):
    if n < 0 or div < 2:
        raise ValueError("n must be >= 0 and div must be >=2!")
    
    d = div
    v = sp.Matrix(div, 1, lambda row, col: 1 if row == 0 else 0)
    if n <= 0:
        return v

    remainder = n % d
    temp_n = n - remainder
    quotient = temp_n // d

    M = CyclicSievingMatrix(div, debug)
    v = M.return_product_matrix()
    if debug:
        M.print_sum_of_columns()
        print("factors for generating functions:", M.factors)
    for i in range(div):
        v[i,0] = v[i,0]**quotient

    # Eliminate roots in polynomial myself since Sympy seems to be unable to eliminate roots in sum.
    v2 = v.T * M.factors
    for i in range(div):
        #print(f"i={i}: v2[0,i].free_symbols: {v2[0,i].free_symbols}")
        v2[0,i] = shrink_poly(v2[0,i], M.zeta, M.d, M.divisors)
        # shrink_poly shall have already eliminated all roots.
        assert(len(v2[0,i].free_symbols)==0)
        v2[0,i] = v2[0,i].subs(M.zeta, M.list_of_roots[1])
        assert(v2[0,i] % d == 0)
        v2[0,i] = v2[0,i] // d
    v = v2.T

    if remainder == 0:
        return v

    M = create_matrix_from_solution_vector_sp(v)
    r = generate_solution_vectors_sp(d)
    v = M * r[remainder-1]
    return v

def test_calculation(n, div):
    print(f"##### div: {div}, n: {n}")
    t = time.process_time()
    r1 = calcluate_sum_subsets_cs(n, div, debug=True)
    time_log = time.process_time() - t
    print(f"# calcluate_sum_subsets_cs: {r1}")
    print(f"# time elapsed={time_log}")
    r2 = calcluate_sum_subsets_logarithmic(n, div)
    print(f"# Difference to result from calcluate_sum_subsets_logarithmic: {(r1-r2).evalf()}")


def test_roots_of_unity():
    expressions = [ "(-1)**(3/5)", "-(-1)**(3/5)", "(-1)**(1/5)", "(-1)**(1/3)", "(-1)**(2/3)", "(-1)**(4/5)", "(-1)**(2/5)"]
    for e in expressions:
        sympy_expr = sp.sympify(e)
        # Output examples:
        # "(-1)**(3/5), arg: atan(sqrt(sqrt(5)/8 + 5/8)/(1/4 - sqrt(5)/4)) + pi, 1.88495559215388"
        # "-(-1)**(3/5), arg: -atan(sqrt(sqrt(5)/8 + 5/8)/(-1/4 + sqrt(5)/4)), -1.25663706143592"
        # "(-1)**(1/5), arg: atan(sqrt(5/8 - sqrt(5)/8)/(1/4 + sqrt(5)/4)), 0.628318530717959"
        # "(-1)**(1/3), arg: pi/3, 1.04719755119660"
        print(f"{e}, arg: {sp.arg(sympy_expr)}, {sp.arg(sympy_expr).evalf()}" )

    root = sp.root(1, 5, 1)
    sum = 0
    for i in range(1, 11):
        t = root**i
        sum += t
        print(f"i: {i}, root**i: {t}, sum: {sum}")

   
def main():
    init_printing(use_unicode=True)
    # test_roots_of_unity()
    test_calculation(2000, 5)


if __name__ == '__main__':
    main()