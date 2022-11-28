# Construct formulas for counting sum of subsets modulo div via cyclic sieving 
# Problem related to youtube video "Olympiad level counting - How many subsets of {1,…,2000} have a sum divisible by 5"
#  by 3Blue1Brown (url: https://www.youtube.com/watch?v=bOXCLR3Wric)

import sympy as sp

from sympy.interactive.printing import init_printing
from sympy.abc import x
import time

from calculate_sum_subsets import *


def shrink_poly(expr, sym, div, divisors, frantic=False):
    # frantic = True: Test aggressively by testing cycle for factors p of d also 
    # with offsets from 1 to p-1.
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
        offsets = [0]
        if frantic:
            offsets = [x for x in range(p)]
        for offset in offsets:
            common_coefficient = None
            for j in range(0, div, p):
                i = j + offset
                # skip 0 so coeffcient for zeta**0 can be negative.
                if i == 0:
                    continue
                index = div - i - 1
                if common_coefficient is None:
                    common_coefficient = new_coeffs[index]
                elif common_coefficient > new_coeffs[index]:
                    common_coefficient = new_coeffs[index]
            if common_coefficient is not None and common_coefficient > 0:
                    for j in range(0, div, p):
                        i = j + offset
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
                product = product * (self.M[i, j] + 1)
            sp.expand(product)
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
    r = generate_solution_vectors_sp(d, remainder)
    v = M * r[remainder-1]
    return v


def calculate_product_matrix(d):
    v = sp.zeros(d, 1)
    
    # First number is always 2**d
    v[0,0] = 2**d

    # If d is power of two, all other numbers are 0
    if(power_of_two(d)):
        return v

    # If p is prime, all other numbers are 2
    divisors = sp.divisors(d)
    if len(divisors) == 2:
        for i in range(1, d):
            v[i,0] = 2
        return v

    trailing = sp.trailing(d)
    pow_of_2 = 2**trailing
    for i in range(1, d):
        # if d has 2 as factor more often than i, then this number is 0
        if trailing >= 1 and i % pow_of_2 != 0:
                v[i,0] = 0
                continue
        
        # All other cases
        v[i,0] = 2**sp.gcd(i, d)
    return v


def calcluate_sum_subsets_gcd(n, div):
    if n < 0 or div < 2:
        raise ValueError("n must be >= 0 and div must be >=2!")
    
    d = div
    v = sp.Matrix(div, 1, lambda row, col: 1 if row == 0 else 0)
    if n <= 0:
        return v

    remainder = n % d
    temp_n = n - remainder
    quotient = temp_n // d
    divisors = sp.divisors(d)

    v = calculate_product_matrix(d)
    for i in range(div):
        v[i,0] = v[i,0]**quotient
     
    # factors to generating functions
    # row: index for functions, col: index of coefficient
    zeta = sp.symbols("zeta")
    gf_factors = sp.Matrix(d, d, lambda row, col: zeta ** (((d-row)*col) % d))

    # Eliminate roots in polynomial myself since Sympy seems to be unable to eliminate roots in sum.
    v2 = v.T * gf_factors
    root = sp.root(1, d, 1)
    for i in range(div):
        v2[0,i] = shrink_poly(v2[0,i], zeta, d, divisors)
        # shrink_poly shall have already eliminated all roots.
        assert(len(v2[0,i].free_symbols)==0)
        v2[0,i] = v2[0,i].subs(zeta, root)
        assert(v2[0,i] % d == 0)
        v2[0,i] = v2[0,i] // d
    v = v2.T

    if remainder == 0:
        return v

    M = create_matrix_from_solution_vector_sp(v)
    r = generate_solution_vectors_sp(d, remainder)
    v = M * r[remainder-1]
    return v


def test_calculation(n, div):
    print(f"##### div: {div}, n: {n}")
    t_gcd = time.process_time()
    r1 = calcluate_sum_subsets_gcd(n, div)
    t_gcd = time.process_time() - t_gcd
    print(f"# calcluate_sum_subsets_gcd: {r1}")
    print(f"# time elapsed={t_gcd}")
    
    t_log = time.process_time()
    r2 = calcluate_sum_subsets_logarithmic(n, div)
    t_log = time.process_time() - t_log
    is_equal = r2.equals(r1)
    print(f"# Same to result from calcluate_sum_subsets_logarithmic? {is_equal}, time elapsed={t_log}")
    return is_equal


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

    root = sp.root(1, 6, 1)
    product = 1
    for i in range(1, 7):
        t = root**(i * 2)
        product *= (t + 1)
        print(f"i: {i}, root**i: {t + 1}")
    print(product.evalf())


def test_calculate_product_matrix(start, end_exclusive):
    for i in range(start, end_exclusive):
        M = CyclicSievingMatrix(i)
        t_PM = time.process_time()
        PM = M.return_product_matrix()
        t_PM = time.process_time() - t_PM
        t_CPM =  time.process_time()
        CPM = calculate_product_matrix(i)
        t_CPM =  time.process_time() - t_CPM
        is_equal = PM.equals(CPM)
        print(f"i: {i}, calculate_product_matrix: {CPM}, equals: {is_equal}, t_CPM: {t_CPM}, t_PM {t_PM}")
        assert(is_equal)

   
def main():
    init_printing(use_unicode=True)
    #test_roots_of_unity()
    #test_calculate_product_matrix(2, 66)
    for i in range(2, 50):
        assert(test_calculation(2000, i))


if __name__ == '__main__':
    main()