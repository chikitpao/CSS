# Construct formulars for counting sum of subsets modulo div via cyclic sieving 
# Problem related to youtube video "Olympiad level counting - How many subsets of {1,…,2000} have a sum divisible by 5"
#  by 3Blue1Brown (url: https://www.youtube.com/watch?v=bOXCLR3Wric)

import sympy as sp

from sympy.interactive.printing import init_printing
from sympy.abc import x

from calculate_sum_subsets import *

def shrink_poly(expr, sym, div):
    new_expr = sp.poly(sp.expand(expr), sym)
    
    # Replace power with (power % div)
    old_coeffs = new_expr.all_coeffs()
    new_coeffs = [0] * div
    old_coeffs_len = len(old_coeffs)
    for i in range(old_coeffs_len):
        new_coeffs[div - 1 - (i % div)] += old_coeffs[old_coeffs_len - 1 - i]

    # Cancel out whole sets of "roots of unity".
    min_coeff = min(new_coeffs)
    if min_coeff > 0:
        #print("min_coeff:", min_coeff)
        for i in range(len(new_coeffs)):
            new_coeffs[i] -= min_coeff
    new_poly = sp.Poly(new_coeffs, sym)
    return new_poly

class CyclicSievingMatrix:
    def __init__(self, d):
        self.d = d
        self.zeta = sp.symbols("zeta")
        # M: Matrix of zetas in as functions of sum 
        # row: index for functions, col: index of coefficient
        self.M = sp.Matrix(d, d, lambda row, col: self.zeta**((row * (col + 1)) % d))
        #print("# Matrix:", self.M)
        self.list_of_roots = [sp.root(1, d, k) for k in range(d)]
        #root = self.list_of_roots[1]
 
        # factors to functions
        # row: index for functions, col: index of coefficient
        self.factors = sp.Matrix(d, d, lambda row, col: self.zeta ** (((d-row)*col) % d))

    def print_sum_of_columns(self):
        # row: index for functions, col: index of coefficient
        v = self.M.T * self.factors
        print("# v: ", v)
        for i in range(0, self.d):
            for j in range(0, self.d):
                v[i, j] = v[i, j].subs(self.zeta,  self.list_of_roots[1]).evalf()
        for i in range(0, self.d):
            print(f"sum of columns after subs and evalf for i={i}: {v[:, i]}")
    
    def return_product_matrix(self):
        # row: index for functions
        v = sp.zeros(self.d, 1)        
        for i in range(0, self.d):
            product = 1
            for j in range(0, self.d):
                product = product * (self.M[i, j] + 1)
            new_poly = shrink_poly(product, self.zeta, self.d)
            #print("product", product, product.subs(self.zeta,  self.list_of_roots[1]).evalf())
            #print("new_poly", new_poly, new_poly.subs(self.zeta,  self.list_of_roots[1]))
            #v[i] = self.replace_poly_with_roots(new_poly)
            v[i] = new_poly.subs(self.zeta,  self.list_of_roots[1]).evalf()
        #print(f"product of rows: {v}")
        return v
        
    def replace_poly_with_roots(self, polynomial):
        coeffs = polynomial.all_coeffs()
        for i in range(1, len(self.list_of_roots)):
            polynomial = polynomial.subs(self.zeta**i, self.list_of_roots[i])
            print(f"polynomial after i={i}: {polynomial}")
        return polynomial


def calcluate_sum_subsets_cs(n, div):
    if n < 0 or div < 2:
        raise ValueError("n must be >= 0 and div must be >=2!")
    
    d = div
    v = sp.Matrix(div, 1, lambda row, col: 1 if row == 0 else 0)
    if n <= 0:
        return v

    remainder = n % d
    temp_n = n - remainder
    quotient = temp_n // d

    M = CyclicSievingMatrix(div)
    v = M.return_product_matrix()
    for i in range(div):
        v[i,0] = v[i,0]**quotient / d
    v2 = v.T * M.factors
    for i in range(div):
        v2[0, i] = v2[0, i].subs(M.zeta, M.list_of_roots[1]).evalf()
    #for i in range(div):
    #    print(f"#{i}: {v2[:,i]}")
    v = v2.T

    if remainder == 0:
        return v

    M = create_matrix_from_solution_vector_sp(v)
    r = generate_solution_vectors_sp(d)
    v = M * r[remainder-1]
    return v

def test():
    div = 15
    n = 2000
    print(f"##### div: {div}, n: {n}")
    r = calcluate_sum_subsets_cs(n, div)
    print(f"calcluate_sum_subsets_cs: {r}")
    r = calcluate_sum_subsets_logarithmic(n, div)
    print(f"calcluate_sum_subsets_logarithmic: {r}")

def main():
    init_printing(use_unicode=True)
    test()


if __name__ == '__main__':
    main()