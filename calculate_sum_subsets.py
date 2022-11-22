# Count / calculate number of sums of subsets modulo d 
# Problem related to youtube video "Olympiad level counting - How many subsets of {1,…,2000} have a sum divisible by 5"
#  by 3Blue1Brown (url: https://www.youtube.com/watch?v=bOXCLR3Wric)
from itertools import chain, combinations
#import numpy as np
import sympy as sp
#from sympy import Symbol 
from sympy import symbols
import time

from sympy.interactive.printing import init_printing


## generate combinations and count
# Source: https://docs.python.org/3/library/itertools.html#itertools-recipes
def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))


def calculate_sum_mod(l, div):
    result = 0
    for e in l:
        result = (result + e) % div
    return result


def try_result(n, div):
    myList = [i for i in range(1, n + 1)]

    #print(myList)
    myPowerSetList = list(powerset(myList))
    mod = [0] * div
    for e in myPowerSetList:
        r = calculate_sum_mod(e, div)
        mod[r] = mod[r] + 1
    return mod


# ## calculate with matrix (numpy)
# def create_matrix_from_solution_vector_np(v):
#     dim = np.shape(v)[0] # 0 -> rows, 1 -> cols
#     M = np.empty((dim, dim), dtype=int)
#     for j in range(0, dim):
#         for k in range(0, dim):
#             M[(j + k) % dim][j] = v[k]
#     return M

            
# def generate_solution_vectors_np(div):
#     result = []
    
#     # Step 1: Create matrix for 0 elements (= Identity matrix)
#     M = np.eye(div, dtype=int)
    
#     # Step 2: Create vector for different values (needs new matrix from last value)
#     element_vector = np.zeros((div, 1), dtype=int)
#     element_vector[0][0] = 1
#     for i in range(1, div + 1):
#         r = i % div
#         element_vector[r][0] += 1
#         # N = M * element_vector
#         # v2 = np.multiply(M, element_vector)
#         v2 = np.matmul(M, element_vector)
#         result.append(v2)
#         element_vector[r][0] -= 1        
#         M = create_matrix_from_solution_vector_np(v2)
#     return result
# ##


## calculate with matrix (sympy)
def create_matrix_from_solution_vector_sp(v):
    dim = v.shape[0] # 0 -> rows, 1 -> cols
    M = sp.zeros(dim, dim)
    for j in range(0, dim):
        for k in range(0, dim):
            M[(j + k) % dim, j] = v[k]
    return M

            
def generate_solution_vectors_sp(div, maxVector=None):
    result = []
    
    vectorCount = div
    if(maxVector is not None and maxVector >= 1 and maxVector <= div):
        vectorCount = maxVector
    
    # Step 1: Create matrix for 0 elements (= Identity matrix)
    M = sp.eye(div)
    
    # Step 2: Create vector for different values (needs new matrix from last value)
    element_vector = sp.zeros(div, 1)
    element_vector[0, 0] = 1
    for i in range(1, vectorCount + 1):
        r = i % div
        element_vector[r, 0] += 1
        v2 = M * element_vector
        result.append(v2)
        element_vector[r, 0] -= 1        
        M = create_matrix_from_solution_vector_sp(v2)
    return result


def calcluate_sum_subsets_linear(n, div):
    if n < 0 or div < 2:
        raise ValueError("n must be >= 0 and div must be >=2!")

    v = sp.Matrix(div, 1, lambda row, col: 1 if row == 0 else 0)
    
    if n <= 0:
        return v
    
    remainder = n % div
    quotient = n // div

    r = generate_solution_vectors_sp(div)
    if quotient == 0:
        return r[remainder-1]

    M = sp.eye(div, div)

    # Normally we shall handle full chunks first and then remainder.
    # For besser program flow, order of handlings is swapped.
    if remainder != 0:
        v = M * r[remainder-1]
        M = create_matrix_from_solution_vector_sp(v)

    r_last = r[div-1]
    for i in range(1, quotient + 1):
        v = M * r_last
        M = create_matrix_from_solution_vector_sp(v)

    return v


def calcluate_sum_subsets_logarithmic(n, div):
    if n < 0 or div < 2:
        raise ValueError("n must be >= 0 and div must be >=2!")

    v = sp.Matrix(div, 1, lambda row, col: 1 if row == 0 else 0)
    
    if n <= 0:
        return v
    
    remainder = n % div
    quotient = n // div

    r = generate_solution_vectors_sp(div)
    if quotient == 0:
        return r[remainder-1]

    M = sp.eye(div, div)

    # Normally we shall handle full chunks first then remainder.
    # For besser program flow order of handlings is swapped.
    if remainder != 0:
        v = M * r[remainder-1]
        M = create_matrix_from_solution_vector_sp(v)

    temp_v = r[div-1]
    temp_quotient = quotient
    temp_M = sp.eye(div, div)
    is_lsb_set = False
    while temp_quotient > 0:
        is_lsb_set = temp_quotient & 1
        temp_v = temp_M * temp_v
        if is_lsb_set:
            M = M * temp_v
            M = create_matrix_from_solution_vector_sp(M)
        temp_M = create_matrix_from_solution_vector_sp(temp_v)
        temp_quotient = temp_quotient // 2

    return M[:, 0]


def calcluate_sum_subsets_constant(n, div): # Use formulas
    if n < 0 or div < 2:
        raise ValueError("n must be >= 0 and div must be >=2!")

    # REMARK CKP: We can't handle all divisors here! 
    v = sp.Matrix(div, 1, lambda row, col: 1 if row == 0 else 0)
    if n <= 0:
        return v

    remainder = n % div
    temp_n = n - remainder
    
    if div == 3:
        if temp_n > 0:
            c3_0 = lambda x: (2**x + 2 * 2**(x//3)) // 3
            c3_1 = lambda x: (2**x - 2**(x//3)) // 3
            c_0 = c3_0(temp_n)
            c_1 = c3_1(temp_n)
            v = sp.Matrix(div, 1, lambda row, col: c_0 if row == 0 else c_1)
        # else v remains the same
    elif div == 5:
        if temp_n > 0:
            c5_0 = lambda x: (2**x + 4 * 2**(x//5)) // 5
            c5_1 = lambda x: (2**x - 2**(x//5)) // 5
            c_0 = c5_0(temp_n)
            c_1 = c5_1(temp_n)
            v = sp.Matrix(div, 1, lambda row, col: c_0 if row == 0 else c_1)
        # else v remains the same
    elif div == 6:
        if temp_n > 0:
            c6_0 = lambda x: (2**x + 2 * 4**(x//6)) // 6
            c6_1 = lambda x: (2**x - 4**(x//6)) // 6
            c_0 = c6_0(temp_n)
            c_1 = c6_1(temp_n)
            v = sp.Matrix(div, 1, lambda row, col: c_0 if row % 3 == 0 else c_1)
        # else v remains the same
    elif div == 7:
        if temp_n > 0:
            c7_0 = lambda x: (2**x + 6 * 2**(x//5)) // 7
            c7_1 = lambda x: (2**x - 2**(x//5)) // 7
            c_0 = c7_0(temp_n)
            c_1 = c7_1(temp_n)
            v = sp.Matrix(div, 1, lambda row, col: c_0 if row == 0 else c_1)
        # else v remains the same
    elif div == 9:
        if temp_n > 0:
            c9_0 = lambda x: (2**x + 2 * 2**(x//9) + 2* 2**(x//9 - 1) * (2 * 4**(x//9) + 4)) // 9
            c9_1 = lambda x: (2**x + 2 * 2**(x//9) - 2**(x//9 - 1) * (2 * 4**(x//9) + 4)) // 9
            c9_3 = lambda x: (2**x - 7 * 2**(x//9) + 2 * 2**(x//9 - 1) * (2 * 4**(x//9) + 4)) // 9
            c_0 = c9_0(temp_n)
            c_1 = c9_1(temp_n)
            c_3 = c9_3(temp_n)
            for i in range(0, div):
                if i == 0:
                    v[i, 0] = c_0
                elif i % 3 == 0:
                    v[i, 0] = c_3
                else:
                    v[i, 0] = c_1
        # else v remains the same
    else:
        if div >= 2 and power_of_two(div):
            if temp_n > 0:
                c2n_0 = lambda x: 2**temp_n // div
                c_0 = c2n_0(temp_n)
                v = sp.Matrix(div, 1, lambda row, col: c_0)
            # else v remains the same
        else:
            raise ValueError(f"div={div} is not supported by calcluate_sum_subsets_constant")

    if remainder == 0:
        return v

    M = create_matrix_from_solution_vector_sp(v)
    r = generate_solution_vectors_sp(div, remainder)
    v = M * r[remainder-1]
    return v

 # Source for "power of two": https://stackoverflow.com/questions/57025836/how-to-check-if-a-given-number-is-a-power-of-two
def power_of_two(x):
    return (x & (x-1) == 0) and x > 0
##


## test try_result
def test_try_result():
    for div in range(2, 21):
        for i in range (1, div + 1):
            print(f"i: {i}, div: {div}, try_result: {try_result(i, div)}")


# ## test numpy functions
# def test_numpy():
#     c6_0 = lambda x: (2**x + 2 * 4**(x//6)) // 6
#     c6_1 = lambda x: (2**x - 4**(x//6)) // 6
#     div = 6
#     r6 = generate_solution_vectors_np(div)
#     print(r6)
#     v6 = r6[div-1]
#     n = div
#     print(f"6: v6: {v6} v6[0]-v6[1]: {v6[0]-v6[1]} n: {n} c6_0(n): {c6_0(n)} c6_1(n) {c6_1(n)}")
#     for i in range(2, 6):
#         M = create_matrix_from_solution_vector_np(v6)
#         v6 = np.matmul(M, r6[div-1])
#         n = div * i
#         print(f"6: i {i} finished: v6: {v6} v6[0][0]-v6[1][0]: {v6[0][0]-v6[1][0]} n: {n} c6_0(n): {c6_0(n)} c6_1(n): {c6_1(n)}")    


## test sympy functions
def test_sympy():
    c6_0 = lambda x: (2**x + 2 * 4**(x//6)) // 6
    c6_1 = lambda x: (2**x - 4**(x//6)) // 6
    div = 6
    r6 = generate_solution_vectors_sp(div)
    print(r6)
    v6 = r6[div-1]
    n = div
    print(f"6: v6: {v6} v6[0]-v6[1]: {v6[0]-v6[1]} n: {n} c6_0(n): {c6_0(n)} c6_1(n) {c6_1(n)}")
    for i in range(2, 12):
        M = create_matrix_from_solution_vector_sp(v6)
        v6 = M * v6  # was r6[div-1]
        n = 2 * n
        print(f"6: i {i} finished: v6: {v6} v6[0,0]-v6[1,0]: {v6[0,0]-v6[1,0]} n: {n} c6_0(n): {c6_0(n)} c6_1(n): {c6_1(n)}")    


def test_power_of_two():
    t = 2**2000
    print(t, power_of_two(t))
    t += 1
    print(t, power_of_two(t))

def test(n, div):
    print(f"### n={n} div={div}")
    
    t = time.process_time()
    r_log = calcluate_sum_subsets_logarithmic(n, div)
    time_log = time.process_time() - t

    t = time.process_time()
    r_lin = calcluate_sum_subsets_linear(n, div)
    time_lin = time.process_time() - t
    
    r = r_log - r_lin
    print(f"result logarithm={r_lin}, elapsed={time_log} difference with linear={r}, elapsed={time_lin}")

    try:
        t = time.process_time()
        r_const = calcluate_sum_subsets_constant(n, div)
        time_con = time.process_time() - t
        r = r_log - r_const
        print(f"result logarithm={r_lin}, elapsed={time_log} difference with constant={r}, elapsed={time_con}")
    except:
        pass


def main():
    init_printing(use_unicode=False, wrap_line=False, num_columns=300)
    #test_numpy()
    #test_sympy()

    #div = 9
    #for i in range (0, 20):
    #    test(i, div)
    
    test(2003, 9)


if __name__ == '__main__':
    main()
