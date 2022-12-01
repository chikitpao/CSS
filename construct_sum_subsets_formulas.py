# Construct formulas for counting sum of subsets modulo div
# Problem related to youtube video "Olympiad level counting - How many subsets of {1,…,2000} have a sum divisible by 5"
#  by 3Blue1Brown (url: https://www.youtube.com/watch?v=bOXCLR3Wric)

import sympy as sp
#from sympy import Symbol
from sympy import symbols

from sympy.interactive.printing import init_printing

from calculate_sum_subsets import *


class KeySymbol:
    def __init__(self, symbol, index):
        self.sym = symbol
        self.members = [symbol]
        self.index = index

    def add_symbol(self, symbol):
        self.members.append(symbol)

    def __str__(self):
        return f"KeySymbol(members: {self.members}, index: {self.index})"

    def __repr__(self):
        return f"KeySymbol(members={self.members}, index={self.index})"


class SymbolTable:
    def __init__(self, div):
        self.div = div
        self.sym_c = symbols('c:' + str(div))
        self.lookup_v_s = {}  # value -> symbol
        self.key_symbols = {}  # key symbol -> KeySymbol objects
        self.sym_c2 = list(self.sym_c)
        self.r_div = generate_solution_vectors_sp(div)
        self.v_div = self.r_div[div-1]
        self.replace_symbols()

    def replace_symbols(self):
        for j in range(0, self.div):
            key = self.v_div[j, 0]
            if key in self.lookup_v_s:
                self.sym_c2[j] = self.lookup_v_s[key]
                self.key_symbols[self.sym_c2[j]].add_symbol(self.sym_c[j])
            else:
                self.lookup_v_s[key] = self.sym_c[j]
                self.key_symbols[self.sym_c[j]] = KeySymbol(self.sym_c[j], j)

    def get_v_div(self):
        return self.v_div

    def return_symbol_matrix(self):
        M = sp.zeros(self.div, self.div)
        for j in range(0, self.div):
            for k in range(0, self.div):
                M[(j + k) % self.div, j] = self.sym_c2[k]
        return M

    def calculate_result(self, v_div_diff):
        print("## SymbolTable.key_symbols: ",
              self.key_symbols, "len:", len(self.key_symbols))
        if len(self.key_symbols) == 1:
            symbol0 = list(self.key_symbols.keys())[0]
            n = self.div
            print(
                f"# Apply rule of \"power of two\" (#1). Formula: {symbol0} = 2**n // d")
            print(f"n = {n} => {symbol0} = {2**n // self.div}")
        elif len(self.key_symbols) == 2 or len(self.key_symbols) == 3:
            key_symbols_keys = list(self.key_symbols)
            diff = None
            index_of_second_symbol_in_diff = None
            for i in range(0, self.div):
                free_symbols = v_div_diff[i, 0].free_symbols
                if len(free_symbols) == 2:
                    symbol_list = list(free_symbols)
                    if symbol_list[0] == key_symbols_keys[0]:  # c_0 minus another c_i
                        diff = symbol_list[0] - symbol_list[1]
                        index_of_second_symbol_in_diff = self.key_symbols[symbol_list[1]].index
                    else:
                        diff = symbol_list[1] - symbol_list[0]
                        index_of_second_symbol_in_diff = self.key_symbols[symbol_list[0]].index
                    break

            if diff is None:
                print("Cannot find difference with 2 variables.")
                return

            print("Difference with 2 terms:",
                  free_symbols, " Difference:", diff)
            # Can be optimized for performance, but require more complex code
            M = sp.Matrix(self.div, 1, lambda row,
                          col: sp.cancel(v_div_diff[row, 0] / diff))
            print(f"v_div_diff / ({diff}): {M}")

            # Assume if there are 3 columns, there will be still columns cannot be canceled out yet and they will still have three variables
            row_index_with_3_symbols = None
            if len(self.key_symbols) == 3:
                for i in range(0, self.div):
                    if len(M[i, 0].free_symbols) == 3:
                        row_index_with_3_symbols = i
                        break
                if row_index_with_3_symbols is None:
                    print("Cannot find remaining row with 3 variables.")
                    return

            # print(v_div_diff[self.key_symbols[key_symbols_keys[0]].index])
            #print(len(M[self.key_symbols[key_symbols_keys[1]].index, 0].free_symbols))
            assert v_div_diff[self.key_symbols[key_symbols_keys[0]].index] == 0
            assert len(M[index_of_second_symbol_in_diff, 0].free_symbols) == 0

            # occurence of c_0
            i = len(self.key_symbols[key_symbols_keys[0]].members)
            j = self.v_div[self.key_symbols[key_symbols_keys[0]].index] - \
                self.v_div[index_of_second_symbol_in_diff]  # difference c_0 - c_i as n == d
            # factor of difference from previous step to current step
            k = M[index_of_second_symbol_in_diff, 0]
            #c_i_times = len(self.key_symbols[key_symbols_keys[1]].members)
            if len(self.key_symbols) == 2:
                print(f"# Apply rule of \"1 difference with 2 terms\" (#2).")
                print(
                    f"i: {i} # occurences of c_0 = mu_0; occurences of c_1 = mu_1")
                mu_0 = i
                mu_1 = self.div - mu_0
                print(f"j: {j} # difference c_0 - c_i as n = d")
                print(
                    f"k: {k} # factor of difference from previous step to current step")
                d = self.div
                n = d   # for example
                print(
                    f"# Formula for {key_symbols_keys[0]}: (2**n + (mu_1 * j * k**((n//d)-1)))// d")
                print(
                    f"# Formula for {key_symbols_keys[1]}:  (2**n - (i * j * k**((n//d)-1))// d")

                def formula2_0(n, d): return (
                    2**n + (mu_1 * j * k**((n//d)-1)))//d
                def formula2_1(n, d): return (
                    2**n - (i * j * k**((n//d)-1))) // d
                print(
                    f"# n = {n} => {key_symbols_keys[0]} = {formula2_0(n, d)}, {key_symbols_keys[1]} = {formula2_1(n, d)}")
            else:  # len(self.key_symbols) == 3:
                print(f"i: {i} # occurences of c_0 = mu_0")
                print(f"j: {j} # difference c_0 - c_i as n = d")
                print(
                    f"k: {k} # factor of difference from previous step to current step")
                mu_0 = i
                mu_1 = len(
                    self.key_symbols[self.sym_c[index_of_second_symbol_in_diff]].members)
                mu_2 = len(
                    self.key_symbols[self.sym_c[row_index_with_3_symbols]].members)
                print(
                    f"mu_1: {mu_1} # occurences of {self.sym_c[index_of_second_symbol_in_diff]}")
                print(
                    f"mu_2: {mu_2} # occurences of {self.sym_c[row_index_with_3_symbols]}")
                o = v_div_diff[row_index_with_3_symbols, 0].coeff(
                    key_symbols_keys[0])
                p = v_div_diff[row_index_with_3_symbols, 0].coeff(
                    self.sym_c[index_of_second_symbol_in_diff])
                q = v_div_diff[row_index_with_3_symbols, 0].coeff(
                    self.sym_c[row_index_with_3_symbols])
                diff_0 = self.v_div[0]-self.v_div[row_index_with_3_symbols]
                print(f"o (coefficient {key_symbols_keys[0]}): {o}")
                print(
                    f"p (coefficient {self.sym_c[index_of_second_symbol_in_diff]}): {p}")
                print(
                    f"q (coefficient {self.sym_c[row_index_with_3_symbols]}): {q}")
                print(
                    f"diff_0 ({key_symbols_keys[0]}-{self.sym_c[row_index_with_3_symbols]} as n = d): {diff_0}")
                # assert coeff_c0 == coeff_c_in_diff    # Condition can be false!
                assert (o + p + q) == 0
                assert o > 0 and p > 0 and q < 0
                #gcd = (-q).gcd(p).gcd(j).gcd(k)
                #assert(gcd > 1)
                #print(f"gcd(-q, p, j, k): {gcd}")
                # assert ((-q) % k) == 0 and (p % k) == 0 and (j % k) == 0 # Condition can be false!
                #print(f"m=1 =>  (-q) - p * j * k**(m-1): { (-q) - p * j * k**(0-1)}, c0-c1: {v_div[0]-v_div[1]}")
                # REMARK: type of e.g. t in "t = 2**2000" -> <class 'int'>
                # REMARK: type of e.g. p -> <class 'sympy.core.numbers.Integer'>
                # TODO CKP: try perfect_power, e.g. perfect_power(p)
                # assert(power_of_two(p)) # TypeError: unsupported operand type(s) for &: ...
                # assert(power_of_two(k)) # TypeError: unsupported operand type(s) for &: ...
                # TODO CKP: try integer_log, e.g. integer_log(-q, 2)
                assert (sp.log(p, 2) >= sp.log(k, 2))

                def formula3_part1(n, d): return j * k**((n//d)-1)
                # formula3_part2 = lambda n, d: self.v_div[row_index_with_3_symbols]

                def formula3_part2(n, d): return diff_0 * (-q)**((n//d)-1) - p * j * (k**((n//d)-2)) * (
                    ((2**(sp.log(-q, 2) - sp.log(k, 2)))**((n//d)-1)) - 1) // (2**(sp.log(-q, 2) - sp.log(k, 2)) - 1)
                print(f"# Apply rule of \"1 difference with 2 terms and a third term with this difference and a third variable\" (#3).")
                print("difference c0 - c of row_index_of_second_symbol_in_diff:")
                print("# f3_1(n, d): j * k**((n//d)-1)")
                print("difference c0 - c of row_index_with_3_symbols: ")
                #print(" m = 1: self.v_div[row_index_with_3_symbols]")
                print("# f3_2(n, d): m >= 1: j * (-q)**((n//d)-1) - p * j * (k**((n//self.div)-2)) * (s from 0 to (m-2) Sum of ((2**(log_2(p) - log_2(k))**s)")
                print("  => diff_0 * (-q)**((n//d)-1) - p * j * (k**((n//self.div)-2)) * ((2**(log_2(-q) - log_2(k))**s) - 1) // (2**(log_2(-q) - log_2(k)) - 1)")
                d = self.div
                # Formulas for coefficients
                print(
                    f"# Formula for {key_symbols_keys[0]}: (2**n + mu_1 *  f3_1(n, d) + mu_2 * f3_2(n, d))// d")
                print(
                    f"# Formula for {self.sym_c[index_of_second_symbol_in_diff]}: (2**n + (mu_1 - d) *  f3_1(n, d) + mu_2 * f3_2(n, d))// d")
                print(
                    f"# Formula for {self.sym_c[row_index_with_3_symbols]}: (2**n + mu_1 * f3_1(n, d) + (mu_2 - d) * f3_2(n, d))// d")

                def formula3_0(n, d): return (
                    2**n + mu_1 * formula3_part1(n, d) + mu_2 * formula3_part2(n, d)) // d

                def formula3_1(n, d): return (
                    2**n + (mu_1 - d) * formula3_part1(n, d) + mu_2 * formula3_part2(n, d)) // d

                def formula3_2(n, d): return (
                    2**n + mu_1 * formula3_part1(n, d) + (mu_2 - d) * formula3_part2(n, d)) // d
                n = d  # for example
                print(
                    f"# n = {n} => {key_symbols_keys[0]} = {formula3_0(n, d)}, {self.sym_c[index_of_second_symbol_in_diff]} = {formula3_1(n, d)}, {self.sym_c[row_index_with_3_symbols]} = {formula3_2(n, d)}")

                # Tests
                for ni in range(1, 6):
                    n = ni * d
                    print(
                        f"# n = {n} => c0-{self.sym_c[index_of_second_symbol_in_diff]}: {formula3_part1(n, d)}, c0-{self.sym_c[row_index_with_3_symbols]}: {formula3_part2(n, d)}")

                M = create_matrix_from_solution_vector_sp(self.v_div)
                v_new = self.v_div
                for var in range(1, 6):
                    print(
                        f"(matric calc.) m: {var}, c0-c1: {v_new[0]-v_new[1]}, c0-{self.sym_c[index_of_second_symbol_in_diff]}: {v_new[0]-v_new[index_of_second_symbol_in_diff]}")
                    v_new = M * self.v_div
                    M = create_matrix_from_solution_vector_sp(v_new)

            #    #<expr>.atoms(Symbol)
            #    temp, with_ci = v_div_diff[i, 0].as_independent((diff), as_Mul=True)
            #    print("temp, with_ci: ", temp, with_ci)
         # Filter rule I'm working on:
        else:  # i.e. len(self.key_symbols) > 3:
            raise NotImplementedError(
                f"## Ignore div={self.div} for now. Key symbols count {len(self.key_symbols)}.")


# test sympy symbols
def test_symbols(div):
    """Construct formulas for current problem"""
    symbolTable = SymbolTable(div)
    v_div = symbolTable.get_v_div()
    M = symbolTable.return_symbol_matrix()

    r = M * v_div
    v_div_diff = sp.Matrix(div, 1, lambda row, col: r[0, 0] - r[row, 0])
    print("##### div:", div)
    print("# v_div: ", v_div)
    print("# v_div_diff: ", v_div_diff)
    print("# r = M * v_div: ", r)
    # print("# M.det(): ", M.det())
    try:
        symbolTable.calculate_result(v_div_diff)
    except NotImplementedError as nie:
        print(str(nie))
    print("")


def test_symbols_2(div):  # print vectors only as symbols
    sym_c = symbols('c:' + str(div))
    v = sp.Matrix(div, 1, lambda row, col: sym_c[row])

    vs = []
    vs.append(v)
    for i in range(0, 4):
        M_i = create_matrix_from_solution_vector_sp(vs[i])
        v_new = M_i * v
        vs.append(v_new)

    print("##### div:", div)
    for i in range(0, 5):
        print("# i:", i, "vs[i]", vs[i])


def test_symbols_3(div):  # print vectors as multiples of initial symbols
    sym_c = symbols('c:' + str(div))
    v = sp.Matrix(div, 1, lambda row, col: sym_c[row])
    v_num = generate_solution_vectors_sp(div)[div-1]

    vs = []
    vs.append(v)
    for i in range(0, 4):
        M_i = create_matrix_from_solution_vector_sp(vs[i])
        vs_new = M_i * v_num
        vs.append(vs_new)

    print("##### div:", div)
    for i in range(0, 5):
        print("# i:", i, "vs[i]", vs[i])


def test_symbols_4():
    symbol_count_to_div = {}
    for div in range(2, 101):
        symbolTable = SymbolTable(div)

        count = len(symbolTable.key_symbols)
        if count in symbol_count_to_div:
            symbol_count_to_div[count].append(div)
        else:
            symbol_count_to_div[count] = [div]

    print(f"symbol_count_to_div: {symbol_count_to_div}")


def main():
    init_printing(use_unicode=False, wrap_line=False, num_columns=300)

    for i in range(2, 101):
        test_symbols(i)
    # test_symbols_2(5)
    # test_symbols_3(15)

    # test_symbols_4()


if __name__ == '__main__':
    main()
