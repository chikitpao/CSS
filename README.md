CSS
===
Count / calculate number of sums of subsets of n numbers modulo d.

Problem related to youtube video "Olympiad level counting - How many subsets of {1,…,2000} have a sum divisible by 5" by 3Blue1Brown (url: https://www.youtube.com/watch?v=bOXCLR3Wric).


Installation
------------
You need to install Python and Sympy to run my python scripts. 

File list
------------
- **doc/Solutions of counting problem.odt**: Mathematical description of solutions to counting problem.
- **doc/Solutions of counting problem.pdf**: Mathematical description of solutions to counting problem, in PDF format.
- **calculate_sum_subsets.py**: Calculate number of sums of subsets of n numbers modulo d.
- **calculate_sum_subsets_cs.py**: Calculate number of sums of subsets of n numbers modulo d via cyclic sieving on a generating function.
- **construct.txt**: Output of script construct_sum_subsets_formulas.py.
- **construct_sum_subsets_formulas.py**: Construct formulas for different divisors (not working for several values of d).
- **construct_sum_subsets_formulas2.py**: Construct formulas for different divisors (for d from 2 to 100).
- **construct2.txt**: Output of script construct_sum_subsets_formulas2.py.
- **sum_subsets_formula.py**: Contains class with formulas. Used by function **calcluate_sum_subsets_constant** in  **calculate_sum_subsets.py**.

Usage
------------
- calculate_sum_subsets.py
    * **calcluate_sum_subsets_linear(n, div)**: Use matrix multiplications to calculate number of sums of subsets of n numbers modulo div. Return a dx1-matrix with count of sums congruent to (0, 1, ...) modulo d.
    * **calcluate_sum_subsets_logarithmic(n, div)**: Like calcluate_sum_subsets_linear, but only has logarithmic complexity instead of linear. Return a dx1-matrix with count of sums congruent to (0, 1, ...) modulo d.
    * **calcluate_sum_subsets_constant(n, div, sum_subset_formulas: SumSubsetsFormulas)**: Use formulas to calculate number of sums of subsets of n numbers modulo div. Return a dx1-matrix with count of sums congruent to (0, 1, ...) modulo d. Needs an object of class SumSubsetsFormulas. Throws ValueError if d is not supported (supports d from 2 to 100).
- calcluate_sum_subsets_cs.py
    * **calcluate_sum_subsets_cs(n, div, debug=False)**: Use cyclic sieving on a generating function to calculate number of sums of subsets of n numbers modulo div. Return a dx1-matrix with count of sums congruent to (0, 1, ...) modulo d. With parameter **debug** is set to True, additional information is shown to the output which is helpful if you want to create formulas for calculation :-).
    * **calcluate_sum_subsets_gcd(n, div)**: Like calcluate_sum_subsets_cs, but is much faster by replacing the most polynomial multiplications with evaluation of behavior regarding roots of unity.

**CAUTION:** calculation will take longer with larger d. It might be caused by calculation with large matrices. For function **calcluate_sum_subsets_cs**, I have to elminiate roots in polynoimal myself since SymPy seems to be unable to eliminate roots in sum. These are calculation times of function **calcluate_sum_subsets_cs** with n = 2000:
* d = 5: time = 0.063 s
* d = 10: time = 0.12 s
* d = 15: time = 0.4 s
* d = 20: time = 1.0 s
* d = 25: time = 2.1 s
* d = 50: time = 64 s

License & Copyright
-------------------
This open source release is licensed under the CC0 license. All trademarks are the property of their respective owners.
