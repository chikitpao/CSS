CSS
===
Count / calculate number of sums of subsets of n numbers modulo d.

Problem related to youtube video "Olympiad level counting - How many subsets of {1,…,2000} have a sum divisible by 5" by 3Blue1Brown (url: https://www.youtube.com/watch?v=bOXCLR3Wric).


Installation
------------
You need to install Python and Sympy to run my python scripts. 

File list
------------
- **calculate_sum_subsets.py**: Calculate number of sums of subsets of n numbers modulo d.
- **construct_sum_subsets_formulas.py**: Construct formulas for different divisors.
- **calculate_sum_subsets_cs.py**: Calculate number of sums of subsets of n numbers modulo d via cyclic sieving.
- **construct.txt**: Output of script construct_sum_subsets_formulas.py.

Usage
------------
- calculate_sum_subsets.py
    * **calcluate_sum_subsets_linear(n, div)**: Use matrix multiplications to calculate number of sums of subsets of n numbers modulo div. Return a (d, 1) matrix with count of sums congruent to (0, 1, ...) modulo d.
    * **calcluate_sum_subsets_logarithmic(n, div)**: Like calcluate_sum_subsets_linear, but only has logarithmic complexity instead of linear. Return a (d, 1) matrix with count of sums congruent to (0, 1, ...) modulo d.
    * **calcluate_sum_subsets_constant(n, div)**: Use formulas to calculate number of sums of subsets of n numbers modulo div. Return a (d, 1) matrix with count of sums congruent to (0, 1, ...) modulo d. Throws ValueError if d is not supported.
- calcluate_sum_subsets_cs.py
    * **calcluate_sum_subsets_cs(n, div)**: Use cyclic sieving to calculate number of sums of subsets of n numbers modulo div. Return a (d, 1) matrix with count of sums congruent to (0, 1, ...) modulo d. *It's not fully functional since floating point values are returned.*


License & Copyright
-------------------
This open source release is licensed under the CC0 license. All trademarks are the property of their respective owners.
