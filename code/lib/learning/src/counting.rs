//! Various utility functions for counting and combinatorics.

/// A function to count the number of ways to choose `r` things from `n` potential options.
///
/// This is taken directly from
/// [here](https://stackoverflow.com/questions/65561566/number-of-combinations-permutations).
pub fn count_combinations(n: usize, r: usize) -> usize {
    if r > n {
        0
    } else {
        (1..=r.min(n - r)).fold(1, |acc, val| acc * (n - val + 1) / val)
    }
}
