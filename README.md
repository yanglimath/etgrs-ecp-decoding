# Decoding (+)-ETGRS codes

A SageMath program comparing our ECP decoder in Algorithm 1 with Beelen-type decoders based on BW, Gao and EEA decoding methods.

The tests use the smallest admissible field and largest permitted length for each of the three constructions in Table 1, with all feasible dimensions:

- `(q, n+1, k) = (8, 7, 3), (81, 9, 3), (25, 21, k)` for `k = 3, ..., 9`.
- Every error weight `w` from `1` to `floor((N - 1 - k)/2)`.
- Random seed `20261017`, ten paired trials per case, one computational thread.

Run with SageMath 10.6:

```sh
sage -python etgrs_decoding.py
```

The program prints the 50 rows of mean decoding times in milliseconds, to four decimal places, and checks that all decoded codewords are correct. Timings vary across runs and machines.
