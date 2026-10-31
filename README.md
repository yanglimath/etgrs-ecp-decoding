# Decoding (+)-ETGRS codes

SageMath programs for *Explicit Error-Correcting Pairs and Some Structural Properties of (+)-ETGRS Codes*.

- `etgrs_decoding.py`: the matrix ECP decoder and Beelen-type decoding with Sage's BW, Gao, or EEA decoder.
- `benchmark.py`: the nine code instances and the timing comparison in the paper.

## Run

Tested with SageMath 10.6 and Python 3.12.5.

```sh
sage -python benchmark.py --output results
```

By default, seed `20261017` gives ten fresh message/error pairs for every code and positive error weight, shared by all four decoders. Use `--seeds 0 1 42 1729 20261017` to include multiple seeds; every requested seed is reported separately and included in the pooled results. The script records individual timings, means, medians, and environment information. It uses one computational thread and times complete decoding calls; fixed preprocessing and input generation are excluded. The script records the seed derivation and field encodings. The reported seed-20261017 timings were extracted from the original five-seed run.

## Use the decoders

`ECPDecoder(S, v, k, eta)` and `BeelenDecoder(S, v, k, eta, method="Gao")` expose `decode(y)`. They return a full codeword or raise Sage's `DecodingError`. The method can be `"BW"`, `"Gao"`, or `"EEA"`; the radius is `floor((len(S)-k)/2)`.

The Beelen loop guesses the coefficient, subtracts its twist contribution, and decodes the full EGRS residual while retaining the received last coordinate. It then checks the coefficient and full distance. The received last symbol is tried first.

The EGRS step uses an exact full-length GRS representation. The three methods call Sage's `BerlekampWelch`, `Gao`, and `KeyEquationSyndrome` decoders. EEA denotes the extended Euclidean algorithm used by the last decoder. No coordinate is deleted. The representation requires `len(S) < |F|`; EEA additionally needs `len(S)+1 < |F|` to choose nonzero evaluation points. All nine instances satisfy this condition. The formula is in `ExtendedGRSDecoder`.

The guessing framework follows Beelen, Puchinger, and Rosenkilde, “Twisted Reed-Solomon codes,” ISIT 2017, pp. 336–340, Section VII.
