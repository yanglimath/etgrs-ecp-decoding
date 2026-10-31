"""Compare ECP and Beelen-type decoding on the nine (+)-ETGRS code instances."""

import hashlib
import os
import random
import statistics
from itertools import chain
from time import perf_counter_ns

for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "SAGE_NUM_THREADS"):
    os.environ[name] = "1"

from sage.all import GF, PolynomialRing, ZZ, codes, diagonal_matrix, matrix, vector
from sage.coding.decoder import DecodingError

SEED = 20261017
TRIALS = 10
METHODS = ("ECP", "BW", "Gao", "EEA")

#(+)-ETGRS codes obtained from Table 1#

CODE_IDS = (103, 203, 303, 304, 305, 306, 307, 308, 309)
FIELDS = {
    8: (2, 3, [1, 1, 0, 1]),
    81: (3, 4, [2, 0, 0, 2, 1]),
    25: (5, 2, [2, 4, 1]),
}
POINTS = {
    1: [1, 2, 4, 3, 6, 7],
    2: [73, 74, 36, 2, 38, 37, 72, 1],
    3: [1, 5, 8, 12, 21, 2, 10, 11, 24, 17, 4, 20, 22, 18, 9, 3, 15, 19, 6, 13],
}
MULTIPLIERS = {
    103: [2, 3, 2, 1, 1, 3],
    203: [70, 36, 62, 37, 39, 73, 11, 1],
    303: [23, 10, 9, 12, 17, 1, 12, 22, 5, 11, 14, 5, 17, 6, 21, 2, 6, 19, 15, 8],
    304: [23, 2, 7, 1, 2, 3, 14, 2, 7, 14, 16, 4, 14, 2, 4, 4, 23, 1, 14, 23],
    305: [23, 21, 20, 19, 5, 4, 22, 6, 21, 24, 14, 13, 10, 22, 15, 3, 11, 12, 13, 12],
    306: [23, 12, 17, 18, 19, 2, 10, 22, 15, 13, 16, 24, 9, 6, 8, 1, 20, 11, 5, 21],
    307: [23, 23, 16, 4, 7, 1, 1, 2, 16, 4, 14, 14, 23, 2, 16, 2, 3, 4, 23, 2],
    308: [23, 8, 10, 11, 24, 3, 9, 6, 13, 5, 16, 11, 20, 22, 18, 4, 13, 18, 21, 10],
    309: [23, 5, 21, 12, 21, 4, 24, 22, 20, 22, 14, 15, 13, 6, 13, 3, 12, 19, 10, 11],
}

def _parameters(S, v, k):
    S = tuple(S)
    if not S:
        raise ValueError("The evaluation set must be nonempty")
    try:
        F = S[0].parent()
    except AttributeError as error:
        raise ValueError("The evaluation points must be finite-field elements") from error
    if not F.is_field() or not F.is_finite():
        raise ValueError("The evaluation points must lie in a finite field")
    try:
        k = int(ZZ(k))
    except (TypeError, ValueError) as error:
        raise ValueError("The dimension must be an integer") from error
    S, v = tuple(F(a) for a in S), tuple(F(b) for b in v)
    if len(set(S)) != len(S) or len(v) != len(S) or any(not b for b in v):
        raise ValueError("Use distinct evaluation points and nonzero multipliers")
    if not 1 <= k <= len(S):
        raise ValueError("The dimension must satisfy 1 <= k <= len(S)")
    return F, S, v, k


def _word(F, y, length):
    y = vector(F, y)
    if len(y) != length:
        raise ValueError("The received vector has the wrong length")
    return y


class ExtendedGRSDecoder:
    """Decode all n+1 EGRS coordinates by an exact projective representation.

    For beta outside S, z_i=(a_i-beta)^(-1), z_infinity=0 and
    V_i=v_i*(a_i-beta)^(k-1), V_infinity=1. The polynomial map
    h(z)=z^(k-1)*g(beta+1/z) preserves every codeword coordinate.
    A uniform translation z -> z-c leaves the codeword set and multipliers
    unchanged. When N=n+1 < |F|, c can avoid all z points, giving nonzero
    evaluation points for Sage's KeyEquationSyndrome decoder. The complete
    received vector is passed unchanged to the built-in decoder.
    The EGRS representation needs n < |F|; EEA additionally needs N < |F|.
    """

    def __init__(self, S, v, k, method="Gao"):
        self.F, self.S, self.v, self.k = _parameters(S, v, k)
        F, S, v, k = self.F, self.S, self.v, self.k
        names = {"BW": "BerlekampWelch", "Gao": "Gao", "EEA": "KeyEquationSyndrome"}
        if method not in names:
            raise ValueError("Choose method='BW', 'Gao', or 'EEA'")
        points = set(S)
        beta = next((a for a in F if a not in points), None)
        if beta is None:
            raise ValueError("This EGRS representation requires a field element outside S")
        self.beta = beta
        unshifted = tuple(1/(a-beta) for a in S) + (F.zero(),)
        point_set = set(unshifted)
        shift = next((a for a in F if a not in point_set), None)
        if shift is None and method == "EEA":
            raise ValueError("EEA requires N < |F| to obtain nonzero evaluation points")
        self.shift = F.zero() if shift is None else shift
        self.evaluation_points = tuple(a-self.shift for a in unshifted)
        self.multipliers = tuple(b*(a-beta)**(k-1)
                                 for a, b in zip(S, v)) + (F.one(),)
        self.length = len(S)+1
        self.t = (self.length-k)//2
        self.method = method
        self.code = codes.GeneralizedReedSolomonCode(
            list(self.evaluation_points), k, list(self.multipliers))
        self.H = self.code.parity_check_matrix()
        self.sage_decoder_name = names[method]
        self._decoder = self.code.decoder(self.sage_decoder_name)

    def decode(self, y):
        """Return a complete EGRS codeword within the guaranteed unique radius."""
        y = _word(self.F, y, self.length)
        try:
            decoded = self._decoder.decode_to_code(y)
        except ZeroDivisionError as error:
            raise DecodingError("The EGRS decoder failed") from error
        except ValueError as error:
            # BW/Gao can reach their encoder with an over-degree polynomial
            # on an uncorrectable input. Other ValueErrors remain visible.
            expected = "The polynomial to encode must have degree at most %s" % (self.k-1)
            if str(error) != expected:
                raise
            raise DecodingError("The EGRS decoder returned an over-degree polynomial") from error
        if not (self.H*decoded).is_zero() or (y-decoded).hamming_weight() > self.t:
            raise DecodingError("No valid EGRS candidate was returned within the radius")
        return decoded


class ECPDecoder:
    """The paper's matrix ECP decoder, with fixed diagonal-matrix assembly."""

    def __init__(self, S, v, k, eta):
        self.F, self.S, self.v, self.k = _parameters(S, v, k)
        F, S, v, k = self.F, self.S, self.v, self.k
        n = len(S)
        eta = F(eta)
        if not eta or not 3 <= k <= n-2:
            raise ValueError("The ECP construction requires eta != 0 and 3 <= k <= n-2")
        self.n, self.length, self.t = n, n+1, (n-k)//2
        weights = []
        for i, a in enumerate(S):
            denominator = F.one()
            for j, b in enumerate(S):
                if i != j:
                    denominator *= a-b
            weights.append(1/denominator)
        powers = [[a**h for a in S] for h in range(n-k+1)]
        self.H = matrix(F, [[weights[j]/v[j]*powers[h][j] for j in range(n)]
                             + [F.zero()] for h in range(n-k+1)])
        self.H[n-k-1, n] = -eta
        self.H[n-k, n] = -1-eta*sum(S, F.zero())
        self.GA = matrix(F, [powers[h]+[F.zero()] for h in range(self.t+1)])
        self.GB = matrix(F, [[weights[j]/v[j]*powers[h][j] for j in range(n)]
                              + [F.zero()] for h in range(self.t)])
        if (n-k) % 2 == 0:
            self.GA[self.t, n] = 1
            self.GB[self.t-1, n] = -eta
        else:
            self.GB[self.t-1, n] = 1
        self._GAT = self.GA.transpose()

    def decode(self, y):
        """Return a complete (+)-ETGRS codeword within t, or report failure."""
        y = _word(self.F, y, self.length)
        syndrome = self.H*y
        if syndrome.is_zero():
            return y
        M = self.GB*diagonal_matrix(y)*self._GAT
        kernel = M.right_kernel()
        if kernel.dimension() == 0:
            raise DecodingError("The ECP locating system has no nonzero solution")
        locator = kernel.basis()[0]*self.GA
        positions = [i for i, value in enumerate(locator) if not value]
        HZ = self.H.matrix_from_columns(positions)
        if not positions or HZ.rank() != len(positions):
            raise DecodingError("The ECP error values are not uniquely determined")
        try:
            values = HZ.solve_right(syndrome)
        except ValueError as error:
            raise DecodingError("The ECP error-value system is inconsistent") from error
        if values.hamming_weight() > self.t:
            raise DecodingError("The ECP candidate exceeds the decoding radius")
        error = vector(self.F, self.length)
        for i, value in zip(positions, values):
            error[i] = value
        return y-error


class BeelenDecoder:
    """Guess the twist coefficient and decode the complete EGRS residual."""

    def __init__(self, S, v, k, eta, method="Gao"):
        self.F, self.S, self.v, self.k = _parameters(S, v, k)
        F, S, v, k = self.F, self.S, self.v, self.k
        eta = F(eta)
        if not eta:
            raise ValueError("The twist parameter must be nonzero")
        self.n, self.length = len(S), len(S)+1
        self.t = (len(S)-k)//2
        self.method = method
        self._twist = vector(F, [eta*b*a**k for a, b in zip(S, v)] + [F.zero()])
        self.egrs = ExtendedGRSDecoder(S, v, k, method)

    def decode(self, y):
        """Try the received last coordinate first and stop at a verified candidate."""
        y = _word(self.F, y, self.length)
        first = y[-1]
        for b in chain((first,), (a for a in self.F if a != first)):
            twist = b*self._twist
            residual = y-twist  # Its final coordinate is still the received y[-1].
            try:
                candidate = self.egrs.decode(residual)
            except DecodingError:
                continue
            if candidate[-1] != b:
                continue
            corrected = candidate+twist
            if (y-corrected).hamming_weight() <= self.t:
                return corrected
        raise DecodingError("No (+)-ETGRS candidate was found within the decoding radius")

def make_case(code_id):
    """Reconstruct the original field representation and construction parameters."""
    family, k = divmod(code_id, 100)
    q, eta_number = {1: (8, 6), 2: (81, 2), 3: (25, 3)}[family]
    p, degree, modulus = FIELDS[q]
    R = PolynomialRing(GF(p), "z")
    F = GF(q, name="a", modulus=R(modulus))
    elements = tuple(sum((F((j // p**i) % p) * F.gen()**i
                          for i in range(degree)), F.zero()) for j in range(q))
    S = [elements[j] for j in POINTS[family]]
    v = [elements[j] for j in MULTIPLIERS[code_id]]
    eta = elements[eta_number]
    n = len(S)
    G = matrix(F, [[v[j] * S[j]**i for j in range(n)] + [F.zero()]
                   for i in range(k - 1)] +
                  [[v[j] * (S[j]**(k - 1) + eta * S[j]**k) for j in range(n)]
                   + [F.one()]])
    assert G.rank() == k and (G * G.transpose()).is_zero()
    return dict(code_id=code_id, family=family, q=q, n=n, k=k, N=n + 1,
                t=(n - k) // 2, F=F, S=S, v=v, eta=eta, G=G, elements=elements)

def received_word(case, weight, codeword_seed, error_seed, infinity_error=None):
    """Draw independent uniform message coefficients and nonzero error values."""
    F, elements = case["F"], case["elements"]
    message_rng, error_rng = random.Random(codeword_seed), random.Random(error_seed)
    message = vector(F, [elements[message_rng.randrange(case["q"])]
                         for _ in range(case["k"])])
    sent = message * case["G"]
    if infinity_error is None:
        positions = error_rng.sample(range(case["N"]), weight)
    elif infinity_error:
        positions = error_rng.sample(range(case["n"]), weight - 1) + [case["n"]]
    else:
        positions = error_rng.sample(range(case["n"]), weight)
    error = vector(F, case["N"])
    for position in positions:
        error[position] = elements[error_rng.randrange(1, case["q"])]
    received = sent + error
    received.set_immutable()
    return sent, received, int(case["n"] in positions)

def sample_seeds(code_id, weight, trial, phase="trial"):
    """Derive the message and error seeds from the fixed master seed."""
    prefix = f"{phase}:{SEED}:{code_id}:{weight}:{trial}"
    return tuple(int.from_bytes(hashlib.sha256(f"{prefix}:{kind}".encode()).digest(),
                                "big") for kind in ("message", "error"))

def run():
    print("  q   N   k   t   w      ECP(ms)     B-BW(ms)    B-Gao(ms)    B-EEA(ms)")
    verified = 0
    inputs = 0
    combinations = 0
    for code_id in CODE_IDS:
        case = make_case(code_id)
        S, v, k, eta = (case[key] for key in ("S", "v", "k", "eta"))
        decoders = {"ECP": ECPDecoder(S, v, k, eta)}
        decoders.update({name: BeelenDecoder(S, v, k, eta, method=name)
                         for name in METHODS[1:]})
        functions = {name: decoder.decode for name, decoder in decoders.items()}
        for weight in range(1, case["t"] + 1):
            # Warm up once without and once with an error at infinity.
            for j in range(2):
                message_seed, error_seed = sample_seeds(code_id, weight, j, "warmup")
                sent, received, _ = received_word(case, weight, message_seed,
                                                  error_seed, infinity_error=bool(j))
                for function in functions.values():
                    if function(received) != sent:
                        raise AssertionError((code_id, weight, "warmup"))
            times = {name: [] for name in METHODS}
            for trial in range(1, TRIALS + 1):
                message_seed, error_seed = sample_seeds(code_id, weight, trial)
                sent, received, _ = received_word(case, weight, message_seed, error_seed)
                shift = (trial - 1 + code_id + weight) % len(METHODS)
                order = METHODS[shift:] + METHODS[:shift]
                for name in order:
                    function = functions[name]
                    start = perf_counter_ns()
                    decoded = function(received)
                    elapsed = perf_counter_ns() - start
                    if decoded != sent or elapsed <= 0:
                        raise AssertionError((code_id, weight, trial, name))
                    times[name].append(elapsed)
                    verified += 1
                inputs += 1
            means = [statistics.fmean(times[name]) / 10**6 for name in METHODS]
            print(f"{case['q']:3d} {case['N']:3d} {k:3d} {case['t']:3d} {weight:3d}"
                  + "".join(f"{value:13.4f}" for value in means), flush=True)
            combinations += 1
    assert combinations == 50 and inputs == 500 and verified == 2000
    print(f"\nAll {verified} decoding calls returned the transmitted codeword.")


if __name__ == "__main__":
    run()

