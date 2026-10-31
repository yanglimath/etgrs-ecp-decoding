"""Matrix ECP and Beelen decoding through complete extended GRS codes."""

from itertools import chain

from sage.all import ZZ, codes, diagonal_matrix, matrix, vector
from sage.coding.decoder import DecodingError


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
