"""Compare ECP and three Beelen-type decoders on the nine paper instances."""

import argparse
import csv
import hashlib
import json
import os
import platform
import random
import statistics
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "SAGE_NUM_THREADS"):
    os.environ[name] = "1"

from sage.all import GF, PolynomialRing, matrix, vector
from sage.env import SAGE_VERSION
from etgrs_decoding import BeelenDecoder, ECPDecoder

DEFAULT_SEEDS = (20261017,)
METHODS = ("ECP", "BW", "Gao", "EEA")
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


def system_info():
    info = dict(sage=SAGE_VERSION, python=sys.version, platform=platform.platform(),
                machine=platform.machine(), processor=platform.processor(),
                cpu_count=os.cpu_count(), threads={key: os.environ[key] for key in
                ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "SAGE_NUM_THREADS")})
    if sys.platform == "darwin":
        for name, key in (("processor", "machdep.cpu.brand_string"),
                          ("model", "hw.model"), ("memory_bytes", "hw.memsize")):
            result = subprocess.run(["/usr/sbin/sysctl", "-n", key],
                                    capture_output=True, text=True)
            if result.returncode == 0:
                info[name] = result.stdout.strip()
    return info


def sample_seeds(seed, code_id, weight, trial, phase="trial"):
    """Derive separate message/error seeds without overlapping adjacent roots."""
    prefix = f"{phase}:{seed}:{code_id}:{weight}:{trial}"
    return tuple(int.from_bytes(hashlib.sha256(f"{prefix}:{kind}".encode()).digest(),
                                "big") for kind in ("message", "error"))


def timing_summary(times):
    row = {name + "_ms": statistics.fmean(times[name]) / 10**6 for name in METHODS}
    row.update({name + "_median_ms": statistics.median(times[name]) / 10**6
                for name in METHODS})
    return row


def run(output, code_ids, trials, seeds=DEFAULT_SEEDS):
    output.mkdir(parents=True, exist_ok=True)
    paths = {name: output / name for name in
             ("individual_times.csv", "summary_by_seed.csv", "summary.csv", "metadata.json")}
    if any(path.exists() for path in paths.values()):
        raise FileExistsError("Use a new output directory for each benchmark run.")
    combinations = sum((len(POINTS[code_id // 100]) - code_id % 100) // 2
                       for code_id in code_ids)
    source_dir = Path(__file__).resolve().parent
    metadata = dict(started_utc=datetime.now(timezone.utc).isoformat(),
                    completed=False, environment=system_info(), code_ids=list(code_ids),
                    seeds=list(seeds), trials_per_seed_per_weight=trials,
                    trials_per_weight=len(seeds) * trials, methods=list(METHODS),
                    seed_policy="All specified seeds are reported; none is selected by measured performance",
                    rng="Python random.Random; finite-field elements use base-p polynomial encoding",
                    seed_rule="SHA256(phase:seed:code_id:w:trial:kind), interpreted as a big-endian integer; phase=trial or warmup, kind=message or error",
                    timing="perf_counter_ns around one complete decoder call; fixed setup, sampling, correctness checks and output excluded",
                    warmup="Two independent words per code and positive weight before all seeds, without/with an infinity error; root=20261017, trials=0,1, phase=warmup",
                    ordering="Seeds rotate by trial; methods rotate by trial, seed index, code ID and error weight",
                    aggregation="Arithmetic mean and median of all per-call times; summary_by_seed separates seeds, summary pools every specified seed",
                    expected_parameter_weight_combinations=combinations,
                    expected_seed_parameter_weight_combinations=combinations * len(seeds),
                    expected_inputs=combinations * len(seeds) * trials,
                    expected_decodings=combinations * len(seeds) * trials * len(METHODS),
                    expected_warmup_decodings=combinations * 2 * len(METHODS),
                    source_sha256={name: hashlib.sha256((source_dir / name).read_bytes()).hexdigest()
                                   for name in ("benchmark.py", "etgrs_decoding.py")},
                    verified_decodings=0, warmup_decodings=0,
                    parameter_weight_combinations=0, seed_parameter_weight_combinations=0)
    paths["metadata.json"].write_text(json.dumps(metadata, indent=2) + "\n")
    parameters = ["code_id", "family", "q", "n", "k", "N", "t"]
    columns = parameters + ["seed", "w", "trial", "codeword_seed", "error_seed",
                            "infinity_error", "method", "elapsed_ns"]
    time_columns = [name + "_ms" for name in METHODS] + [name + "_median_ms" for name in METHODS]
    summary_columns = parameters + ["w", "trials"] + time_columns
    seed_columns = parameters + ["seed", "w", "trials"] + time_columns
    try:
        with paths["individual_times.csv"].open("w", newline="") as raw_file, \
                paths["summary_by_seed.csv"].open("w", newline="") as seed_file, \
                paths["summary.csv"].open("w", newline="") as summary_file:
            raw_writer = csv.DictWriter(raw_file, fieldnames=columns)
            seed_writer = csv.DictWriter(seed_file, fieldnames=seed_columns)
            summary_writer = csv.DictWriter(summary_file, fieldnames=summary_columns)
            raw_writer.writeheader()
            seed_writer.writeheader()
            summary_writer.writeheader()
            for code_id in code_ids:
                case = make_case(code_id)
                S, v, k, eta = (case[key] for key in ("S", "v", "k", "eta"))
                decoders = {"ECP": ECPDecoder(S, v, k, eta)}
                decoders.update({name: BeelenDecoder(S, v, k, eta, method=name)
                                 for name in METHODS[1:]})
                functions = {name: decoder.decode for name, decoder in decoders.items()}
                params = {key: case[key] for key in parameters}
                for weight in range(1, case["t"] + 1):
                    for j in range(2):
                        codeword_seed, error_seed = sample_seeds(20261017, code_id, weight,
                                                                j, phase="warmup")
                        sent, received, _ = received_word(case, weight, codeword_seed,
                                                        error_seed, infinity_error=bool(j))
                        for function in functions.values():
                            if function(received) != sent:
                                raise AssertionError((code_id, weight, "warmup"))
                            metadata["warmup_decodings"] += 1
                    times = {seed: {name: [] for name in METHODS} for seed in seeds}
                    for trial in range(1, trials + 1):
                        for offset in range(len(seeds)):
                            seed_index = (offset + trial - 1) % len(seeds)
                            seed = seeds[seed_index]
                            codeword_seed, error_seed = sample_seeds(seed, code_id, weight, trial)
                            sent, received, flag = received_word(case, weight, codeword_seed,
                                                                error_seed)
                            shift = (trial - 1 + seed_index + code_id + weight) % len(METHODS)
                            order = METHODS[shift:] + METHODS[:shift]
                            for name in order:
                                function = functions[name]
                                start = time.perf_counter_ns()
                                decoded = function(received)
                                elapsed = time.perf_counter_ns() - start
                                if decoded != sent or elapsed <= 0:
                                    raise AssertionError((code_id, weight, seed, trial, name))
                                times[seed][name].append(elapsed)
                                metadata["verified_decodings"] += 1
                                raw_writer.writerow(dict(params, seed=seed, w=weight, trial=trial,
                                    codeword_seed=codeword_seed, error_seed=error_seed,
                                    infinity_error=flag, method=name, elapsed_ns=elapsed))
                    for seed in seeds:
                        seed_writer.writerow(dict(params, seed=seed, w=weight, trials=trials,
                                                  **timing_summary(times[seed])))
                        metadata["seed_parameter_weight_combinations"] += 1
                    pooled = {name: [elapsed for seed in seeds for elapsed in times[seed][name]]
                              for name in METHODS}
                    summary_writer.writerow(dict(params, w=weight, trials=trials * len(seeds),
                                                 **timing_summary(pooled)))
                    metadata["parameter_weight_combinations"] += 1
                    raw_file.flush()
                    seed_file.flush()
                    summary_file.flush()
                print(f"{code_id}: completed {case['t']} error weights for {len(seeds)} seeds", flush=True)
        assert metadata["verified_decodings"] == metadata["expected_decodings"]
        assert metadata["seed_parameter_weight_combinations"] == metadata["expected_seed_parameter_weight_combinations"]
        metadata["completed"] = True
    finally:
        metadata["finished_utc"] = datetime.now(timezone.utc).isoformat()
        paths["metadata.json"].write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--trials", type=int, default=10, help="trials per seed and error weight")
    parser.add_argument("--seeds", nargs="+", type=int, default=DEFAULT_SEEDS,
                        help="root seeds; all supplied seeds are reported")
    parser.add_argument("--codes", default=",".join(map(str, CODE_IDS)),
                        help="comma-separated code IDs")
    args = parser.parse_args()
    code_ids = tuple(int(value) for value in args.codes.split(","))
    seeds = tuple(args.seeds)
    if args.trials < 1 or len(set(code_ids)) != len(code_ids) or \
            not code_ids or any(code_id not in CODE_IDS for code_id in code_ids) or \
            not seeds or len(set(seeds)) != len(seeds):
        parser.error("Use positive trials, distinct seeds and distinct code IDs from the nine listed cases.")
    run(args.output, code_ids, args.trials, seeds)
