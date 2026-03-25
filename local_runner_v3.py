#!/usr/bin/env python3
import argparse
import subprocess
import sys
import tempfile
import time
from pathlib import Path


ROOT = Path("/home/kamisama/huawei_test")
DEFAULT_SOLUTION = ROOT / "demos" / "Solution.cpp"
DEFAULT_INPUT = ROOT / "data" / "practice_1.in"
DEFAULT_OUTPUT = ROOT / "results" / "practice_1.v3.out"


def detect_lang(source: Path) -> str:
    suffix = source.suffix.lower()
    if suffix == ".c":
        return "c"
    if suffix == ".cpp":
        return "cpp"
    if suffix == ".java":
        return "java"
    if suffix == ".py":
        return "python"
    raise ValueError(f"unsupported source extension: {source}")


def compile_source(source: Path, build_dir: Path):
    lang = detect_lang(source)
    if lang == "c":
        exe = build_dir / "Solution"
        cmd = ["gcc", "-std=c17", str(source), "-o", str(exe), "-O2", "-lpthread", "-lm"]
        subprocess.run(cmd, check=True, text=True)
        return lang, [str(exe)], cmd
    if lang == "cpp":
        exe = build_dir / "Solution"
        cmd = ["g++", "-std=c++17", str(source), "-o", str(exe), "-O2", "-lpthread"]
        subprocess.run(cmd, check=True, text=True)
        return lang, [str(exe)], cmd
    if lang == "java":
        cmd = ["javac", "-implicit:none", "-d", str(build_dir), str(source)]
        subprocess.run(cmd, check=True, text=True)
        return lang, ["java", "-cp", str(build_dir), "Solution"], cmd
    cmd = ["python3", str(source)]
    return lang, cmd, None


def parse_static_case(input_path: Path):
    tokens = input_path.read_text(encoding="utf-8").split()
    idx = 0

    def need(count: int):
        if idx + count > len(tokens):
            raise ValueError(f"input truncated near token index {idx} in {input_path}")

    need(2)
    n_a = int(tokens[idx])
    n_b = int(tokens[idx + 1])
    idx += 2

    need(2 * n_a)
    poly_a = []
    for _ in range(n_a):
        poly_a.append((float(tokens[idx]), float(tokens[idx + 1])))
        idx += 2

    need(2 * n_b)
    poly_b = []
    for _ in range(n_b):
        poly_b.append((float(tokens[idx]), float(tokens[idx + 1])))
        idx += 2

    need(1)
    k = int(tokens[idx])
    idx += 1

    need(2 * k)
    queries = []
    for _ in range(k):
        queries.append((float(tokens[idx]), float(tokens[idx + 1])))
        idx += 2

    if idx != len(tokens):
        raise ValueError(f"unexpected extra tokens in {input_path}: parsed {idx}, total {len(tokens)}")

    return poly_a, poly_b, queries


def format_payload(poly_a, poly_b, queries):
    parts = [f"{len(poly_a)} {len(poly_b)}\n"]
    parts.append(" ".join(f"{x:.5f} {y:.5f}" for x, y in poly_a) + "\n")
    parts.append(" ".join(f"{x:.5f} {y:.5f}" for x, y in poly_b) + "\n")
    parts.append("OK\n")
    parts.append(f"{len(queries)}\n")
    for x, y in queries:
        parts.append(f"{x:.5f} {y:.5f}\n")
    parts.append("OK\n")
    return "".join(parts)


def parse_solver_stdout(stdout: str, query_count: int):
    lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    expected = query_count + 3
    if len(lines) != expected:
        raise ValueError(f"protocol line count mismatch: expected {expected}, got {len(lines)}")
    if lines[0] != "OK":
        raise ValueError(f"missing preprocess OK, got first line: {lines[0] if lines else '<empty>'}")
    if lines[1] != str(query_count):
        raise ValueError(f"reported answer count mismatch: expected {query_count}, got {lines[1]}")
    if lines[-1] != "OK":
        raise ValueError(f"missing final OK, got last line: {lines[-1]}")

    answers = []
    for line in lines[2:-1]:
        cols = line.split()
        if len(cols) != 2:
            raise ValueError(f"bad answer line: {line}")
        answers.append((float(cols[0]), float(cols[1])))
    if len(answers) != query_count:
        raise ValueError(f"answer count mismatch after parse: expected {query_count}, got {len(answers)}")
    return lines, answers


def main():
    parser = argparse.ArgumentParser(description="Local protocol-v3 runner for Huawei challenge static inputs.")
    parser.add_argument("solution", nargs="?", default=str(DEFAULT_SOLUTION))
    parser.add_argument("input", nargs="?", default=str(DEFAULT_INPUT))
    parser.add_argument("output", nargs="?", default=str(DEFAULT_OUTPUT))
    parser.add_argument("--timeout-seconds", type=float, default=620.0)
    parser.add_argument("--keep-stderr", action="store_true", help="write solver stderr to <output>.stderr")
    args = parser.parse_args()

    solution = Path(args.solution).resolve()
    input_path = Path(args.input).resolve()
    output_path = Path(args.output).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    poly_a, poly_b, queries = parse_static_case(input_path)
    payload = format_payload(poly_a, poly_b, queries)

    with tempfile.TemporaryDirectory(prefix="local_runner_v3_") as tmp:
        build_dir = Path(tmp)
        lang, run_cmd, compile_cmd = compile_source(solution, build_dir)

        start_ns = time.perf_counter_ns()
        proc = subprocess.run(
            run_cmd,
            input=payload,
            capture_output=True,
            text=True,
            check=False,
            timeout=args.timeout_seconds,
        )
        elapsed_ns = time.perf_counter_ns() - start_ns

        if args.keep_stderr or proc.stderr.strip():
            output_path.with_suffix(output_path.suffix + ".stderr").write_text(proc.stderr, encoding="utf-8")

        if proc.returncode != 0:
            output_path.write_text(proc.stdout, encoding="utf-8")
            print(f"Run failed with exit code {proc.returncode}", file=sys.stderr)
            if proc.stderr.strip():
                print(proc.stderr.strip(), file=sys.stderr)
            return proc.returncode

        try:
            lines, _answers = parse_solver_stdout(proc.stdout, len(queries))
        except Exception as exc:
            output_path.write_text(proc.stdout, encoding="utf-8")
            print(f"Protocol check failed: {exc}", file=sys.stderr)
            return 3

        answer_stage = "\n".join(lines[1:]) + "\n"
        output_path.write_text(answer_stage, encoding="utf-8")

        print(f"Language: {lang}")
        if compile_cmd is not None:
            print("Compile:", " ".join(compile_cmd))
        print("Run:", " ".join(run_cmd))
        print(f"Input: {input_path}")
        print(f"Output: {output_path}")
        print(f"Queries: {len(queries)}")
        print(f"Elapsed: {elapsed_ns} ns")
        print("Protocol: OK")
        return 0


if __name__ == "__main__":
    sys.exit(main())
