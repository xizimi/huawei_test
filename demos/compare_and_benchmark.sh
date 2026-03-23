#!/usr/bin/env bash
# Usage: ./scripts/compare_and_benchmark.sh [input_for_benchmark]
# Default benchmark input: data/practice_1.in

ROOT_DIR="/home/zhaobowen/project/c++_project/Linux/"
DEMOS="$ROOT_DIR/demos"
DATA="$ROOT_DIR/data"
RESULTS="$ROOT_DIR/results"
RUN_SCRIPT="$ROOT_DIR/run_and_test.sh"

BENCH_INPUT="${1:-$DATA/practice_1.in}"

mkdir -p "$RESULTS/compare_tmp"

SOL_C="$DEMOS/Solution.c"
SOL_CPP="$DEMOS/Solution.cpp"

echo "Comparing outputs using run_and_test.sh (Interactor/Runner)"
chmod +x "$RUN_SCRIPT" || true
total=0
match=0
for infile in "$DATA"/*.in; do
    ((total++))
    base=$(basename "$infile")
    out_ref="$RESULTS/compare_tmp/${base}.c.out"
    out_cmp="$RESULTS/compare_tmp/${base}.cpp.out"

    echo "Running reference (C) on $base..."
    "$RUN_SCRIPT" "$SOL_C" "$infile" "$out_ref" > "$RESULTS/compare_tmp/${base}.c.log" 2>&1 || echo "Reference run returned non-zero for $base"

    echo "Running candidate (C++) on $base..."
    "$RUN_SCRIPT" "$SOL_CPP" "$infile" "$out_cmp" > "$RESULTS/compare_tmp/${base}.cpp.log" 2>&1 || echo "Candidate run returned non-zero for $base"

    if diff -q "$out_ref" "$out_cmp" >/dev/null 2>&1; then
        ((match++))
    else
        echo "Mismatch for $base: see ${base}.c.out vs ${base}.cpp.out"
    fi
done

accuracy=0
if [ $total -gt 0 ]; then
    accuracy=$(awk "BEGIN {printf \"%.2f\", ($match/$total)*100}")
fi

echo "Accuracy: $match / $total = $accuracy%"

echo "\nRunning $RUN_SCRIPT with Solution.cpp $BENCH_INPUT 10 times to collect Interactor total time..."
times_file="$RESULTS/compare_tmp/interactor_times.txt"
> "$times_file"

# ensure run script is executable
chmod +x "$RUN_SCRIPT" || true

for i in $(seq 1 10); do
    echo "Run #$i..."
    # Run and capture both stdout and stderr
    output=$("$RUN_SCRIPT" "$SOL_CPP" "$BENCH_INPUT" "$RESULTS/compare_tmp/run_${i}.out" 2>&1)
    # Try to extract total time in ns from lines containing "total time"
    # Example line: [INFO]Interactor.cpp:229|main|total time 5503665 ns.
    t=$(echo "$output" | grep -oE "total time [0-9]+ ns" | grep -oE "[0-9]+")
    if [ -z "$t" ]; then
        echo "Warning: could not extract time from run #$i" >&2
        echo "NaN" >> "$times_file"
    else
        echo "$t" >> "$times_file"
    fi
    # small pause
    sleep 0.2
done

# Calculate average ignoring NaN
sum=0
count=0
while read -r line; do
    if [[ "$line" =~ ^[0-9]+$ ]]; then
        sum=$((sum + line))
        count=$((count + 1))
    fi
done < "$times_file"

if [ $count -gt 0 ]; then
    avg=$(awk "BEGIN {printf \"%.2f\", $sum/$count}")
    echo "Average Interactor total time over $count runs: $avg ns"
else
    echo "No valid time samples collected." >&2
fi

echo "Results saved under $RESULTS/compare_tmp/"

exit 0
