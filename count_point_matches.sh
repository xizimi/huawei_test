#!/usr/bin/env bash
# Compare per-point outputs between C (reference) and C++ (candidate)
# Usage: ./scripts/count_point_matches.sh [tolerance]

ROOT_DIR="/home/zhaobowen/project/c++_project/Linux/"
RESULTS="$ROOT_DIR/results/compare_tmp"
TOL=${1:-0.0001}

if [ ! -d "$RESULTS" ]; then
    echo "Results directory not found: $RESULTS" >&2
    exit 1
fi

total_points=0
total_matched=0

printf "Comparing with tolerance = %s\n\n" "$TOL"

for ref in "$RESULTS"/*.c.out; do
    base=$(basename "$ref" .c.out)
    cand="$RESULTS/${base}.cpp.out"
    if [ ! -f "$cand" ]; then
        echo "Candidate output missing for $base" >&2
        continue
    fi

    # Skip header line if it's a single integer (number of points)
    header_ref=$(head -n1 "$ref")
    header_cand=$(head -n1 "$cand")
    skip_header=0
    if [[ "$header_ref" =~ ^[0-9]+$ ]] && [[ "$header_cand" =~ ^[0-9]+$ ]]; then
        data_ref_cmd="tail -n +2 '$ref'"
        data_cand_cmd="tail -n +2 '$cand'"
    else
        data_ref_cmd="cat '$ref'"
        data_cand_cmd="cat '$cand'"
    fi

    # combine lines and compare two floats per line
    matched=0
    points=0
    eval "$data_ref_cmd" | paste -d '|' - <(eval "$data_cand_cmd") | \
    awk -v tol="$TOL" -F'|' '
    {
        points++
        split($1,a,"[ \t]+")
        split($2,b,"[ \t]+")
        x1=a[1]+0; y1=a[2]+0; x2=b[1]+0; y2=b[2]+0;
        if ( (x1==x2 && y1==y2) || ( (sqrt((x1-x2)*(x1-x2))<=tol) && (sqrt((y1-y2)*(y1-y2))<=tol) ) ) matched++;
    }
    END { printf("%d %d\n", points, matched) }'
    ret=$?

    if [ $ret -ne 0 ]; then
        echo "Error processing $base" >&2
        continue
    fi

    # capture the last printed line (points matched)
    last=$(eval "$data_ref_cmd" | paste -d '|' - <(eval "$data_cand_cmd") | \
    awk -v tol="$TOL" -F'|' '{points++; split($1,a,"[ \t]+"); split($2,b,"[ \t]+"); x1=a[1]+0; y1=a[2]+0; x2=b[1]+0; y2=b[2]+0; if ((x1==x2 && y1==y2) || ((sqrt((x1-x2)^2)<=tol) && ((sqrt((y1-y2)^2))<=tol))) matched++; } END { printf("%d %d\n", points, matched) }')

    pts=$(echo "$last" | awk '{print $1}')
    m=$(echo "$last" | awk '{print $2}')

    total_points=$((total_points + pts))
    total_matched=$((total_matched + m))

    perc=$(awk "BEGIN{if($pts>0) printf(\"%.2f\", ($m/$pts)*100); else print \"0.00\"}")
    printf "%s: %d / %d = %s%%\n" "$base" "$m" "$pts" "$perc"
done

overall_perc=$(awk "BEGIN{if($total_points>0) printf(\"%.2f\", ($total_matched/$total_points)*100); else print \"0.00\"}")
printf "\nOverall: %d / %d = %s%%\n" "$total_matched" "$total_points" "$overall_perc"

exit 0
