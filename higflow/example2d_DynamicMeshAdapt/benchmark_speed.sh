#!/bin/bash

RUNS=5
PIN=0

OUTFILE="benchmark_results.txt"
> "$OUTFILE"

BASE_IN="mult newt-newt De1=1.0 beta1=0.5 Re=0.005 dt=0.001 numsteps=1001 dts=1.0 dtp=0.01 shearing_droplet"

# config format: adapt mesh h_min label_extra
CONFIGS=(
    "true  square_61_br  0.00409842  freq=15"
    "true  square_61_br  0.00409842  freq=5"
    "false square_241_br 0.00414944  uniform"
)

for config in "${CONFIGS[@]}"; do
    adapt=$(echo "$config" | awk '{print $1}')
    mesh=$(echo "$config"  | awk '{print $2}')
    elem=$(echo "$config"  | awk '{print $3}')
    extra=$(echo "$config" | awk '{print $4}')

    if [ "$extra" = "uniform" ]; then
        label="Uniform Fine (square_241_br)"
        extra_params="adapt=false"
    else
        freq="${extra#freq=}"
        label="Adaptive freq=${freq} (square_61_br)"
        extra_params="adapt=true adapt_freq=${freq}"
    fi

    total_sum=0
    iter_sum=0
    valid_runs=0

    echo "======================================================================"
    echo "  Case: $label"
    echo "  Mesh: $mesh  |  $extra_params  |  h_min=$elem"
    echo "======================================================================"

    for i in $(seq 1 $RUNS); do
        ename="bench_$(date +%Y%m%d_%H%M%S)_r${i}"
        echo "  Run $i/$RUNS (ENAME=$ename) ..."

        if [ $i -eq 1 ]; then
            output=$(taskset -c $PIN make build_run \
                IN="$BASE_IN $extra_params" \
                MESH="$mesh" NP=1 ENAME="$ename" 2>&1)
            ret=$?
        else
            output=$(taskset -c $PIN make run \
                IN="$BASE_IN $extra_params" \
                MESH="$mesh" NP=1 ENAME="$ename" 2>&1)
            ret=$?
        fi

        if [ $ret -ne 0 ]; then
            echo "    !! make failed (code $ret). Skipping."
            continue
        fi

        total_t=$(echo "$output" | grep "total time"         | tail -1 | grep -oP '\d+\.\d+' | head -1)
        iter_t=$( echo "$output" | grep "time of iteration"   | tail -1 | grep -oP '\d+\.\d+' | head -1)

        if [ -z "$total_t" ] || [ -z "$iter_t" ]; then
            echo "    !! Failed to extract timings. Skipping."
            continue
        fi

        echo "    total=$total_t s  |  iter=$iter_t s"
        echo "  Run $i: total=$total_t  iter=$iter_t" >> "$OUTFILE"

        total_sum=$(echo "$total_sum + $total_t" | bc -l)
        iter_sum=$( echo "$iter_sum  + $iter_t"  | bc -l)
        valid_runs=$((valid_runs + 1))
    done

    if [ $valid_runs -gt 0 ]; then
        total_avg=$(echo "scale=6; $total_sum / $valid_runs" | bc -l)
        iter_avg=$( echo "scale=6; $iter_sum  / $valid_runs" | bc -l)
    else
        total_avg=0
        iter_avg=0
    fi

    echo "  >> Average: total=$total_avg s  |  iter=$iter_avg s"
    echo ""

    # store for table
    labels[$run_idx]="$label"
    element_sizes[$run_idx]="$elem"
    total_avgs[$run_idx]=$total_avg
    iter_avgs[$run_idx]=$iter_avg
    run_idx=$((run_idx + 1))
done

echo ""
echo "======================================================================"
echo "                      SPEEDUP TABLE"
echo "======================================================================"
echo ""

# Reference: uniform case is the last one (index 2)
ref_idx=2
base_total=${total_avgs[$ref_idx]}

{
echo "Smallest element (h_min):"
printf "  - Adaptive  (square_61_br):  %s\n" "${element_sizes[0]}"
printf "  - Uniform   (square_241_br): %s\n" "${element_sizes[$ref_idx]}"
echo ""
printf "%-45s %-18s %-18s %-12s\n" "Case" "Total (s)" "Iter (s)" "Speedup"
echo "------------------------------------------------------------------------------------------"
for i in $(seq 0 $((run_idx - 1))); do
    if [ "$i" -eq "$ref_idx" ]; then
        speedup_str="1.000 (ref)"
    else
        speedup=$(echo "scale=3; $base_total / ${total_avgs[$i]}" | bc -l)
        speedup_str=$(printf "%.3f" "$speedup")
    fi
    printf "%-45s %-18.6f %-18.6f %-12s\n" "${labels[$i]}" "${total_avgs[$i]}" "${iter_avgs[$i]}" "$speedup_str"
done
echo "------------------------------------------------------------------------------------------"
echo "Speedup = UniformFine_time / Case_time"
} | tee -a "$OUTFILE"

echo ""
echo "Results saved to: $OUTFILE"
