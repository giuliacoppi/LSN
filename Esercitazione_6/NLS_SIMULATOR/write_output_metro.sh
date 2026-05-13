#!/bin/bash

mkdir -p ./DATA_ANALYSIS

# output files names
OUTPUT_U="./DATA_ANALYSIS/metropolis_total_energy.dat"
OUTPUT_CV="./DATA_ANALYSIS/metropolis_specific_heat.dat"
OUTPUT_CHI="./DATA_ANALYSIS/metropolis_susceptibility.dat"
OUTPUT_MAGNET="./DATA_ANALYSIS/metropolis_magnetization.dat"

printf "# %-8s %-15s %-15s\n" "T" "Average" "Error" > "$OUTPUT_U"
printf "# %-8s %-15s %-15s\n" "T" "Average" "Error" > "$OUTPUT_CV"
printf "# %-8s %-15s %-15s\n" "T" "Average" "Error" > "$OUTPUT_CHI"
printf "# %-8s %-15s %-15s\n" "T" "Average" "Error" > "$OUTPUT_MAGNET"

for T in $(seq 0.5 0.1 2.0)
do
    FILE_U="./RESULTS/METROPOLIS/T_$T/total_energy.dat"
    FILE_CV="./RESULTS/METROPOLIS/T_$T/specific_heat.dat"
    FILE_CHI="./RESULTS/METROPOLIS/T_$T/susceptibility.dat"
    FILE_M="./RESULTS/METROPOLIS/MAGNET/T_$T/magnetization.dat"
    
    if [ -f "$FILE_U" ] && [ -f "$FILE_CV" ] && [ -f "$FILE_CHI" ] && [ -f "$FILE_M" ]; then

    	tail -n 1 "$FILE_U" | awk -v t="$T" '{printf "%-10s %-15s %-15s\n", t, $3, $4}' >> "$OUTPUT_U"
    	tail -n 1 "$FILE_CV" | awk -v t="$T" '{printf "%-10s %-15s %-15s\n", t, $3, $4}' >> "$OUTPUT_CV"
    	tail -n 1 "$FILE_CHI" | awk -v t="$T" '{printf "%-10s %-15s %-15s\n", t, $3, $4}' >> "$OUTPUT_CHI"
    	tail -n 1 "$FILE_M" | awk -v t="$T" '{printf "%-10s %-15s %-15s\n", t, $3, $4}' >> "$OUTPUT_MAGNET"
	
    else
	echo "Missinng one or more files for T=$T"
    fi
done

echo "------------------------------------------------"
echo "Generated files: $OUTPUT_U, $OUTPUT_CV, $OUTPUT_CHI, $OUTPUT_MAGNET"
