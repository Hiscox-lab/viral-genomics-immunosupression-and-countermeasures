#!/bin/bash

# Define thresholds
THRESHOLDS=(10 50 100 250 1000)

# Output file
OUTPUT_FILE="coverage_results.csv"

# Create CSV header
echo "Filename,Threshold,Breadth_of_Coverage,Average_Depth" > $OUTPUT_FILE

THRESHOLDS=(10 50 100 200 500 1000)

# Loop through all BAM files in the current directory
for BAM_FILE in *.bam; do
    echo "Processing $BAM_FILE..."

    # Generate per-base depth using samtools
    samtools depth $BAM_FILE > temp_depth.txt

    # Get the total number of reference bases and the sum of depths
    TOTAL_BASES=$(awk '{count++} END {print count}' temp_depth.txt)
    SUM_DEPTH=$(awk '{sum+=$3} END {print sum}' temp_depth.txt)

    # Calculate average depth
    if [[ $TOTAL_BASES -gt 0 ]]; then
        AVERAGE_DEPTH=$(echo "scale=2; $SUM_DEPTH / $TOTAL_BASES" | bc)
    else
        AVERAGE_DEPTH=0
    fi

    # Calculate breadth of coverage for each threshold
    for THRESHOLD in "${THRESHOLDS[@]}"; do
        COVERED_BASES=$(awk -v threshold=$THRESHOLD '$3 >= threshold {count++} END {print count}' temp_depth.txt)
        if [[ $TOTAL_BASES -gt 0 ]]; then
            BREADTH=$(echo "scale=4; $COVERED_BASES / $TOTAL_BASES * 100" | bc)
        else
            BREADTH=0
        fi
        
        # Append results to CSV
        echo "$BAM_FILE,$THRESHOLD,$BREADTH,$AVERAGE_DEPTH" >> $OUTPUT_FILE
    done

    # Clean up temporary depth file
    rm temp_depth.txt
done

echo "Coverage analysis completed. Results saved to $OUTPUT_FILE."
