#!/bin/bash

# List of read files
read_files=(
    "../00-reads/hifi/CG_F_HIFI_None_reads.fastq.gz"
    "../00-reads/hifi/CM_M_HIFI_None_reads.fastq.gz")

    
# Loop through each file and generate a separate CSV file for each
for file in "${read_files[@]}"; do
    # Extract the filename without extension for output naming
    filename=$(basename "$file" .fastq.gz)

    # Define the output CSV file for the current FASTQ file
    output_file="${filename}_read_lengths_and_phred_scores.csv"

    echo "Processing $file ..."

    # Add a header to the output CSV file
#    print Read \t Length_(bp) \t Mean_PHRED_Score \t GC_content" > "$output_file"
# Add a header to the output CSV file
echo -e "Read\tLength_(bp)\tGC_content\tMean_PHRED_Score" > "$output_file"

    # Calculate read lengths and mean PHRED scores using seqkit
#    seqkit fx2tab --name --only-id --length --avg-qual -g -G "$file" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}'  >> "$output_file"
     seqkit fx2tab --name --only-id --length --avg-qual -g "$file" | awk '{print $1"\t"$2"\t"$3"\t"$4}'  >> "$output_file"
    echo "Results saved to $output_file"
done

echo "Done!"