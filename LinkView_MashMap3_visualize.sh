#!/bin/bash

# Usage: ./script.sh -q <query.fa> -r <reference.fa>

# Parse command-line arguments
while getopts ":q:r:" opt; do
    case ${opt} in
        q )
            QUERY=${OPTARG}
            ;;
        r )
            REFERENCE=${OPTARG}
            ;;
        \? )
            echo "Invalid option: -$OPTARG" 1>&2
            exit 1
            ;;
        : )
            echo "Invalid option: -$OPTARG requires an argument" 1>&2
            exit 1
            ;;
    esac
done

# Ensure both inputs are provided
if [ -z "${QUERY}" ] || [ -z "${REFERENCE}" ]; then
    echo "Usage: $0 -q <query.fa> -r <reference.fa>"
    exit 1
fi

# Generate output file prefixes
QUERY_PREFIX=$(basename "${QUERY}" .fa)
REFERENCE_PREFIX=$(basename "${REFERENCE}" .fa)

# Step 1: Run wfmash
echo "Running wfmash..."
wfmash "${REFERENCE}" "${QUERY}" -s 50k -l 150k -p 90 -n 1 -H 1 -m -t 100 > "${REFERENCE_PREFIX}_vs_${QUERY_PREFIX}.out" -N

# Step 2: Extract unique contig names
echo "Extracting unique contig names..."
grep ">" "${REFERENCE}" | cut -d ">" -f2 > "${REFERENCE_PREFIX}_contig_list.txt"

# Step 3: Create subset mashmap outputs
echo "Creating subset mashmap outputs..."
while read -r contig; do
    grep "${contig}" "${REFERENCE_PREFIX}_vs_${QUERY_PREFIX}.out" > "${contig}_vs_qry.subset.mashmapout"
done < "${REFERENCE_PREFIX}_contig_list.txt"

# Step 4: Combine mashmap outputs by chromosome index
echo "Combining mashmap outputs..."
grep ">" "${QUERY}" | cut -f5,6 -d "_" | sort -u > "${REFERENCE_PREFIX}_Chr_indexes.txt"

while read -r index; do
    cat *"${index}"*mashmapout > "${index}_combined.results"
done < "${REFERENCE_PREFIX}_Chr_indexes.txt"

# Step 5: Sort and process combined results for plotting
echo "Sorting and processing results..."
while read -r Chr; do
    cat "${Chr}_combined.results" | cut -f1,3,4,5,6,7,8,9,10 | cut -f2 -d ":" | sort -k5,5 -k8,8n > "table${Chr}"
done < "${REFERENCE_PREFIX}_Chr_indexes.txt"

# Prepare input for LINKVIEW (mapping and karyo files)
echo "Preparing LINKVIEW inputs..."
for file in *table*; do
    cat "$file" | cut -f1,2,3,5,7,8 > "linkview_${file}.tsv"
    cat "$file" | cut -f4,5 > "linkkaryo_${file}.tsv"
done

# Step 6: Run LINKVIEW with mapping table
echo "Running LINKVIEW with mapping table..."
for file in *linkview*; do
    # Extract unique values from column 4
    cat "$file" | cut -f4 | sort -u > "${file%.tsv}.unique_idx"

    # Process each unique value
    while read -r line; do
        # Filter the file by the current line value
        grep "$line" "$file" > "${file%.tsv}.${line}.tsv"

        # Run LINKVIEW for the filtered file
        python LINKVIEW/LINKVIEW.py -t 0 --svg2png cairosvg -o "${file%.tsv}_${line}.png" "${file%.tsv}.${line}.tsv"
    done < "${file%.tsv}.unique_idx"
done

# Step 7: Run LINKVIEW with karyo file
echo "Running LINKVIEW with karyo file..."
for file in *linkview*; do
    # Extract unique values from column 4
    cat "$file" | cut -f4 | sort -u > "${file%.tsv}.unique_idx"

    # Process each unique value
    while read -r line; do
        # Filter the file by the current line value
        grep "$line" "$file" > "${file%.tsv}.${line}.tsv"

        # Get corresponding karyo file
        karyo_file="linkkaryo_${file}.tsv"

        # Run LINKVIEW for the filtered file
        python LINKVIEW/LINKVIEW.py -t 0 --chro_len "$karyo_file" --svg2png cairosvg -o "${file%.tsv}_${line}.png" "${file%.tsv}.${line}.tsv"
    done < "${file%.tsv}.unique_idx"