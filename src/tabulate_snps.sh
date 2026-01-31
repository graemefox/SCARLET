#!/bin/bash

input_file="$1"
prefix="$2"

# Get headers (remove quotes)
header_line=$(head -n 1 "$input_file" | sed 's/^"//;s/"$//' | tr '","' '\n')
# Read into array
IFS=$'\n' read -d '' -r -a headers <<< "$header_line"

# Process data lines
tail -n +2 "$input_file" | while IFS=',' read -r -a fields; do
    # Strip quotes from each field
    for i in "${!fields[@]}"; do
        fields[$i]=$(echo "${fields[$i]}" | sed 's/^"//;s/"$//; s/,/ /g')
    done

    # Print prefix
    printf "%s" "$prefix"

    # Print key-value pairs
    for i in "${!headers[@]}"; do
        printf "%s %s" "${headers[$i]}" "${fields[$i]}"
    done

    printf "\n"
done
