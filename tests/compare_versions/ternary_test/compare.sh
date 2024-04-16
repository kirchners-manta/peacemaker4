#!/bin/bash

# Change directory to the peacemaker4 results
cd ../../../working-code/test-calc/ternary/

# Set the output file to a specific path (change this path to your desired location)
output_file=../../../tests/compare_versions/ternary_test/differences.txt

# Remove the existing differences.txt file
rm -f "$output_file"

# Run the calculation 
../../pm_meson/build/peacemaker qce.toml cmw.clusterset.toml

# Look for the result files
for file in ../../../peacemaker3/test-calc/ternary/*.dat ; do

    filename=$(basename "$file")

    echo "Differences for file: $filename" >> "$output_file"

    # compare the results
    diff "$file" "$filename" >> "$output_file"
done

echo "Comparison completed. Differences saved in differences.txt"
