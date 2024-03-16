#!/bin/bash

for file in "$@"; do
    if [ -f "$file" ]; then
        filename=$(basename "$file")

        fourth_column=$(awk 'NR==1{print $4; exit}' "$file")

        new_filename="${fourth_column}.txt"

        cp "$file" "$new_filename"

        echo "done done done"
    else
        echo "error"
    fi
done
