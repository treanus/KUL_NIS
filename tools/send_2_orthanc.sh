#!/bin/bash
shopt -s expand_aliases  # Enable alias expansion
source ~/.bashrc         # Or ~/.bash_aliases if you store aliases there

# Set your password here
PASSWORD="changeme"
AET="changeme"
AEC="changeme"
IP="changeme"
PORT="changeme"

# Directory containing 7z files (change if needed)
DIR="./zips"

# Process each .7z file in the directory
for FILE in "$DIR"/*.zip; do
    [ -e "$FILE" ] || continue  # Skip if no .7z files exist
    
    echo "Processing: $FILE"
    BASENAME="$(basename "$FILE" .zip)"
    OUTPUT_DIR="$DIR/$BASENAME"
    mkdir -p "$OUTPUT_DIR"
    
    # Try extracting without password
    7z x "$FILE" -o"$OUTPUT_DIR" -p"" -y &>/dev/null
    
    if [ $? -ne 0 ]; then
        echo "Extraction without password failed. Retrying with password..."
        7z x "$FILE" -o"$OUTPUT_DIR" -p"$PASSWORD" -y &>/dev/null
    fi
    
    if [ $? -eq 0 ]; then
        echo "Extraction successful: $OUTPUT_DIR"
        sleep 5

        # Perform operations on extracted files here
        # Example: list contents
        ls "$OUTPUT_DIR"/*
        dcmsend --scan-directories --recurse -aet $AET -aec $AEC $IP $PORT -v "$OUTPUT_DIR"
        sleep 5

        # Clean up extracted files
        rm -rf "$OUTPUT_DIR"
        echo "Cleaned up: $OUTPUT_DIR"
    else
        echo "Failed to extract: $FILE"
    fi

done
