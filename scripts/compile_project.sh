#!/bin/bash

# Check for correct number of arguments
if [ "$#" -ne 6 ]; then
    echo "Usage: $0 N_X N_Y N_Z SPACING_X SPACING_Y SPACING_Z"
    exit 1
fi

# Arguments
N_X="$1"
N_Y="$2"
N_Z="$3"
SPACING_X="$4"
SPACING_Y="$5"
SPACING_Z="$6"

HEADER_FILE="simulation_constants.hpp"

# Enter include/ to modify the header
echo "Current directory: $(pwd)"
cd include || { echo "include directory not found!"; exit 1; }

# Check if header file exists
if [ ! -f "$HEADER_FILE" ]; then
    echo "Header file $HEADER_FILE not found in include/"
    exit 1
fi

# Backup the original file
cp "$HEADER_FILE" "${HEADER_FILE}.bak"

# Perform replacements
sed -i -E "s/(PHANTOM_N_X\s*=\s*)[0-9]+;/\1$N_X;/" "$HEADER_FILE"
sed -i -E "s/(PHANTOM_N_Y\s*=\s*)[0-9]+;/\1$N_Y;/" "$HEADER_FILE"
sed -i -E "s/(PHANTOM_N_Z\s*=\s*)[0-9]+;/\1$N_Z;/" "$HEADER_FILE"
sed -i -E "s/(PHANTOM_SPACING_X\s*=\s*)[0-9.]+f?;/\1$SPACING_X;/" "$HEADER_FILE"
sed -i -E "s/(PHANTOM_SPACING_Y\s*=\s*)[0-9.]+f?;/\1$SPACING_Y;/" "$HEADER_FILE"
sed -i -E "s/(PHANTOM_SPACING_Z\s*=\s*)[0-9.]+f?;/\1$SPACING_Z;/" "$HEADER_FILE"

echo "Updated $HEADER_FILE successfully."

# Go back to project root
cd ..

# Enter build directory
cd build || { echo "build directory not found!"; exit 1; }

# Run make
echo "Building project..."
make

cd ..

echo "Build finished."
