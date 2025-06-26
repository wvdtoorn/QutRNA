#!/bin/bash

# Minimal QutRNA Workflow Runner
# This script runs the minimal QutRNA workflow

set -e  # Exit on any error

echo "Starting Minimal QutRNA Workflow..."

# Check if required files exist
if [ ! -f "Snakefile" ]; then
    echo "Error: Snakefile not found!"
    exit 1
fi

if [ ! -f "minimal_config.yaml" ]; then
    echo "Error: minimal_config.yaml not found!"
    exit 1
fi

# Run the workflow with both config files
echo "Running Snakemake workflow..."
snakemake -c 1 --configfile=minimal_config.yaml

echo "Minimal QutRNA workflow completed successfully!"
echo "Check the 'config[qutrna][output_dir]/results/' directory for output files." 