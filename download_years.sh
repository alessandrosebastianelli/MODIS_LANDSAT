#!/bin/bash

# Create log directories if they do not exist
mkdir -p logs/out logs/err

first=true

# Submit jobs for years 2007–2021
for year in {2007..2021}; do
    job_name="pipeline_${year}"

    echo "Submitting job ${job_name}"

    if [ "$first" = true ]; then
        # First iteration: no --skip-static
        bsub \
            -J "${job_name}" \
            -q s_long \
            -M 100G \
            -P R000 \
            -oo "logs/out/${job_name}.out" \
            -eo "logs/err/${job_name}.err" \
            python run_pipeline_yearly.py --year ${year}

        first=false
    else
        # All following iterations: add --skip-static
        bsub \
            -J "${job_name}" \
            -q s_long \
            -M 100G \
            -P R000 \
            -oo "logs/out/${job_name}.out" \
            -eo "logs/err/${job_name}.err" \
            python run_pipeline_yearly.py --year ${year} --skip-static --skip-validation
    fi
done