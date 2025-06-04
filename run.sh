# master script for ATP concentration simulation
# requirements: GCC for C compilation

# takes a command-line argument -- a single string, concatenated with commas or other non-whitespace symbols, determining which aspects of the pipeline will be run
# the options are:

# default     -- default model structure
# smaller     -- smaller cell size
# thinner     -- thinner vertical cell section
# fatter      -- fatter vertical cell section
# more        -- more mitochondria
# fewer       -- fewer mitochondria
# coarser     -- coarser grid elements in numerical scheme

if [ $# -eq 0 ]
then
    echo "No modules selected! Not running anything. See script preamble for options."
    exit 1
fi

commandstr=$1

# compile source code
gcc -o3 mito-agents.c -lm -o mito-agents.ce

if [[ $commandstr == *default* ]]; then
    # default experimental setup
    ./mito-agents.ce 50 10 100 2 0 > tmp-0 &
    ./mito-agents.ce 50 10 100 2 1 > tmp-1 &
    ./mito-agents.ce 50 10 100 2 2 > tmp-2 &
    ./mito-agents.ce 50 10 100 2 3 > tmp-3 &
    ./mito-agents.ce 50 10 100 2 4 5e-6 > tmp-4 &
    ./mito-agents.ce 50 10 100 2 5 5e-6 > tmp-5 &
    ./mito-agents.ce 50 10 100 2 6 1 > tmp-6 &
    ./mito-agents.ce 50 10 100 2 7 1 > tmp-7 &
    ./mito-agents.ce 50 10 100 2 8 0.5 > tmp-8 &
    ./mito-agents.ce 50 10 100 2 9 0.5 > tmp-9 &
fi

if [[ $commandstr == *smaller* ]]; then
    # smaller cell size
    ./mito-agents.ce 20 10 100 2 0 > tmp-0a &
    ./mito-agents.ce 20 10 100 2 1 > tmp-1a &
    ./mito-agents.ce 20 10 100 2 2 > tmp-2a &
    ./mito-agents.ce 20 10 100 2 3 > tmp-3a &
    ./mito-agents.ce 20 10 100 2 4 5e-6 > tmp-4a &
    ./mito-agents.ce 20 10 100 2 5 5e-6 > tmp-5a &
    ./mito-agents.ce 20 10 100 2 6 1 > tmp-6a &
    ./mito-agents.ce 20 10 100 2 7 1 > tmp-7a &
    ./mito-agents.ce 20 10 100 2 8 0.5 > tmp-8a &
    ./mito-agents.ce 20 10 100 2 9 0.5 > tmp-9a &
fi

if [[ $commandstr == *thinner* ]]; then
    # thinner z section
    ./mito-agents.ce 50 2 100 2 0 > tmp-0b &
    ./mito-agents.ce 50 2 100 2 1 > tmp-1b &
    ./mito-agents.ce 50 2 100 2 2 > tmp-2b &
    ./mito-agents.ce 50 2 100 2 3 > tmp-3b &
    ./mito-agents.ce 50 2 100 2 4 5e-6 > tmp-4b &
    ./mito-agents.ce 50 2 100 2 5 5e-6 > tmp-5b &
    ./mito-agents.ce 50 2 100 2 6 1 > tmp-6b &
    ./mito-agents.ce 50 2 100 2 7 1 > tmp-7b &
    ./mito-agents.ce 50 2 100 2 8 0.5 > tmp-8b &
    ./mito-agents.ce 50 2 100 2 9 0.5 > tmp-9b &
fi

if [[ $commandstr == *fatter* ]]; then
    # fatter z section
    ./mito-agents.ce 50 20 100 2 0 > tmp-0c &
    ./mito-agents.ce 50 20 100 2 1 > tmp-1c &
    ./mito-agents.ce 50 20 100 2 2 > tmp-2c &
    ./mito-agents.ce 50 20 100 2 3 > tmp-3c &
    ./mito-agents.ce 50 20 100 2 4 5e-6 > tmp-4c &
    ./mito-agents.ce 50 20 100 2 5 5e-6 > tmp-5c &
    ./mito-agents.ce 50 20 100 2 6 1 > tmp-6c &
    ./mito-agents.ce 50 20 100 2 7 1 > tmp-7c &
    ./mito-agents.ce 50 20 100 2 8 0.5 > tmp-8c &
    ./mito-agents.ce 50 20 100 2 9 0.5 > tmp-9c &
fi

if [[ $commandstr == *more* ]]; then
    # more mitochondria
    ./mito-agents.ce 50 10 200 2 0 > tmp-0d &
    ./mito-agents.ce 50 10 200 2 1 > tmp-1d &
    ./mito-agents.ce 50 10 200 2 2 > tmp-2d &
    ./mito-agents.ce 50 10 200 2 3 > tmp-3d &
    ./mito-agents.ce 50 10 200 2 4 5e-6 > tmp-4d &
    ./mito-agents.ce 50 10 200 2 5 5e-6 > tmp-5d &
    ./mito-agents.ce 50 10 200 2 6 1 > tmp-6d &
    ./mito-agents.ce 50 10 200 2 7 1 > tmp-7d &
    ./mito-agents.ce 50 10 200 2 8 0.5 > tmp-8d &
    ./mito-agents.ce 50 10 200 2 9 0.5 > tmp-9d &
fi

if [[ $commandstr == *fewer* ]]; then
    # more mitochondria
    ./mito-agents.ce 50 10 50 2 0 > tmp-0e &
    ./mito-agents.ce 50 10 50 2 1 > tmp-1e &
    ./mito-agents.ce 50 10 50 2 2 > tmp-2e &
    ./mito-agents.ce 50 10 50 2 3 > tmp-3e &
    ./mito-agents.ce 50 10 50 2 4 5e-6 > tmp-4e &
    ./mito-agents.ce 50 10 50 2 5 5e-6 > tmp-5e &
    ./mito-agents.ce 50 10 50 2 6 1 > tmp-6e &
    ./mito-agents.ce 50 10 50 2 7 1 > tmp-7e &
    ./mito-agents.ce 50 10 50 2 8 0.5 > tmp-8e &
    ./mito-agents.ce 50 10 50 2 9 0.5 > tmp-9e &
fi


if [[ $commandstr == *coarser* ]]; then
    # coarser simulation grid
    ./mito-agents.ce 20 10 100 1 0 > tmp-0f &
    ./mito-agents.ce 20 10 100 1 1 > tmp-1f &
    ./mito-agents.ce 20 10 100 1 2 > tmp-2f &
    ./mito-agents.ce 20 10 100 1 3 > tmp-3f &
    ./mito-agents.ce 20 10 100 1 4 5e-6 > tmp-4f &
    ./mito-agents.ce 20 10 100 1 5 5e-6 > tmp-5f &
    ./mito-agents.ce 20 10 100 1 6 1 > tmp-6f &
    ./mito-agents.ce 20 10 100 1 7 1 > tmp-7f &
    ./mito-agents.ce 20 10 100 1 8 0.5 > tmp-8f &
    ./mito-agents.ce 20 10 100 1 9 0.5 > tmp-9f &
fi

# you can run the R script on the fly as results are being produced, but it obviously won't provide a full analysis until the code is finished running
