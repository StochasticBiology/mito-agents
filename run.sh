# compile source code
gcc -o3 mito-agents.c -lm -o mito-agents.ce

# run different experiments in parallel
./mito-agents.ce 0 > tmp-0 &
./mito-agents.ce 1 > tmp-1 &
./mito-agents.ce 2 > tmp-2 &
./mito-agents.ce 3 > tmp-3 &

# you can run the R script on the fly as results are being produced, but it obviously won't provide a full analysis until the code is finished running
