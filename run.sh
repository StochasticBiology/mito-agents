# compile source code
gcc -o3 mito-agents.c -lm -o mito-agents.ce

# run different experiments in parallel
./mito-agents.ce 20 0 > tmp-0a &
./mito-agents.ce 20 1 > tmp-1a &
./mito-agents.ce 20 2 > tmp-2a &
./mito-agents.ce 20 3 > tmp-3a &
./mito-agents.ce 20 4 1e-5 > tmp-4a &
./mito-agents.ce 20 5 1e-5 > tmp-5a &
./mito-agents.ce 20 6 1 > tmp-6a &
./mito-agents.ce 20 7 1 > tmp-7a &

./mito-agents.ce 50 0 > tmp-0 &
./mito-agents.ce 50 1 > tmp-1 &
./mito-agents.ce 50 2 > tmp-2 &
./mito-agents.ce 50 3 > tmp-3 &
./mito-agents.ce 50 4 1e-5 > tmp-4 &
./mito-agents.ce 50 5 1e-5 > tmp-5 &
./mito-agents.ce 50 6 1 > tmp-6 &
./mito-agents.ce 50 7 1 > tmp-7 &

# you can run the R script on the fly as results are being produced, but it obviously won't provide a full analysis until the code is finished running
