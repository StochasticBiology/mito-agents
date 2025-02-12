# compile source code
gcc -o3 mito-agents.c -lm -o mito-agents.ce

# run different experiments in parallel
./mito-agents.ce 20 2 0 > tmp-0a &
./mito-agents.ce 20 2 1 > tmp-1a &
./mito-agents.ce 20 2 2 > tmp-2a &
./mito-agents.ce 20 2 3 > tmp-3a &
./mito-agents.ce 20 2 4 5e-6 > tmp-4a &
./mito-agents.ce 20 2 5 5e-6 > tmp-5a &
./mito-agents.ce 20 2 6 1 > tmp-6a &
./mito-agents.ce 20 2 7 1 > tmp-7a &

./mito-agents.ce 50 2 0 > tmp-0 &
./mito-agents.ce 50 2 1 > tmp-1 &
./mito-agents.ce 50 2 2 > tmp-2 &
./mito-agents.ce 50 2 3 > tmp-3 &
./mito-agents.ce 50 2 4 5e-6 > tmp-4 &
./mito-agents.ce 50 2 5 5e-6 > tmp-5 &
./mito-agents.ce 50 2 6 1 > tmp-6 &
./mito-agents.ce 50 2 7 1 > tmp-7 &

./mito-agents.ce 20 1 0 > tmp-0c &
./mito-agents.ce 20 1 1 > tmp-1c &
./mito-agents.ce 20 1 2 > tmp-2c &
./mito-agents.ce 20 1 3 > tmp-3c &
./mito-agents.ce 20 1 4 5e-6 > tmp-4c &
./mito-agents.ce 20 1 5 5e-6 > tmp-5c &
./mito-agents.ce 20 1 6 1 > tmp-6c &
./mito-agents.ce 20 1 7 1 > tmp-7c &

# you can run the R script on the fly as results are being produced, but it obviously won't provide a full analysis until the code is finished running
