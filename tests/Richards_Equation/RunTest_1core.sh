# Number of OpenMP threads
export OMP_NUM_THREADS=1

# parallel = 0 (Serial), parallel = 1 (OpenMP)
export parallel=1

cp '/Users/elliotcarr/Dropbox/Documents/Research/DECRA/Code/twoscale_dev/twoscale/bin/run' .

./run problem_options.txt