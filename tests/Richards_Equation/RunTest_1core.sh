# Number of OpenMP threads
export OMP_NUM_THREADS=1

# parallel = 0 (Serial), parallel = 1 (OpenMP)
export parallel=1

cd ..
cd ..
cp bin/run tests/Richards_Equation
echo Exectuable copied from bin to current directory
cd tests/Richards_Equation

./run problem_options2.txt