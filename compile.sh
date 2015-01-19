clear

# Delete binary
if [ -a bin/run ]; then

rm bin/run
echo Deleted old binary.

fi

# Compile the code using g++
echo Compiling the code. 

mkdir bin

g++ ./src/*.cpp -o ./bin/run -O2 -framework Accelerate -fopenmp -I include

#g++ main.cpp macro_mesh_properties.cpp micro_mesh_properties.cpp \
#read_effective_conductivity.cpp exprem.cpp phipade.cpp set_exprem_options.cpp \
#Gfunc_serial.cpp Gfunc_parallel.cpp save_solution.cpp -o ./bin/run -O2 \
#-framework Accelerate -fopenmp -I include

# Catch compile error otherwise run the code
if [ $? -eq 0 ]; then

echo Code compiled successfully.

# Number of OpenMP threads
export OMP_NUM_THREADS=1

# parallel = 0 (Serial), parallel = 1 (OpenMP)
export parallel=1

else

echo Code failed to compile.

fi