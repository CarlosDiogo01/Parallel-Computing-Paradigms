#/bin/bash
COMPILER=mpicc
PROCESSES_NUMBER=2
MPI_FILE=tp2_hybrid

echo "Start compiling ...."
make all
cd bin

echo "EXECUTING MPIRUN"
mpirun -np $PROCESSES_NUMBER $MPI_FILE
cd ..
make clean

