#!/bin/sh


#PBS -N tp2_hybrid
#PBS -l walltime=02:00:00
#PBS -q mei

#PBS -m abe
#PBS -M carlos.sa01@gmail.com
#PBS -M a69855@alunos.uminho.pt


#PBS -l nodes=2:r641:ppn=32

cd ~/PCP/Code_With_Hybrid_Solution/
module purge
module load gnu/4.9.0
module load gnu/4.9.3
module load gnu/openmpi_eth/1.8.4
ompi_info --param mpi all

read -r node_info<$PBS_NODEFILE

echo "Allocated computing node: $node_info"
echo "Eth Network"

max_ppn=64
sample_size=5
max_matrix=4096
max_omp_threads=32

for (( matrix_size=1024, vec_size=1024; matrix_size <= $max_matrix; matrix_size+=matrix_size, vec_size+=vec_size ))
	do
		echo "Running for COO Matrix Size: $matrix_size x $matrix_size"
		echo "Running for Vector Size : $vec_size"
		for (( ppn=2; ppn <= max_ppn; ppn+=2 )) 
		do
			echo "Running $sample_size * ppn: $ppn processes"
			for (( omp_threads=2; omp_threads <= max_omp_threads; omp_threads+=2 ))
			do
				echo "Running each process with: $omp_threads threads OpenMP"
				for (( seq_num=1; seq_num <= $sample_size; ++seq_num ))
				do
					mpirun -np $ppn --map-by core -mca btl self,sm,tcp --report-bindings bin/tp2_hybrid $matrix_size $matrix_size $vec_size $omp_threads "--map-by core" "eth" 2
				done
			done
		done
	done
echo "Done..."

