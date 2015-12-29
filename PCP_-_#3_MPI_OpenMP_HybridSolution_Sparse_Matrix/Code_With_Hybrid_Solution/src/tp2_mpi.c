#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include <omp.h>

/********************
Desenvolvido por: 
Carlos Sá- A59905
Ana Sousa - A69855
*********************/


/****** Generating COO Matrix *******/
double** geraMatriz (int nLinhas, int nCols, int nnz) {
	int i, j, x, y, aux, cont;
	double **mat = (double **) calloc (nLinhas, sizeof(double*));
	double **coo = (double **) calloc (nLinhas, sizeof(double*));
	double valRand;

	for (i=0; i<nLinhas; i++) {
		mat[i] = (double*) calloc (nCols, sizeof(double));
	}

	srand(1);
	for (i=0; i<nnz; i += aux) {
		x = rand() % nLinhas;
		y = rand() % nCols;
		if (mat[x][y] == 0.0) {
			while ((valRand = (rand() % 50) * 0.25) == 0) {;}
			mat[x][y] = valRand; 
			aux = 1;
		} else {
			aux = 0;
		}
	}

	for (i=0; i<nLinhas; i++) {
		cont = 0;

		for(j=0; j<nCols; j++) { if (mat[i][j] != 0.0) {cont++;} }
		if (cont>0) {
			coo[i] = (double*) calloc ((cont*2)+1, sizeof(double));
			coo[i][0] = cont * 2;
			aux = 1;
			for (j=0; j<nCols; j++) {
				if (mat[i][j] != 0.0) {
					coo[i][aux] = j;
					coo[i][aux+1] = mat[i][j];
					aux += 2;

				}
			}
		} else {
			coo[i] = (double*) calloc (1, sizeof(double));
			coo[i][0] = 0;
		}
		free(mat[i]);
	}
	free(mat);
	
	return coo;
}

/***** Generating the Vector *****/
double* geraVetor (int nCols) {
	int i;
	double *vect;

	vect = (double*) calloc (nCols, sizeof(double));
	srand(1);
	for(i=0; i<nCols; i++) {
		vect[i] = (rand() % 20) * 0.25;
	}
	
	return vect;
}

void forceLoadRam(){
	int i;
  	double forceloadram [30000000];
  	for (i = 0; i < 30000000; ++i)
    	forceloadram[i] = i;
}

int main(int argc, char *argv[]) {		
	int nLinhas, nCols, nLinhasVect, i, j, incr, maxElem;
	double startTime, endTime, startAux, endAux, tbcast, tcomp, trecv, tsend, wallTime, n;
	double **coo, *vect, *result, *linha, res;
	int MASTER = 0;

	/*** OpenMP Addictions for Hybrid Solution ***/
	int total_threads_per_process = atoi(argv[4]);		/** Receiving the Number of OpenMP Threads for each MPI Process */



	/* MPI */
	MPI_Status status;	
	int idProc;		/* ID de cada processo */
	int tag;
	int totalProcs;		/* Numero total de processos */
	int procDest;		/* Processo que vai tratar a linha corrente */

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &idProc);
	MPI_Comm_size(MPI_COMM_WORLD, &totalProcs);

	nLinhas = atoi(argv[1]);
	nCols = atoi(argv[2]);
	nLinhasVect = atoi(argv[3]);
	incr = (int)(nLinhas/(totalProcs-1));
	maxElem = (nCols*2) + 1;
	tcomp = 0; tsend = 0; trecv = 0; tbcast = 0;
	wallTime = 0; startTime = endTime = 0;

	linha = (double*) calloc (maxElem, sizeof(double));
	result = (double*) calloc (nLinhas, sizeof(double));
	vect = (double*) calloc (nLinhasVect, sizeof(double));


	/***** Filling the COO Matrix and Vector ******/
	if (idProc == MASTER){
		/* Preenchimento Matriz e Vetor */
		unsigned int nnz = (nLinhas * nCols) * 0.25;
		coo = geraMatriz(nLinhas, nCols, nnz);
		vect = geraVetor(nLinhasVect);
	}

	/***** Sending generated Vector to other processes *****/
	startAux = MPI_Wtime();
	MPI_Bcast(vect, nLinhasVect, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	endAux = MPI_Wtime();
	tbcast = endAux - startAux;
//	forceLoadRam();
	MPI_Barrier(MPI_COMM_WORLD);

	if (idProc == MASTER){
		/** START MEASURING MPI AND OPENMP WORK **/
		startTime = MPI_Wtime(); //walltime MPI
		startAux = MPI_Wtime(); //tsend MPI
		procDest = 1;
		int val_n;		//used to "bypass" n value. OpenMP does not accept variable n in for loop condition
		for (i=0; ( (i<nLinhas) && (procDest < totalProcs) ) ; i+= incr){
			n = i+incr;
			val_n = n;
			//#pragma omp parallel num_threads(total_threads_per_process)
			//{
			//#pragma omp for nowait
			for (j=i; j<val_n; j++) {
				tag = j;
				MPI_Send(coo[j], (int)(coo[j][0] + 1), MPI_DOUBLE, procDest, tag, MPI_COMM_WORLD);
			}
			//}
			procDest++;
		}	
		endAux = MPI_Wtime(); 
		tsend += (endAux - startAux);
		
		//Computacao das linhas que sobraram
		startAux = MPI_Wtime(); //tcomp
		int line_left_i;		//used to "bypass" i. OpenMP needs a initialization variable in for loop
		int val;			//used to "bypass" n value. OpenMP does not accept variable n in for loop condition
		for(line_left_i=i; line_left_i<nLinhas; line_left_i++){	
			res = 0.0;
			n = (int)coo[line_left_i][0];
			val = n;
			#pragma omp parallel num_threads(total_threads_per_process)
			{
			#pragma omp for reduction(+:res)
			for(j=1; j<val; j+=2){
				res += coo[line_left_i][j+1] * vect[(int)coo[line_left_i][j]];
			}
			}
			result[line_left_i] = res;
		}
		endAux = MPI_Wtime(); 
		tcomp = endAux - startAux;
		
		/** Walltime MPI **/
		endTime = MPI_Wtime();
		wallTime = endTime - startTime;
	}

	//Caso nao seja o MASTER
	else{
		/** Start measuring MPI times **/
		startTime = MPI_Wtime();
		for(i=0; i<incr; i++){
			startAux = MPI_Wtime();
			MPI_Recv(linha, maxElem, MPI_DOUBLE, MASTER, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			tag = status.MPI_TAG;
			endAux = MPI_Wtime();
			trecv += (endAux - startAux);
			startAux = MPI_Wtime();
			res = 0.0;
			int valn;
			n = (int)linha[0];
			valn = n;
			#pragma omp parallel num_threads(total_threads_per_process)
			{
			#pragma omp for schedule(dynamic) reduction(+:res)
			for(j=1; j <= valn; j+=2){
				res += linha[j+1] * vect[(int)linha[j]];
			}
			result[tag] = res;
			endAux = MPI_Wtime();
			tcomp += (endAux - startAux);
			}
		}
		/** Walltime MPI **/
		endTime = MPI_Wtime();
		wallTime = endTime - startTime;
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	double result_global[nLinhas], wallTime_global, trecv_global, tcomp_global;
	
	MPI_Reduce(result, result_global, nLinhas, MPI_DOUBLE, MPI_MAX, MASTER, MPI_COMM_WORLD);
	MPI_Reduce(&wallTime, &wallTime_global, 1, MPI_DOUBLE, MPI_MAX, MASTER, MPI_COMM_WORLD);
	MPI_Reduce(&trecv, &trecv_global, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
	MPI_Reduce(&tcomp, &tcomp_global, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	if(idProc == MASTER) {
		/* Visualizacao dos resultados
		printf("---------Matriz formato COO---------\n");
		for(i=0; i<nLinhas; i++){
			if (coo[i]) {
				n = coo[i][0];
				for(j=1; j<=n; j+=2) {
					printf("%f %f %f\n", (double)i, (double)coo[i][j], (double)coo[i][j+1]);
				}
			}
		}
		
		printf("\n---------Vetor---------\n");
		for(i=0; i<nLinhasVect; i++) {
			printf("---%f ", vect[i]);
		}
		printf("\n");

		printf("-----------Resultado Multiplicacao---------\n");
		for(i=0; i<nLinhas; i++)
			printf("---%f ", result_global[i]);
		printf("\n\n");
		
		**/
		printf("%.12f, %.12f, %.12f, %.12f, %.12f\n", tbcast*1000, tsend*1000, trecv_global*1000, tcomp_global*1000, wallTime_global*1000);

		/* Libertar a memoria alocada */
		 for(i=0; i<nLinhas; i++) {
			free(coo[i]);
		}
		free(coo);
	}

	free(vect);
	free(linha);
	free(result);

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

return 0;
}
