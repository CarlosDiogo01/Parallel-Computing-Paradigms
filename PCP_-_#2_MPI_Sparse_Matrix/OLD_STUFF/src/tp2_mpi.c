#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

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

void forceLoadram(){
	int i;
  	double forceloadram [30000];
  	for (i = 0; i < 30000; ++i)
    	forceloadram[i] = i;
}

int main(int argc, char *argv[]) {		
	int nLinhas, nCols, nLinhasVect, i, j, incr, maxElem;
	double startTime, endTime, startAux, endAux, tcomp, trecv, tsend, wallTime, n;
	double **coo, *vect, *result, *linha, res;
	int MASTER = 0;

	/* MPI */
	MPI_Status status;	
	int idProc;			/* ID de cada processo */
	int tag;
	int totalProcs;		/* Numero total de processos */
	int procDest;		/* Processo que vai tratar a linha corrente */
	int acabou; 		/* Para terminar a execucao do processo */

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &idProc);
	MPI_Comm_size(MPI_COMM_WORLD, &totalProcs);

	nLinhas = atoi(argv[1]);
	nCols = atoi(argv[2]);
	nLinhasVect = atoi(argv[3]);
	acabou = nLinhas + 1;
	incr = (int)(nLinhas/(totalProcs-1));
	maxElem = (nCols*2) + 1;
	tcomp = 0; tsend = 0; trecv = 0; 
	wallTime = 0; startTime = endTime = 0;

	linha = (double*) calloc (maxElem, sizeof(double));
	result = (double*) calloc (nLinhas, sizeof(double));
	vect = (double*) calloc (nLinhasVect, sizeof(double));

	if (idProc == MASTER){
		/* Preenchimento Matriz e Vetor */
		unsigned int nnz = (nLinhas * nCols) * 0.25;
		coo = geraMatriz(nLinhas, nCols, nnz);
		vect = geraVetor(nLinhasVect);
	}

	/* Enviar o vetor gerado no MASTER aos outros processos */
	MPI_Bcast(vect, nLinhasVect, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	if (idProc == MASTER){
		forceLoadram();
		startTime = MPI_Wtime();
		procDest = 1;
		for (i=0; (i<nLinhas) && (procDest<totalProcs); i+= incr){
			n = i+incr;
			for (j=i; j<n && j<nLinhas; j++) {
				tag = j;
				MPI_Send(coo[j], (int)(coo[j][0] + 1), MPI_DOUBLE, procDest, tag, MPI_COMM_WORLD);
			}
			procDest++;
		}
		
		procDest--;
		/* Caso nao sejam enviadas todas as linhas no ciclo anterior */
		while(i < nLinhas){
			tag = i;
			MPI_Send(coo[i], (int)(coo[i][0] + 1), MPI_DOUBLE, procDest, tag, MPI_COMM_WORLD);
			i++;
		}

		/* Enviar msg para os processos terminarem */
		for(i=1; i<totalProcs; i++) {
			tag = acabou;
			linha[0] = -1;
			MPI_Send(&linha,1, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
		}

		endTime = MPI_Wtime();
		wallTime = endTime - startTime;
		tsend = wallTime;
		//printf("Wall time proc %d: %.12f\n", idProc, wallTime*100);
	}

	//Caso nao seja o MASTER
	else{
		startTime = startAux = MPI_Wtime();
		MPI_Recv(linha, maxElem, MPI_DOUBLE, MASTER, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		endAux = MPI_Wtime();
		trecv += endAux - startAux;

		tag = status.MPI_TAG;
		while(tag != acabou){
			startAux = MPI_Wtime();
			res = 0.0;
			n = (int)linha[0];
			for(i=1; i<=n; i+=2){
				res += linha[i+1] * vect[(int)linha[i]];
			}
			result[tag] = res;
			endAux = MPI_Wtime();
			tcomp += endAux - startAux;

			startAux = MPI_Wtime();
			MPI_Recv(linha, maxElem, MPI_DOUBLE, MASTER, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			endAux = MPI_Wtime();
			trecv += endAux - startAux;
			tag = status.MPI_TAG;
		}
		endTime = MPI_Wtime();
		wallTime = endTime - startTime;
		//printf("Wall time proc %d: %.12f\n", idProc, wallTime*1000);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	double result_global[nLinhas], wallTime_global, trecv_global, tcomp_global;
	
	MPI_Reduce(result, result_global, nLinhas, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
	MPI_Reduce(&wallTime, &wallTime_global, 1, MPI_DOUBLE, MPI_MAX, MASTER, MPI_COMM_WORLD);
	MPI_Reduce(&trecv, &trecv_global, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
	MPI_Reduce(&tcomp, &tcomp_global, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	if(idProc == MASTER) {
		/* Visualizacao dos resultados */
		/*printf("---------Matriz formato COO---------\n");
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
		*/
		
		printf("%.12f, %.12f, %.12f, %.12f\n", tsend*1000, trecv_global*1000, tcomp_global*1000, wallTime_global*1000);
		
		/* Libertar a memoria alocada */
		for(i=0; i<nLinhas; i++) {
			free(coo[i]);
		}
		free(coo);
	}

	free(linha);
	free(vect);
	free(result);

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}

