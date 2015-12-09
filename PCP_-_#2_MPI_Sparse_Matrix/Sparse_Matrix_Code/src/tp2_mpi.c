#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>
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

int main(int argc, char *argv[]) {		
	int nLinhas, nCols, nLinhasVect, i, j, incr;
	double startTime, finalTime, n;
	double **coo, *vect, *result, *linha;
	int MASTER = 0;

	/* MPI */
	MPI_Status status;	
	int idProc;		/* ID de cada processo */
	int tag;
	int totalProcs;		/* Numero total de processos */
	int nrLinha;		/* Numero de linha a ser processada */
	int procDest;		/* Processo que vai tratar a linha corrente */
	int acabou; 		/* Para terminar a execucao do processo */
	double msg; 		/* Para a comunicação entre o MASTER e os outros processos */

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &idProc);
	MPI_Comm_size(MPI_COMM_WORLD, &totalProcs);

	nLinhas = 5;	//atoi(argv[1]);
	nCols = 5; 	//atoi(argv[2]);
	nLinhasVect = 5;	//atoi(argv[3]);
	acabou = nLinhas + 1;
	incr = (int)(nLinhas/(totalProcs-1));

	printf("Incremento %d\n", incr);
	if (idProc == MASTER){
		/* Preenchimento Matriz */
		unsigned int nnz = (nLinhas * nCols) * 0.25;
		coo = geraMatriz(nLinhas, nCols, nnz);

		/* Preenchimento Vector */
		vect = geraVetor(nLinhasVect);

		result = (double*) calloc (nLinhas, sizeof(double));
		linha = (double*) calloc (150, sizeof(double));

		/*
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
		printf("\n");*/

	} else {
		vect = (double*) calloc (nLinhasVect, sizeof(double));
		linha = (double*) calloc (150, sizeof(double));
	}

	/* Enviar o vetor aos outros processos */
	MPI_Bcast (vect, nLinhasVect, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	
	if (idProc == MASTER){
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
		while(i < nLinhas){
			tag = i;
			MPI_Send(coo[i], (int)(coo[i][0] + 1), MPI_DOUBLE, procDest, tag, MPI_COMM_WORLD);
			i++;
		}

		for(i=1; i<totalProcs; i++) {
			tag = acabou;
			linha[0] = -1;
			MPI_Send(&linha,1, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
		}

		nrLinha = 0;
		while(1){
			MPI_Recv(&msg, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			tag = status.MPI_TAG;
			result[tag] = msg;
			printf("Linha %d resultado %f\n", tag, msg);
			nrLinha++;
			if ( nrLinha == nLinhas){ break; }
		}

		/* Libertar a memoria alocada */
		for(i=0; i<nLinhas; i++) {
			free(coo[i]);
		}
		free(coo);
		free(linha);
		free(vect);
	}

	//Caso nao seja o MASTER
	else{
		double res;
		MPI_Recv(linha, 150, MPI_DOUBLE, MASTER, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		tag = status.MPI_TAG;
		while(tag != acabou){
			res = 0.0;
			n = (int)linha[0];

			for(i=1; i<=n; i+=2){
				printf("(%f * %f) + ", linha[i+1], vect[(int)linha[i]]);
				res += linha[i+1] * vect[(int)linha[i]];
			}
			printf(" = %f\n", res);
			MPI_Send(&res, 1, MPI_DOUBLE, MASTER, tag, MPI_COMM_WORLD);
			MPI_Recv(linha, 150, MPI_DOUBLE, MASTER, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			tag = status.MPI_TAG;
		}
		
		free(linha);
		free(vect);
	}
	
	if (idProc == MASTER){
		//printf("Time seq: %.12f\n",(finalTime - startTime)*1000);
		//printf("PAI acabou\n");
		
		/* Visualizacao dos resultados */
		printf("-----------Resultado Multiplicacao---------\n");
		for(i=0; i<nLinhas; i++)
			printf("---%f ", result[i]);
		
		printf("\n");

		free(result);
	}

	MPI_Finalize();
	return 0;
}

