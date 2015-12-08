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
	unsigned int i, j, n, nnz, nLinhas, nCols, nLinhasVect;
	double startTime, finalTime;
	double **coo, *vect, *result;
	double soma;

	/* MPI */
	MPI_Status status;	
	int idProc;		/* ID de cada processo */
	int tag;
	double res;		/* Resultado de fazer a multiplicacao de uma linha pela coluna do vector */
	int totalProcs;		/* Numero total de processos */
	int nrLinha;		/* Numero de linha a ser processada */
	int procDest;		/* Processo que vai tratar a linha corrente */
	int nrProc;		/* Numero de processos correntes */
	int acabou; 		/* Para terminar a execucao do processo */
	double msg; 		/* Para a comunicação entre o MASTER e os outros processos */
	double *linha;
	int MASTER = 0;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &idProc);
	MPI_Comm_size(MPI_COMM_WORLD, &totalProcs);

	nLinhas = 50;	//atoi(argv[1]);
	nCols = 50; 	//atoi(argv[2]);
	nLinhasVect = 50;	//atoi(argv[3]);
	acabou = nLinhas + 1;

	if (idProc == MASTER){
		/* Preenchimento Matriz */
		nnz = (nLinhas * nCols) * 0.25;
		coo = geraMatriz(nLinhas, nCols, nnz);

		/* Preenchimento Vector */
		vect = geraVetor(nLinhasVect);

		result = (double*) calloc (nLinhas, sizeof(double));
		linha = (double*) calloc (150, sizeof(double));

		/*printf("---------Matriz formato COO---------\n");
		for(i=0; i<nLinhas; i++){
			if (coo[i]) {
				n = coo[i][0];
				for(j=1; j<=n; j+=2) {
					printf("%f %f %f\n", (double)i, (double)coo[i][j], (double)coo[i][j+1]);
				}
			}
		}*/
		
		printf("\n---------Vetor---------\n");
		for(i=0; i<nLinhasVect; i++) {
			printf("%f ", vect[i]);
		}

	} else {
		vect = (double*) calloc (nLinhasVect, sizeof(double));
		linha = (double*) calloc (150, sizeof(double));
	}

	/* Enviar o vetor aos outros processos */
	MPI_Bcast (vect, nLinhasVect, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	
	if (idProc == MASTER){
		nrLinha = 0;
		for (i=1; i<= totalProcs-1; i++){
			procDest = i;
			tag = nrLinha;
			MPI_Send(coo[nrLinha], (int)coo[nrLinha][0], MPI_DOUBLE, procDest, tag, MPI_COMM_WORLD);
			nrLinha++;
		}

		nrProc = totalProcs-1;
		while(1){
			MPI_Recv(&msg, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_SOURCE, MPI_COMM_WORLD, &status);
			tag = status.MPI_TAG;
			result[tag] = msg;
			printf("PAI: recebeu do processo %i res %d linha %i\n", status.MPI_SOURCE, result[tag], tag);
			if ( nrLinha < nLinhas){
				tag = nrLinha;
				procDest = status.MPI_SOURCE;
				MPI_Send(coo[nrLinha], (int)coo[nrLinha][0], MPI_DOUBLE, procDest, tag, MPI_COMM_WORLD);	//ana nao tem a certeza
				nrLinha++;
			}
			else{
				nrProc--;
				tag = acabou;
				procDest = status.MPI_SOURCE;
				linha[0] = -1.0;
				MPI_Send(&linha, 1, MPI_DOUBLE, procDest, tag, MPI_COMM_WORLD);
				printf("PAI: Fez Send ACABOU!\n");
				if (nrProc == 0) break;
			}
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
		MPI_Recv(linha, 150, MPI_DOUBLE, MASTER, MPI_ANY_SOURCE, MPI_COMM_WORLD, &status);
		tag = status.MPI_TAG;
		while(tag != acabou){
			res = 0.0;
			for(i=1; i<=(int)linha[0]; i+=2){
				res += linha[i+1] * vect[(int)linha[i]]	;
			}
			//printf("FILHO %i: enviou res %d da LINHA MAT %i-------\n", idProc, res, tag);
			MPI_Send(&res, 1, MPI_DOUBLE, MASTER, tag, MPI_COMM_WORLD);
			MPI_Recv(linha, 150, MPI_DOUBLE, MASTER, MPI_ANY_SOURCE, MPI_COMM_WORLD, &status);
			tag = status.MPI_TAG;
		}
		printf("FILHO %i ACABOU\n", idProc);
		free(linha);
		free(vect);
	}

	/* Multiplicacao */
	//startTime = omp_get_wtime();
	/*for(i=0; i<nLinhas; i++) {
		n = (int)coo[i][0];
		soma=0;
		for(j=1; j<=n; j+=2){
			soma += coo[i][j+1]*vect[0][(int)coo[i][j]];
		}
		result[i] = soma;
	}*/
	//finalTime = omp_get_wtime();

	
	if (idProc == MASTER){
		//printf("Time seq: %.12f\n",(finalTime - startTime)*1000);
		printf("PAI acabou\n");
		
		/* Visualizacao dos resultados */
		/*	
		printf("-----------Resultado Multiplicacao---------\n");
		for(i=0; i<nLinhas; i++)
			printf("%f\n", result[i]);
		*/

		free(result);
	}

	MPI_Finalize();
	return 0;
}

