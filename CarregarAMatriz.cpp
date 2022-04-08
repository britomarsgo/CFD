#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cstdlib>
#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <string>
#include <sys/time.h>
#include <string.h>


double** carregaMatriz(char* nomeArquivo, int IMAX, int JMAX){

// declara uma matriz de ponteiros
double **m;

// declara uma variavel file handle para um arquivo de texto
std::ifstream txtFile;

int l=IMAX,c=JMAX;

//Abrindo o arquivo [matriz.dat]
//txtFile.open("C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto3/Matriz.txt");

  char *str; /* o ponteiro para o espaço alocado */
  /* aloco um único byte na memória */
  str = (char *)malloc(100);
  strcpy (str,"C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/");
  strcat (str,nomeArquivo);
  strcat(str,".txt");

  txtFile.open(str);


m = new double* [l];
    for (int i = 0; i < l; i++){
        m[i] = new double[c];
    }

for (int i=0; i<l;i++){

    for (int j=0; j<c;j++){
    txtFile >> m[i][j];

    }

}

return m;

}

