/*
CC-297 - Projeto 3 - Solução do Potencial Trânsonico Completo
Programa Desenvolvido por : Alisson Vinicius Brito Lopes
*/

// Include das principais bibliotecas do C/C++
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cstdlib>
#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <string>
#include <sys/time.h>
#include "CarregarAMatriz.h"

# define M_PI 3.14159265358979323846  /* pi */

using namespace std;

class WriteFile{

public:

static void writeMatrix(double **M, int qtdlinhas, int qtdColunas, char nomeArquivo[]){

     ofstream write;
     write.open(nomeArquivo);

     for (int i = 0; i < qtdlinhas; i++ ){

        write << "\n";

        for (int j = 0; j < qtdColunas; j++){

             write.width(10);
             write << M[i][j] << "        ";

        }
    }
             write.close();
    }


static void writeVetAgrup(double *V1, double *V2, int qtdlinhas, char nomeArquivo[]){

    double MatAgrVet[qtdlinhas][2];

    for(int i=0; i<qtdlinhas;i++){
        MatAgrVet[i][0]=V1[i];
    }

    for(int i=0; i<qtdlinhas;i++){
        MatAgrVet[i][1]=V2[i];
    }

     ofstream write;
     write.open(nomeArquivo);
       //variables="x" "y"
  //zone I=          93  ,J=          15  ,K=           1  F=point

  write << " variables = ''x'' ''y'' " << "\n";
  write << "zone I=          93  ,J=          15  ,K=           1  F=point"<<endl;

     for (int i = 0; i < qtdlinhas; i++ ){

        write << "\n";

        for (int j = 0; j < 2; j++){

             write.width(11);
             write << MatAgrVet[i][j] << " ";

        }
    }
             write.close();

}


static void writeVetAgrup1(double *V1, double *V2, double *V3, int qtdlinhas, char nomeArquivo[]){

    double MatAgrVet[qtdlinhas][3]; // alterar a quantidade de colunas para salvar

    for(int i=0; i<qtdlinhas;i++){
        MatAgrVet[i][0]=V1[i];
    }

    for(int i=0; i<qtdlinhas;i++){
        MatAgrVet[i][1]=V2[i];
    }

    for(int i=0; i<qtdlinhas;i++){
        MatAgrVet[i][2]=V3[i];

        //cout << MatAgrVet[i][1] << endl;
    }

     ofstream write;
     write.open(nomeArquivo);
       //variables="x" "y" "Phi"
  //zone I=          93  ,J=          15  ,K=           1  F=point

  write << " variables = ''x'' ''y'' ''Mach'' " << "\n";
  write << "zone I=          93  ,J=          15  ,K=           1  F=point"<<endl;

     for (int i = 0; i < qtdlinhas; i++ ){

        write << "\n";

        for (int j = 0; j < 3; j++){

             write.width(11);
             write << MatAgrVet[i][j] << " ";

        }
    }
             write.close();

}


static void writeVetor(double *V, int qtdColunas, char nomeArquivo[]){

     ofstream write;
     write.open(nomeArquivo);

     for (int i = 0; i < qtdColunas; i++ ){

            write << V[i] <<"\n";

        }
             write.close();
}

};

class Vetor{

    public:

    static void imprimir(double *V, int qtdColunas){

    for (int i = 0; i < qtdColunas; i++ ){

        cout << V[i] << "\t";
    }
        cout << endl;
    }


};

class Matrix{

public:

   static void  imprimir(double **M,int qtdLinhas, int qtdColunas){ // Nome de objeto minusculo

        for ( int i = 0; i < qtdLinhas; i++){

                Vetor::imprimir(M[i],qtdColunas);
        }
    }

    static void imprimirTecPlot(double **M,int qtdLinhas, int qtdColunas){

        for (int j = 0; j < qtdColunas; j++ ){
            for (int i = 0; i < qtdLinhas; i++){
                    cout << M[i][j] << "\n";
            }
         }
    }

};

class MatrixTri{

public:

   static void  imprimir(double **M,int qtdLinhas, int qtdColunas, int nDim ){ // Nome de objeto minusculo

        for ( int i = 0; i < qtdLinhas; i++){

                Vetor::imprimir(M[i],qtdColunas);
        }
    }

    static void imprimirTecPlot(double **M,int qtdLinhas, int qtdColunas){

        for (int j = 0; j < qtdColunas; j++ ){
            for (int i = 0; i < qtdLinhas; i++){
                    cout << M[i][j] << "\n";
            }
         }
    }

};

class PotFlow02 {

    public:

    int IMAX, JMAX, ILE, ITE;
    double XSF, YSF;

// Criando ponteiro duplos e simples para Alocação dinâmica de memoria e "zeramento" de vetores e matrizes na memoria do computador

    double **f, **g,**DELTAPhi,**Phi,**LPhi, **Cij;

    double *RHS, *AT, *BT, *CT, *DT, *ATX, *BTX, *CTX, *VT,*ATY, *BTY, *CTY;

    double **A, **B, **C, **D,**AMX, **AMY;

    double **x,**y;

    double *Rx,*Ry;

    double ***rho,***rhobarra,***rhotil,***A1, ***A2, ***A3, ***J, ***U, ***V, ***Phiqsi, ***Phieta, ***Bbarra, ***Btil ;

    double **etax,**qsix,**etay, **qsiy, **v, **u, **Ucampo,**rhocampo, **cp, **MACH,*cplinha;


void inicializaMatriz_Ucampo(int IMAX, int JMAX){ // MATRIX DE RESIDUOS

  int i, j;

  Ucampo = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        Ucampo[i] = new double[JMAX];
    }

    // Inicializando a matriz com as condições de contorno de : Entrada, superior e saída

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 Ucampo[i][j] = 0.0;
         }
    }
}

void inicializaMatriz_rhocampo(int IMAX, int JMAX){ // MATRIX DE RESIDUOS

  int i, j;

  rhocampo = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        rhocampo[i] = new double[JMAX];
    }

    // Inicializando a matriz com as condições de contorno de : Entrada, superior e saída

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 rhocampo[i][j] = 0.0;
         }
    }
}

void inicializaMatriz_Bbarra(int IMAX, int JMAX, int nDim){

  int i, j, k;

 Bbarra = new double**[IMAX];

  for (int i = 0; i<IMAX; i++) {
      Bbarra[i] = new double*[JMAX];
      for (int j = 0; j<JMAX; j++) {
          Bbarra[i][j] = new double[nDim];
          for (int k = 0; k<nDim; k++) {
               Bbarra[i][j][k] = 0.0;
          }
       }
  }

}


void inicializaMatriz_Btil(int IMAX, int JMAX, int nDim){

  int i, j, k;

 Btil = new double**[IMAX];

  for (int i = 0; i<IMAX; i++) {
      Btil[i] = new double*[JMAX];
      for (int j = 0; j<JMAX; j++) {
          Btil[i][j] = new double[nDim];
          for (int k = 0; k<nDim; k++) {
               Btil[i][j][k] = 0.0;
          }
       }
  }

}

void inicializaMatriz_Phi(int IMAX, int JMAX,double Uinf,double gamma, double theta){ // Phi[i][0] = Uinf*(x[i][j]*cos(theta) + y[i][j]*sin(theta))

  int i, j;

  Phi = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        Phi [i] = new double[JMAX];
    }

    // Inicializando a matriz com as condições de contorno de : Entrada, superior e saída

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 Phi [i][j] = 0.0;

         }
    }
}

void inicializaMatriz_Phiqsi(int IMAX, int JMAX, int nDim,double Uinf,double gamma){

  int i, j, k;

 Phiqsi= new double**[IMAX];

  for (int i = 0; i<IMAX; i++) {
      Phiqsi[i] = new double*[JMAX];
      for (int j = 0; j<JMAX; j++) {
          Phiqsi[i][j] = new double[nDim];
          for (int k = 0; k<nDim; k++) {
               Phiqsi[i][j][k] = 0.0;
          }
       }
  }

}

void inicializaMatriz_Phieta(int IMAX, int JMAX, int nDim,double Uinf,double gamma){

  int i, j, k;

 Phieta= new double**[IMAX];

  for (int i = 0; i<IMAX; i++) {
      Phieta[i] = new double*[JMAX];
      for (int j = 0; j<JMAX; j++) {
          Phieta[i][j] = new double[nDim];
          for (int k = 0; k<nDim; k++) {
               Phieta[i][j][k] = 0.0;
          }
       }
  }

}


void inicializaMatriz_rho(int IMAX, int JMAX, int nDim){

  int i, j, k;

 rho = new double**[IMAX];

  for (int i = 0; i<IMAX; i++) {
      rho[i] = new double*[JMAX];
      for (int j = 0; j<JMAX; j++) {
          rho[i][j] = new double[nDim];
          for (int k = 0; k<nDim; k++) {
                rho[i][j][k] = 0.0;
          }
       }
  }

}

void inicializaMatriz_rhobarra(int IMAX, int JMAX, int nDim){

  int i, j, k;

 rhobarra = new double**[IMAX];

  for (int i = 0; i<IMAX; i++) {
      rhobarra[i] = new double*[JMAX];
      for (int j = 0; j<JMAX; j++) {
          rhobarra[i][j] = new double[nDim];
          for (int k = 0; k<nDim; k++) {
                rhobarra[i][j][k] = 0.0;
          }
       }
  }

}

void inicializaMatriz_rhotil(int IMAX, int JMAX, int nDim){

  int i, j, k;

 rhotil = new double**[IMAX];

  for (int i = 0; i<IMAX; i++) {
      rhotil[i] = new double*[JMAX];
      for (int j = 0; j<JMAX; j++) {
          rhotil[i][j] = new double[nDim];
          for (int k = 0; k<nDim; k++) {
                rhotil[i][j][k] = 0.0;
          }
       }
  }

}


void inicializaMatriz_LPhi(int IMAX, int JMAX){

  int i, j;

  LPhi = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        LPhi [i] = new double[JMAX];
    }

    // Inicializando a matriz com as condições de contorno de : Entrada, superior e saída

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 LPhi [i][j] = 0.0;
         }
    }
}

void inicializaMatriz_Cij(int IMAX, int JMAX){ // UPDATE Zerando e ja assumindo chute inicial

  int i, j;

  Cij = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        Cij [i] = new double[JMAX];
    }

    // Inicializando a matriz com as condições de contorno de : Entrada, superior e saída

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                Cij [i][j] = 0.0;
             }
    }
}


void inicializaMatriz_Rx( int JMAX){
   int j;

    Rx = new double [JMAX];

    for (j = 0; j < JMAX; j++){

           Rx[j] = 0.0;
    }
}

void inicializaMatriz_Ry( int JMAX){
   int j;

    Ry = new double [JMAX];

    for (j = 0; j < JMAX; j++){

           Ry[j] = 0.0;
    }
}


void inicializaMatriz_x(int IMAX, int JMAX){ // MATRIX DE RESIDUOS

  int i, j;

  x = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        x[i] = new double[JMAX];
    }

    // Inicializando a matriz com as condições de contorno de : Entrada, superior e saída

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 x[i][j] = 0.0;
         }
    }
}

void inicializaMatriz_y(int IMAX, int JMAX){ // MATRIX DE RESIDUOS

  int i, j;

  y = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        y[i] = new double[JMAX];
    }

    // Inicializando a matriz com as condições de contorno de : Entrada, superior e saída

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 y[i][j] = 0.0;
         }
    }
}


// Inicializando e alocando espaço para as variaveis calculadas no algoritmo de Thomas
void inicializaMatriz_AT( int IMAX){
   int j;

    AT = new double [IMAX];

    for (j = 0; j < IMAX; j++){

           AT[j] = 0.0;
    }
}

void inicializaMatriz_BT( int IMAX){
   int j;

    BT = new double [IMAX];

    for (j = 0; j < IMAX; j++){

           BT[j] = 0.0;
    }
}

void inicializaMatriz_CT( int IMAX){
   int j;

    CT = new double [IMAX];

    for (j = 0; j < IMAX; j++){

           CT[j] = 0.0;
    }
}

void inicializaMatriz_DT( int IMAX){
   int j;

    DT = new double [IMAX];

    for (j = 0; j <IMAX; j++){

           DT[j] = 0.0;
    }
}


void inicializaMatriz_ATX( int IMAX){
   int j;

    ATX = new double [IMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (j = 0; j < IMAX; j++){

           ATX[j] = 0.0;
    }
}

void inicializaMatriz_BTX( int IMAX){
   int j;

    BTX = new double [IMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (j = 0; j < IMAX; j++){

           BTX[j] = 0.0;
    }
}

void inicializaMatriz_CTX( int IMAX){
   int j;

    CTX = new double [IMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (j = 0; j < IMAX; j++){

           CTX[j] = 0.0;
    }
}

void inicializaMatriz_VT( int IMAX){
   int j;

    VT = new double [IMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (j = 0; j <IMAX; j++){

           VT[j] = 0.0;
    }
}


void inicializaMatriz_ATY( int IMAX){
   int j;

    ATY = new double [IMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (j = 0; j < IMAX; j++){

           ATY[j] = 0.0;
    }
}

void inicializaMatriz_BTY( int IMAX){
   int j;

    BTY = new double [IMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (j = 0; j < IMAX; j++){

           BTY[j] = 0.0;
    }
}

void inicializaMatriz_CTY( int IMAX){
   int j;

    CTY = new double [IMAX]; // VER AO CERTO ESSA QUANTIDADE DE POSIÇÕES

    for (j = 0; j < IMAX; j++){

           CTY[j] = 0.0;
    }
}

// Matriz dos Coeficientes da Equação 1
void inicializaMatriz_A(int IMAX, int JMAX){

  int i, j;

  A = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        A[i] = new double[JMAX];
    }

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 A[i][j] = 0.0;
         }
    }
}

void inicializaMatriz_B(int IMAX, int JMAX){

  int i, j;

  B = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        B[i] = new double[JMAX];
    }

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 B[i][j] = 0.0;
         }
    }
}

void inicializaMatriz_C(int IMAX, int JMAX){

  int i, j;

  C = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        C[i] = new double[JMAX];
    }

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 C[i][j] = 0.0;
         }
    }
}

void inicializaMatriz_D(int IMAX, int JMAX){

  int i, j;

  D = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        D[i] = new double[JMAX];
    }

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 D[i][j] = 0.0;
         }
    }
}


void inicializaMatriz_AMX(int IMAX, int JMAX){ // UPDATE Zerando e ja assumindo chute inicial

  int i, j;

  AMX = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        AMX[i] = new double[IMAX];
    }

    // Inicializando a matriz com as condições de contorno de : Entrada, superior e saída

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < IMAX; j++){
                AMX[i][j] = 0.0;
             }
    }
}

void inicializaMatriz_AMY(int IMAX, int JMAX){ // UPDATE Zerando e ja assumindo chute inicial

  int i, j;

  AMY = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        AMY[i] = new double[IMAX];
    }

    // Inicializando a matriz com as condições de contorno de : Entrada, superior e saída

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < IMAX; j++){
                AMY[i][j] = 0.0;
             }
    }
}


void inicializaMatriz_f(int IMAX, int JMAX){ // MATRIX DE RESIDUOS

  int i, j;

  f = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        f[i] = new double[JMAX];
    }

    // Inicializando a matriz com as condições de contorno de : Entrada, superior e saída

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 f[i][j] = 0.0;
         }
    }
}

void inicializaMatriz_g(int IMAX, int JMAX){ // MATRIX DE RESIDUOS

  int i, j;

  g = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        g[i] = new double[JMAX];
    }

    // Inicializando a matriz com as condições de contorno de : Entrada, superior e saída

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 g[i][j] = 0.0;
         }
    }
}

void inicializaMatriz_DELTAPhi(int IMAX, int JMAX){ // MATRIX DE RESIDUOS

  int i, j;

  DELTAPhi = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        DELTAPhi[i] = new double[JMAX];
    }

    // Inicializando a matriz com as condições de contorno de : Entrada, superior e saída

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                 DELTAPhi[i][j] = 0.0;
         }
    }
}

// INICIALIZA TERMOS DE MÉTRICA A1 , A2 E A3
void inicializaMatriz_A1(int IMAX, int JMAX, int nDim){

  int i, j, k;

 A1 = new double**[IMAX];

  for (int i = 0; i<IMAX; i++) {
      A1[i] = new double*[JMAX];
      for (int j = 0; j<JMAX; j++) {
          A1[i][j] = new double[nDim];
          for (int k = 0; k<nDim; k++) {
               A1[i][j][k] = 0.0;
          }
       }
  }

}

void inicializaMatriz_A2(int IMAX, int JMAX, int nDim){

  int i, j, k;

 A2 = new double**[IMAX];

  for (int i = 0; i<IMAX; i++) {
      A2[i] = new double*[JMAX];
      for (int j = 0; j<JMAX; j++) {
          A2[i][j] = new double[nDim];
          for (int k = 0; k<nDim; k++) {
               A2[i][j][k] = 0.0;
          }
       }
  }

}

void inicializaMatriz_A3(int IMAX, int JMAX, int nDim){

  int i, j, k;

 A3 = new double**[IMAX];

  for (int i = 0; i<IMAX; i++) {
      A3[i] = new double*[JMAX];
      for (int j = 0; j<JMAX; j++) {
          A3[i][j] = new double[nDim];
          for (int k = 0; k<nDim; k++) {
               A3[i][j][k] = 0.0;
          }
       }
  }

}


void inicializaMatriz_J(int IMAX, int JMAX, int nDim){

  int i, j, k;

 J = new double**[IMAX];

  for (int i = 0; i<IMAX; i++) {
      J[i] = new double*[JMAX];
      for (int j = 0; j<JMAX; j++) {
          J[i][j] = new double[nDim];
          for (int k = 0; k<nDim; k++) {
               J[i][j][k] = 0.0;
          }
       }
  }

}

void inicializaMatriz_U(int IMAX, int JMAX, int nDim){

  int i, j, k;

 U = new double**[IMAX];

  for (int i = 0; i<IMAX; i++) {
      U[i] = new double*[JMAX];
      for (int j = 0; j<JMAX; j++) {
          U[i][j] = new double[nDim];
          for (int k = 0; k<nDim; k++) {
               U[i][j][k] = 0.0;
          }
       }
  }

}

void inicializaMatriz_V(int IMAX, int JMAX, int nDim){

  int i, j, k;

 V = new double**[IMAX];

  for (int i = 0; i<IMAX; i++) {
      V[i] = new double*[JMAX];
      for (int j = 0; j<JMAX; j++) {
          V[i][j] = new double[nDim];
          for (int k = 0; k<nDim; k++) {
               V[i][j][k] = 0.0;
          }
       }
  }

}

// LOCAL MACH NUMBER
void inicializaMatriz_MACH(int IMAX, int JMAX){

  int i, j;

  MACH = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        MACH[i] = new double[JMAX];
    }

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                MACH[i][j] = 0.0;
         }
    }
}


// cp escoamento
void inicializaMatriz_cp(int IMAX, int JMAX){

  int i, j;

  cp = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        cp[i] = new double[JMAX];
    }

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                cp[i][j] = 0.0;
         }
    }
}


// CALCULO PARA AS PRORIEDADES NO POTENCIAL INTEIRO
void inicializaMatriz_u(int IMAX, int JMAX){

  int i, j;

  u = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        u[i] = new double[JMAX];
    }

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                u[i][j] = 0.0;
         }
    }
}

void inicializaMatriz_v(int IMAX, int JMAX){

  int i, j;

  v = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        v[i] = new double[JMAX];
    }

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                v[i][j] = 0.0;
         }
    }
}


void inicializaMatriz_etax(int IMAX, int JMAX){

  int i, j;

  etax = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        etax[i] = new double[JMAX];
    }

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                etax[i][j] = 0.0;
         }
    }
}

void inicializaMatriz_etay(int IMAX, int JMAX){

  int i, j;

  etay = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        etay[i] = new double[JMAX];
    }

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                etay[i][j] = 0.0;
         }
    }
}

void inicializaMatriz_qsix(int IMAX, int JMAX){

  int i, j;

  qsix = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        qsix[i] = new double[JMAX];
    }

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                qsix[i][j] = 0.0;
         }
    }
}

void inicializaMatriz_qsiy(int IMAX, int JMAX){

  int i, j;

  qsiy = new double* [IMAX];

    for (i = 0; i < IMAX; i++){
        qsiy[i] = new double[JMAX];
    }

    for (i = 0; i < IMAX; i++){
         for (j = 0; j < JMAX; j++){
                qsiy[i][j] = 0.0;
         }
    }
}

void inicializaMatriz_cplinha( int IMAX){
   int j;

    cplinha = new double [IMAX];

    for (j = 0; j < IMAX; j++){

          cplinha[j] = 0.0;
    }
}


void inicializarMatrizes(int IMAX, int JMAX,int nDim, double Uinf,double gamma, double theta){

    inicializaMatriz_rho(IMAX,JMAX,nDim);

    // Variaveis A, B, C, D, P, Q da Equação 1 Escritas no espaço Computacional
    inicializaMatriz_A(IMAX,JMAX);
    inicializaMatriz_B(IMAX,JMAX);
    inicializaMatriz_C(IMAX,JMAX);
    inicializaMatriz_D(IMAX,JMAX);

    // Inicialização para o Algoritmo de THOMAS
    inicializaMatriz_AT(IMAX);
    inicializaMatriz_BT(IMAX);
    inicializaMatriz_CT(IMAX);
    inicializaMatriz_DT(IMAX);

    // Inicialização para periodica em X
    inicializaMatriz_ATX(IMAX);
    inicializaMatriz_BTX(IMAX);
    inicializaMatriz_CTX(IMAX);
    inicializaMatriz_VT(IMAX);

    // Inicialização para periodica em Y
    inicializaMatriz_ATY(IMAX);
    inicializaMatriz_BTY(IMAX);
    inicializaMatriz_CTY(IMAX);

    inicializaMatriz_x(IMAX,JMAX);
    inicializaMatriz_y(IMAX,JMAX);

    // Inicializando variaveis para Geração da PotFlow02 Local de referencia
    inicializaMatriz_Rx(JMAX);
    inicializaMatriz_Ry(JMAX);

    inicializaMatriz_AMX(IMAX,JMAX);
    inicializaMatriz_AMY(IMAX,JMAX);

    inicializaMatriz_f(IMAX,JMAX);
    inicializaMatriz_g(IMAX,JMAX);

    inicializaMatriz_DELTAPhi(IMAX,JMAX); // PARA CORREÇÃO DO POTENCIAL

    // POTENCIAL EM DOMINIO COMPUTACIONAL
    inicializaMatriz_Phi(IMAX,JMAX,Uinf,gamma,theta);
    inicializaMatriz_Phiqsi(IMAX, JMAX, nDim,Uinf,gamma);
    inicializaMatriz_Phieta(IMAX, JMAX, nDim,Uinf,gamma);

    // Matrizes de 3 Dimensões
    inicializaMatriz_rho(IMAX,JMAX,nDim);
    inicializaMatriz_rhobarra(IMAX,JMAX,nDim);
    inicializaMatriz_rhotil(IMAX,JMAX,nDim);

      // Matrizes de 3 Dimensões
    inicializaMatriz_Bbarra(IMAX,JMAX,nDim);
    inicializaMatriz_Btil(IMAX,JMAX,nDim);

    inicializaMatriz_A1(IMAX,JMAX,nDim);
    inicializaMatriz_A2(IMAX,JMAX,nDim);
    inicializaMatriz_A3(IMAX,JMAX,nDim);
    inicializaMatriz_J(IMAX,JMAX,nDim);
    inicializaMatriz_U(IMAX,JMAX,nDim);
    inicializaMatriz_V(IMAX,JMAX,nDim);

   // INICIALIZAÇÃO DOS TERMOS DO POTENCIAL
    inicializaMatriz_LPhi(IMAX,JMAX);
    inicializaMatriz_Cij(IMAX,JMAX);

    // PROPRIEDADES DO POTENCIAL
    inicializaMatriz_etax(IMAX,JMAX);
    inicializaMatriz_qsix(IMAX,JMAX);
    inicializaMatriz_etay(IMAX,JMAX);
    inicializaMatriz_qsiy(IMAX,JMAX);

    inicializaMatriz_v(IMAX,JMAX);
    inicializaMatriz_u(IMAX,JMAX);

    // LOCAL MACH NUMBER
    inicializaMatriz_MACH(IMAX,JMAX);
    inicializaMatriz_cp(IMAX,JMAX);
    inicializaMatriz_Ucampo(IMAX,JMAX);
    inicializaMatriz_rhocampo(IMAX,JMAX);


    inicializaMatriz_cplinha(IMAX);


 }

};

void TridiagPhi(double* a, double* b, double* c, double* dv, double *x, int n);
void TDMAPHI(double* a, double* b, double* c, double* d, double *x, int n);
double MaxResiduoPhi (PotFlow02 Pflow,int IMAX,int JMAX,int kiter, double resmaxX);
double MaxResiduoY (PotFlow02 Pflow,int IMAX,int JMAX,int kiter, double resmaxY);
void TransMatrixVetor(double** Mx,double** My, int IMAX, int JMAX, double* xVet, double* yVet);
void TransMatrixVetorTriplo(double** Mx,double** My,double** Mz, int IMAX, int JMAX, double* xVet, double* yVet, double* zVet);
double valorAbsoluto(double x);
PotFlow02 JacMetr (PotFlow02 Pflow,int i, int j, int IMAX, int JMAX, double omega, double alfa, double xqsi, double yqsi, double xeta, double yeta, double gamma, double gammadiv, double ainf, double c);
PotFlow02 AF2 (PotFlow02 Pflow,int i, int j, int IMAX, int JMAX, double omega, double alfa, double xqsi, double yqsi, double xeta, double yeta, double gamma, double gammadiv, double ainf);
PotFlow02 AF2alfa (PotFlow02 Pflow,int i, int j, int IMAX, int JMAX, double omega,double xqsi, double yqsi, double xeta, double yeta, double gamma, double gammadiv, double ainf, double alfa, int Melements,double alfaLow,double alfaHigh);
double vimaismeio(PotFlow02 Pflow,int i, int jj,double c);
double vjmaismeio(PotFlow02 Pflow,int i, int jj,double c, int JMAX);
PotFlow02 VelContra(PotFlow02 Pflow,int i, int j, int IMAX, int JMAX, double omega, double alfa, double xqsi, double yqsi, double xeta, double yeta, double gamma, double gammadiv, double ainf, double c);
PotFlow02 Densidades (PotFlow02 Pflow,int i, int j, int IMAX, int JMAX, double omega, double alfa, double xqsi, double yqsi, double xeta, double yeta, double gamma, double gammadiv, double ainf, double c, double rhoinf);
double Rparametro(PotFlow02 Pflow,int i, int jj);
double Sparametro(PotFlow02 Pflow,int i, int jj);
PotFlow02 Densidadetilbarra (PotFlow02 Pflow,int i, int j, int IMAX, int JMAX, double omega, double alfa, double xqsi, double yqsi, double xeta, double yeta, double gamma, double gammadiv, double ainf, double c, double rhoinf);
PotFlow02 Residuo (PotFlow02 Pflow,int i, int j, int IMAX, int JMAX, double omega, double alfa, double xqsi, double yqsi, double xeta, double yeta, double gamma, double gammadiv, double ainf);
PotFlow02 PropriedadesEscoamento (PotFlow02 Pflow,int i, int j, int IMAX, int JMAX, double Uinf, double ainf);
void AlfaSequence (int Melements, double alfaLow, double alfaHigh,double *alfaa);

int main(int argc, char*argv[]){

int i,j,jj,kiter=1; // i direção csi e j direção eta
double theta = 0.0; // Para Implementação de ângulo de ataque

int IMAX = 93; // Número maximo de pontos na direção Csi >> i default 93 e 15
int JMAX = 15; // Número maximo de pontos na direção Eta >> j
int nDim = 3; // Número de dimensões do Vetor/Matrix

clock_t TempoInicial;
TempoInicial = clock();

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         GRID PARAMETERS                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

//double res;
double  resmaxX = 0.0, resmaxY = 0.0;
int maxit = 2000;

double xqsi, yqsi, xeta, yeta;
// Vetores para Transformação da Matrix para imprimir no formato Fortran (permite impressão tecplot)
double* xVet = new double[IMAX*JMAX];
double* yVet = new double[IMAX*JMAX];
double* zVet = new double[IMAX*JMAX];

double *resmaxiterPhi, RMAXPHI = 1.0;
resmaxiterPhi = new double [maxit];

double r = 1.9; // Fator de Relaxação do SLOR

// Iterations parameters AF2
double Machinf = 0.84;
double c = 1.2;
double alfa = 2.0;
double omega = 1.8;



// Convergence Criterion
double eps = 1.0*pow(10,-12);
double ainf = 340.00;
double gamma = 1.40;
double Uinf = pow(((gamma + 1.0)/(gamma - 1.0 + (2.0/(Machinf*Machinf)))),0.5);
double rhoinf = pow(2.0/(Machinf*Machinf*(gamma-1.0)+ 2.0),(1.0/(gamma-1.0)));
double gammadiv = ((gamma - 1.0)/(gamma + 1.0));

PotFlow02 Pflow;
Pflow.inicializarMatrizes(IMAX,JMAX,nDim,Uinf,gamma,theta);

double** x_aux;
double** y_aux;


// Perfil Biconvexo Malha Eliptica
x_aux = carregaMatriz("Projeto2/Graficos/MalhaElipticaX",IMAX,JMAX);
y_aux = carregaMatriz("Projeto2/Graficos/MalhaElipticaY",IMAX,JMAX);


// Perfil Biconvexo Malha PARABOLICA
//x_aux = carregaMatriz("Projeto2/Graficos/MalhaParabolicaBiconvexoX",IMAX,JMAX);
//y_aux = carregaMatriz("Projeto2/Graficos/MalhaParabolicaBiconvexoY",IMAX,JMAX);

// Perfil NACA0012 malha Eliptica
//x_aux = carregaMatriz("Projeto2/Graficos/MalhaElipticaNACAX",IMAX,JMAX);
//y_aux = carregaMatriz("Projeto2/Graficos/MalhaElipticaNACAY",IMAX,JMAX);


// Carrega a MALHA COMPUTACIONAL DO PRONETO 2 PARA AS COORDENADAS X E Y
for (int i=0;i<IMAX;i++){
        int kk = JMAX-1;
   for (int j=0;j<JMAX;j++){
        Pflow.x[i][kk] = x_aux[i][j];
        Pflow.y[i][kk] = -y_aux[i][j];

        //cout << Pflow.x[i][kk] << " " << Pflow.y[i][kk] << endl;
        kk = kk-1;
    }
}

//cout << Pflow.x[0][JMAX-1] << "  " << Pflow.x[1][JMAX-1] << endl;
//cout << Pflow.x[0][JMAX-1] - Pflow.x[1][JMAX-1] << endl;
//cout << 1.0/(Pflow.x[0][JMAX-1] - Pflow.x[1][JMAX-1]) << endl;

//Imprimi Malha Eliptica Biconvexo
//TransMatrixVetor(Pflow.x,Pflow.y,IMAX,JMAX,xVet,yVet);
//char MalhaElipticaBicProj3[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto3/Graficos/MalhaElipticaBicProj3.plt";
//WriteFile::writeVetAgrup(xVet,yVet,(IMAX*JMAX),MalhaElipticaBicProj3);

//Matrix::imprimir(Pflow.x,IMAX,JMAX);
//cout << "\n\n\n\n\n" << endl;
//Matrix::imprimir(Pflow.y,IMAX,JMAX);


// CARREGA O CHUTE INICIAL PARA O POTENCIAL EM TODO O DOMINIO
for (j =0; j<=JMAX-1;j++){
    for (i=0; i<=IMAX-1; i++){
        Pflow.Phi[i][j] = Uinf*(Pflow.x[i][j]*cos(theta) + Pflow.y[i][j]*sin(theta));
        //cout << Pflow.Phi[i][j] << endl;
    }
}

//Matrix::imprimir(Pflow.Phi,IMAX,JMAX);
//char Potencial[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto3/Graficos/Potencial.plt";
//WriteFile::writeMatrix(Pflow.Phi,IMAX,JMAX,Potencial);

//cout << "\n" << endl;
//cout << theta << endl;
//cout << Uinf << "\n " << endl;

// Chamada das Metricas

Pflow = JacMetr (Pflow,i,j,IMAX,JMAX,omega,alfa,xqsi, yqsi,xeta, yeta,gamma, gammadiv, ainf,c);

    int Melements = 5;
    double alfaLow = 1.0;
    double alfaHigh = 460.00; // Utilizando o Menor intervalo conforme recomendando por Jameson

    double *alfaa;
    alfaa = new double [Melements];
    int kkk =0;
    AlfaSequence(Melements,alfaLow,alfaHigh,alfaa);


// LOOP PARA CALCULO E CORREÇÃO DO POTENCIAL
for (kiter = 0; kiter < maxit; kiter ++){

    Pflow = VelContra(Pflow,i,jj,IMAX,JMAX,omega,alfa,xqsi,yqsi,xeta,yeta,gamma,gammadiv,ainf,c);
    Pflow = Densidades (Pflow,i,jj,IMAX,JMAX,omega,alfa,xqsi,yqsi,xeta,yeta,gamma,gammadiv,ainf,c,rhoinf);
    Pflow = Densidadetilbarra(Pflow,i,jj,IMAX,JMAX,omega,alfa,xqsi,yqsi,xeta,yeta,gamma,gammadiv,ainf,c,rhoinf);
    Pflow = Residuo(Pflow,i,j,IMAX,JMAX,omega,alfa,xqsi, yqsi,xeta, yeta,gamma, gammadiv, ainf);
    //Pflow = AF2 (Pflow,i,j,IMAX,JMAX,omega,alfa,xqsi, yqsi,xeta, yeta,gamma, gammadiv, ainf);
    Pflow = AF2alfa (Pflow,i, j, IMAX, JMAX,omega,xqsi, yqsi, xeta, yeta, gamma,gammadiv,ainf,alfaa[kkk],Melements,alfaLow,alfaHigh);

    RMAXPHI = MaxResiduoPhi(Pflow,IMAX,JMAX,kiter,resmaxX);
    resmaxiterPhi[kiter] =log10(RMAXPHI);
    cout << resmaxiterPhi[kiter] << endl;

    if ( kkk>Melements-2){

        kkk =0;

    }else{

       kkk = kkk+1;
    }

}

Pflow = PropriedadesEscoamento(Pflow,i,j,IMAX,JMAX,Uinf,ainf);

///////////////////////////////////////////
/////////// BICONVEXO SUBCRITICO //////////
///////////////////////////////////////////

//char ResiduoSubBiconvexo[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto3/Graficos/Biconvexo/ResiduoSubM085a1c1.txt";
//WriteFile::writeVetor(resmaxiterPhi,maxit,ResiduoSubBiconvexo);

//char ResiduoSubBiconvexo[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto3/Graficos/Biconvexo/ResiduoSuba22.txt";
//WriteFile::writeVetor(resmaxiterPhi,maxit,ResiduoSubBiconvexo);
//
//char CplinhaSubBiconvexo[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto3/Graficos/Biconvexo/cplinhaSub.txt";
//WriteFile::writeVetor(Pflow.cplinha,IMAX,CplinhaSubBiconvexo);

//// CALCULO DO CAMPO DE PRESSÃO
//TransMatrixVetorTriplo(Pflow.x,Pflow.y,Pflow.cp,IMAX,JMAX,xVet,yVet,zVet);
//char CampocpSubBiconvexo[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto3/Graficos/Biconvexo/CampocpBiSub.plt";
//WriteFile::writeVetAgrup1(xVet,yVet,zVet,(IMAX*JMAX),CampocpSubBiconvexo);

//
//// CALCULO DO CAMPO DE DENSIDADES
//TransMatrixVetorTriplo(Pflow.x,Pflow.y,Pflow.rhocampo,IMAX,JMAX,xVet,yVet,zVet);
//char CampoRhoSubBiconvexo[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto3/Graficos/Biconvexo/CampoRhoBiSub.plt";
//WriteFile::writeVetAgrup1(xVet,yVet,zVet,(IMAX*JMAX),CampoRhoSubBiconvexo);
////

//// CALCULO DO CAMPO DE VELOCIDADES
//TransMatrixVetorTriplo(Pflow.x,Pflow.y,Pflow.Ucampo,IMAX,JMAX,xVet,yVet,zVet);
//char CampoVelocidadeSubBiconvexo[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto3/Graficos/Biconvexo/CampoVelocidadeBiSub.plt";
//WriteFile::writeVetAgrup1(xVet,yVet,zVet,(IMAX*JMAX),CampoVelocidadeSubBiconvexo);


//// CALCULO DO CAMPO DE MACH
//TransMatrixVetorTriplo(Pflow.x,Pflow.y,Pflow.MACH,IMAX,JMAX,xVet,yVet,zVet);
//char CampoMachSubBiconvexo[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto3/Graficos/Biconvexo/CampoMachBiSubParabolica.plt";
//WriteFile::writeVetAgrup1(xVet,yVet,zVet,(IMAX*JMAX),CampoMachSubBiconvexo);

////////////////////////////////////////////
/////////// BICONVEXO SUPERCRITICO /////////
////////////////////////////////////////////

char ResiduoBiconvexo[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto3/Graficos/Biconvexo/ResiduoSupalfa460.txt";
WriteFile::writeVetor(resmaxiterPhi,maxit,ResiduoBiconvexo);
//
//char CplinhaBiconvexo[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto3/Graficos/Biconvexo/cplinhaSup.txt";
//WriteFile::writeVetor(Pflow.cplinha,IMAX,CplinhaBiconvexo);
//
// //CALCULO DO CAMPO DE PRESSÃO
//TransMatrixVetorTriplo(Pflow.x,Pflow.y,Pflow.cp,IMAX,JMAX,xVet,yVet,zVet);
//char CampocpSuperBiconvexo[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto3/Graficos/Biconvexo/CampocpBiSuper.plt";
//WriteFile::writeVetAgrup1(xVet,yVet,zVet,(IMAX*JMAX),CampocpSuperBiconvexo);
////
// //CALCULO DO CAMPO DE DENSIDADES
//TransMatrixVetorTriplo(Pflow.x,Pflow.y,Pflow.rhocampo,IMAX,JMAX,xVet,yVet,zVet);
//char CampoRhoSuperBiconvexo[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto3/Graficos/Biconvexo/CampoRhoBiSuper.plt";
//WriteFile::writeVetAgrup1(xVet,yVet,zVet,(IMAX*JMAX),CampoRhoSuperBiconvexo);
//
//
//// CALCULO DO CAMPO DE VELOCIDADES
//TransMatrixVetorTriplo(Pflow.x,Pflow.y,Pflow.Ucampo,IMAX,JMAX,xVet,yVet,zVet);
//char CampoVelocidadeSuperBiconvexo[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto3/Graficos/Biconvexo/CampoVelocidadeBiSuper.plt";
//WriteFile::writeVetAgrup1(xVet,yVet,zVet,(IMAX*JMAX),CampoVelocidadeSuperBiconvexo);


//// CALCULO DO CAMPO DE MACH
//TransMatrixVetorTriplo(Pflow.x,Pflow.y,Pflow.MACH,IMAX,JMAX,xVet,yVet,zVet);
//char CampoMachSuperBiconvexo[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto3/Graficos/Biconvexo/CampoMachBiSuper2.plt";
//WriteFile::writeVetAgrup1(xVet,yVet,zVet,(IMAX*JMAX),CampoMachSuperBiconvexo);


///////////////////////////////////////////
/////////// NACA0012 SUBCRITICO //////////
///////////////////////////////////////////


//char ResiduosubNACA[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto3/Graficos/NACA/Residuo.txt";
//WriteFile::writeVetor(resmaxiterPhi,maxit,ResiduosubNACA);
//
//char CplinhasubNACA[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto3/Graficos/NACA/cplinhaNACA.txt";
//WriteFile::writeVetor(Pflow.cplinha,IMAX,CplinhasubNACA);
//
//// CALCULO DO CAMPO DE PRESSÃO
//TransMatrixVetorTriplo(Pflow.x,Pflow.y,Pflow.cp,IMAX,JMAX,xVet,yVet,zVet);
//char CampocpSubNACA[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto3/Graficos/CampocpSubNACA.plt";
//WriteFile::writeVetAgrup1(xVet,yVet,zVet,(IMAX*JMAX),CampocpSubNACA);

 //CALCULO DO CAMPO DE DENSIDADES
//TransMatrixVetorTriplo(Pflow.x,Pflow.y,Pflow.rhocampo,IMAX,JMAX,xVet,yVet,zVet);
//char CampoRhoSubNACAconvexo[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto3/Graficos/Biconvexo/CampoRhoNACASub.plt";
//WriteFile::writeVetAgrup1(xVet,yVet,zVet,(IMAX*JMAX),CampoRhoSubNACAconvexo);


//// CALCULO DO CAMPO DE VELOCIDADES
//TransMatrixVetorTriplo(Pflow.x,Pflow.y,Pflow.Ucampo,IMAX,JMAX,xVet,yVet,zVet);
//char CampoVelocidadeSuperNACA[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto3/Graficos/Biconvexo/CampoVelocidadeNACASuper.plt";
//WriteFile::writeVetAgrup1(xVet,yVet,zVet,(IMAX*JMAX),CampoVelocidadeSuperNACA);

///////////////////////////////////////////
/////////// NACA0012 SUPERCRITICO /////////
///////////////////////////////////////////

//char ResiduoNACA[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto3/Graficos/NACA/Residuo.txt";
//WriteFile::writeVetor(resmaxiterPhi,maxit,ResiduoNACA);
//
//char CplinhaNACA[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto3/Graficos/NACA/cplinhaNACAc18.txt";
//WriteFile::writeVetor(Pflow.cplinha,IMAX,CplinhaNACA);
//
//// CALCULO DO CAMPO DE PRESSÃO
//TransMatrixVetorTriplo(Pflow.x,Pflow.y,Pflow.cp,IMAX,JMAX,xVet,yVet,zVet);
//char CampocpSuperNACA[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto3/Graficos/NACA/CampocpSuperNACA.plt";
//WriteFile::writeVetAgrup1(xVet,yVet,zVet,(IMAX*JMAX),CampocpSuperNACA);

//// CALCULO DO CAMPO DE DENSIDADES
//TransMatrixVetorTriplo(Pflow.x,Pflow.y,Pflow.rhocampo,IMAX,JMAX,xVet,yVet,zVet);
//char CampoRhoSuperNACAconvexo[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto3/Graficos/NACA/CampoRhoNACASuper.plt";
//WriteFile::writeVetAgrup1(xVet,yVet,zVet,(IMAX*JMAX),CampoRhoSuperNACAconvexo);
//
//// CALCULO DO CAMPO DE VELOCIDADES
//TransMatrixVetorTriplo(Pflow.x,Pflow.y,Pflow.Ucampo,IMAX,JMAX,xVet,yVet,zVet);
//char CampoVelocidadeSuperNACA[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto3/Graficos/NACA/CampoVelocidadeNACASuper.plt";
//WriteFile::writeVetAgrup1(xVet,yVet,zVet,(IMAX*JMAX),CampoVelocidadeSuperNACA);

//// CALCULO DO CAMPO DE MACH
//TransMatrixVetorTriplo(Pflow.x,Pflow.y,Pflow.MACH,IMAX,JMAX,xVet,yVet,zVet);
//char CampoMachSuperNACA[] = "C:/Users/Alisson/Desktop/Doutorado ITA-IEAV/CC297 AZEVEDO/Projeto3/Graficos/Biconvexo/CampoMachNACASuper.plt";
//WriteFile::writeVetAgrup1(xVet,yVet,zVet,(IMAX*JMAX),CampoMachSuperNACA);

double tempo_total_s = (clock()- TempoInicial)/(double)CLOCKS_PER_SEC;
cout << "Tempo Total de Execucao (segundos) = " << tempo_total_s << endl;
cout << "Numero Total de Iterações = " <<kiter<<endl;
cout << " Tempo por Iteração" << tempo_total_s /kiter << endl;

return 0;

};


PotFlow02 JacMetr (PotFlow02 Pflow,int i, int j, int IMAX, int JMAX, double omega, double alfa, double xqsi, double yqsi, double xeta, double yeta, double gamma, double gammadiv, double ainf, double c){

    double xqsi_imeio, xqsi_jmeio, yqsi_imeio, yqsi_jmeio, xeta_imeio, xeta_jmeio, yeta_imeio, yeta_jmeio;
    double qsix_imeio, qsix_jmeio,qsiy_imeio, qsiy_jmeio;
    double etax_imeio, etax_jmeio, etay_imeio, etay_jmeio;

    int jj,r,s;

for (jj = 0; jj<=JMAX-2; jj++){ // DA FRONTEIRA EXTERNA ATÉ JMAX-2 (ANTES DA PAREDE)

    for (i = 0; i<=IMAX-2; i++){

    if (i == 0){
    // i,j
    xqsi = 0.5*(Pflow.x[i+1][jj] - Pflow.x[IMAX-2][jj]); // Pegar IMAX-2 pois i = 0 !!!!
    yqsi = 0.5*(Pflow.y[i+1][jj] - Pflow.y[IMAX-2][jj]); // Pegar IMAX-2 pois i = 0 !!!!
    // i,j+1/2
    xqsi_jmeio = 0.25*(Pflow.x[i+1][jj+1] + Pflow.x[i+1][jj] - Pflow.x[IMAX-2][jj+1] - Pflow.x[IMAX-2][jj]); // Pegar IMAX-2 pois i = 0 !!!!
    yqsi_jmeio = 0.25*(Pflow.y[i+1][jj+1] + Pflow.y[i+1][jj] - Pflow.y[IMAX-2][jj+1] - Pflow.y[IMAX-2][jj]); // Pegar IMAX-2 pois i = 0 !!!!

    }else{
    // i,j
    xqsi = 0.5*(Pflow.x[i+1][jj] - Pflow.x[i-1][jj]); // centrado normal
    yqsi = 0.5*(Pflow.y[i+1][jj] - Pflow.y[i-1][jj]); // centrado normal
    // i,j+1/2
    xqsi_jmeio = 0.25*(Pflow.x[i+1][jj+1] + Pflow.x[i+1][jj] - Pflow.x[i-1][jj+1] - Pflow.x[i-1][jj]);
    yqsi_jmeio = 0.25*(Pflow.y[i+1][jj+1] + Pflow.y[i+1][jj] - Pflow.y[i-1][jj+1] - Pflow.y[i-1][jj]);
    }

    if (jj == 0){
    // i,j
    xeta = Pflow.x[i][jj+1] - Pflow.x[i][jj];
    yeta = Pflow.y[i][jj+1] - Pflow.y[i][jj];
    // i+1/2,j
    xeta_imeio = 0.5*( (Pflow.x[i][jj+1] - Pflow.x[i][jj])+ (Pflow.x[i+1][jj+1] - Pflow.x[i+1][jj]) );
    yeta_imeio = 0.5*( (Pflow.y[i][jj+1] - Pflow.y[i][jj])+ (Pflow.y[i+1][jj+1] - Pflow.y[i+1][jj]) );

    }else{
    // i,j
    xeta = 0.5*(Pflow.x[i][jj+1] - Pflow.x[i][jj-1]);
    yeta = 0.5*(Pflow.y[i][jj+1] - Pflow.y[i][jj-1]);
    // i+1/2,j
    xeta_imeio = 0.25*(Pflow.x[i+1][jj+1] + Pflow.x[i][jj+1] - Pflow.x[i+1][jj-1] - Pflow.x[i][jj-1]);
    yeta_imeio = 0.25*(Pflow.y[i+1][jj+1] + Pflow.y[i][jj+1] - Pflow.y[i+1][jj-1] - Pflow.y[i][jj-1]);

    }

    //i+1/2,j
    xqsi_imeio = (Pflow.x[i+1][jj] - Pflow.x[i][jj]);
    yqsi_imeio = (Pflow.y[i+1][jj] - Pflow.y[i][jj]);
    // i,j+1/2
    xeta_jmeio = (Pflow.x[i][jj+1] - Pflow.x[i][jj]);
    yeta_jmeio = (Pflow.y[i][jj+1] - Pflow.y[i][jj]);

    Pflow.J[i][jj][0] = 1.0/(xqsi*yeta - xeta*yqsi); // i,j
    Pflow.J[i][jj][1] = 1.0/(xqsi_imeio*yeta_imeio - xeta_imeio*yqsi_imeio); // i+1/2,j
    Pflow.J[i][jj][2] = 1.0/(xqsi_jmeio*yeta_jmeio - xeta_jmeio*yqsi_jmeio); // i,j+1/2

    // Termos da Matrix de Metrica
    Pflow.qsix[i][jj] = Pflow.J[i][jj][0]*yeta;
    qsix_imeio = Pflow.J[i][jj][1]*yeta_imeio;
    qsix_jmeio  = Pflow.J[i][jj][2]*yeta_jmeio;

    Pflow.qsiy[i][jj] = - Pflow.J[i][jj][0]*xeta;
    qsiy_imeio = - Pflow.J[i][jj][1]*xeta_imeio;
    qsiy_jmeio  = - Pflow.J[i][jj][2]*xeta_jmeio;

    Pflow.etax[i][jj] = - Pflow.J[i][jj][0]*yqsi;
    etax_imeio = - Pflow.J[i][jj][1]*yqsi_imeio;
    etax_jmeio = - Pflow.J[i][jj][2]*yqsi_jmeio;

    Pflow.etay[i][jj] = Pflow.J[i][jj][0]*xqsi;
    etay_imeio = Pflow.J[i][jj][1]*xqsi_imeio;
    etay_jmeio = Pflow.J[i][jj][2]*xqsi_jmeio;

    // Coeficientes de Metrica
    Pflow.A1[i][jj][0] = Pflow.qsix[i][jj]*Pflow.qsix[i][jj] + Pflow.qsiy[i][jj]*Pflow.qsiy[i][jj];
    Pflow.A1[i][jj][1] = qsix_imeio*qsix_imeio + qsiy_imeio*qsiy_imeio;
    Pflow.A1[i][jj][2] = qsix_jmeio*qsix_jmeio + qsiy_jmeio*qsiy_jmeio;

    Pflow.A2[i][jj][0] = Pflow.qsix[i][jj]*Pflow.etax[i][jj] + Pflow.qsiy[i][jj]*Pflow.etay[i][jj];
    Pflow.A2[i][jj][1] = qsix_imeio*etax_imeio + qsiy_imeio*etay_imeio;
    Pflow.A2[i][jj][2] = qsix_jmeio*etax_jmeio + qsiy_jmeio*etay_jmeio;

    Pflow.A3[i][jj][0] = Pflow.etax[i][jj]*Pflow.etax[i][jj] + Pflow.etay[i][jj]*Pflow.etay[i][jj];
    Pflow.A3[i][jj][1] = etax_imeio*etax_imeio + etay_imeio*etay_imeio;
    Pflow.A3[i][jj][2] = etax_jmeio*etax_jmeio + etay_jmeio*etay_jmeio;

   }

    Pflow.qsix[IMAX-1][jj] = Pflow.qsix[0][jj];
    Pflow.etax[IMAX-1][jj] = Pflow.etax[0][jj];

    Pflow.qsiy[IMAX-1][jj] = Pflow.qsiy[0][jj];
    Pflow.etay[IMAX-1][jj] = Pflow.etay[0][jj];

    // ATRIBUINDO VALORES PARA O JACOBIANO NA ESTEIRA EM i=0 E i=imax-1
    Pflow.J[IMAX-1][jj][0] = Pflow.J[0][jj][0];  // i,j
    Pflow.J[IMAX-1][jj][1] = Pflow.J[0][jj][1]; // i+1/2,j
    Pflow.J[IMAX-1][jj][2] = Pflow.J[0][jj][2]; // i,j+1/2

    Pflow.A1[IMAX-1][jj][0] = Pflow.A1[0][jj][0];
    Pflow.A1[IMAX-1][jj][1] = Pflow.A1[0][jj][1];
    Pflow.A1[IMAX-1][jj][2] = Pflow.A1[0][jj][2];

    Pflow.A2[IMAX-1][jj][0] = Pflow.A2[0][jj][0];
    Pflow.A2[IMAX-1][jj][1] = Pflow.A2[0][jj][1];
    Pflow.A2[IMAX-1][jj][2] = Pflow.A2[0][jj][2];

    Pflow.A3[IMAX-1][jj][0] = Pflow.A3[0][jj][0];
    Pflow.A3[IMAX-1][jj][1] = Pflow.A3[0][jj][1];
    Pflow.A3[IMAX-1][jj][2] = Pflow.A3[0][jj][2];

}

// CALCULO DE TODOS OS TERMOS EM JMAX-1 (NA PAREDE !!!)
jj = JMAX-1; // Na parede não precisa calcular os pontos em i,j+1/2

for (i = 0; i<=IMAX-2; i++){

    if (i==0){
    xqsi = 0.5*(Pflow.x[i+1][jj] - Pflow.x[IMAX-2][jj]);
    yqsi = 0.5*(Pflow.y[i+1][jj] - Pflow.y[IMAX-2][jj]);

    }else{
    // i,j usando Backward
    xqsi = 0.5*(Pflow.x[i+1][jj] - Pflow.x[i-1][jj]);
    yqsi = 0.5*(Pflow.y[i+1][jj] - Pflow.y[i-1][jj]);
    }

    // i,j usando Backward
    xeta = 0.5*(Pflow.x[i][jj] - Pflow.x[i][jj-1]);
    yeta = 0.5*(Pflow.y[i][jj] - Pflow.y[i][jj-1]);

        // i,j usando one-sided de três pontos
    //xeta = 0.5*(3*Pflow.x[i][jj] - 4*Pflow.x[i][jj-1] + Pflow.x[i][jj-2]); //operador one-sided de três pontos
    //yeta = 0.5*(3*Pflow.y[i][jj] - 4*Pflow.y[i][jj-1] + Pflow.y[i][jj-2]); //operador one-sided de três pontos

    // i+1/2,j
    xqsi_imeio = (Pflow.x[i+1][jj] - Pflow.x[i][jj]);
    yqsi_imeio = (Pflow.y[i+1][jj] - Pflow.y[i][jj]);

    xeta_imeio = Pflow.x[i][jj] - Pflow.x[i][jj-1];
    yeta_imeio = Pflow.y[i][jj] - Pflow.y[i][jj-1];

    // operando one-sided de três pontos
//    xeta_imeio = 0.5*(3*Pflow.x[i][jj] - 4*Pflow.x[i][jj-1] + Pflow.x[i][jj-2]); //operador one-sided de três pontos
//    yeta_imeio = 0.5*(3*Pflow.y[i][jj] - 4*Pflow.y[i][jj-1] + Pflow.y[i][jj-2]); //operador one-sided de três pontos

    // i,j+1/2
    // Não necessita Calcular essa derivada

    Pflow.J[i][jj][0] = 1.0/(xqsi*yeta - xeta*yqsi); // i,j
    Pflow.J[i][jj][1] = 1.0/(xqsi_imeio*yeta_imeio - xeta_imeio*yqsi_imeio); // i+1/2,j
    // Termos da Matrix de Metrica
    Pflow.qsix[i][jj] = Pflow.J[i][jj][0]*yeta;
    qsix_imeio = Pflow.J[i][jj][1]*yeta_imeio;
    //qsix_jmeio  = Pflow.J[i][jj][2]*yeta_jmeio;
    Pflow.qsiy[i][jj] = - Pflow.J[i][jj][0]*xeta;
    qsiy_imeio = - Pflow.J[i][jj][1]*xeta_imeio;
    //qsiy_jmeio  = - Pflow.J[i][jj][2]*xeta_jmeio;

    Pflow.etax[i][jj] = - Pflow.J[i][jj][0]*yqsi;
    etax_imeio = - Pflow.J[i][jj][1]*yqsi_imeio;
    //etax_jmeio = - Pflow.J[i][jj][2]*yqsi_jmeio;

    Pflow.etay[i][jj] = Pflow.J[i][jj][0]*xqsi;
    etay_imeio = Pflow.J[i][jj][1]*xqsi_imeio;
    //etay_jmeio = Pflow.J[i][jj][2]*xqsi_jmeio;

    // Coeficientes de Metrica
    Pflow.A1[i][jj][0] = Pflow.qsix[i][jj]*Pflow.qsix[i][jj] + Pflow.qsiy[i][jj]*Pflow.qsiy[i][jj];
    Pflow.A1[i][jj][1] = qsix_imeio*qsix_imeio + qsiy_imeio*qsiy_imeio;
    //Pflow.A1[i][jj][2] = qsix_jmeio*qsix_jmeio + qsiy_jmeio*qsiy_jmeio;

    Pflow.A2[i][jj][0] = Pflow.qsix[i][jj]*Pflow.etax[i][jj] + Pflow.qsiy[i][jj]*Pflow.etay[i][jj];
    Pflow.A2[i][jj][1] = qsix_imeio*etax_imeio + qsiy_imeio*etay_imeio;
    //Pflow.A2[i][jj][2] = qsix_jmeio*etax_jmeio + qsiy_jmeio*etay_jmeio;

    Pflow.A3[i][jj][0] = Pflow.etax[i][jj]*Pflow.etax[i][jj] + Pflow.etay[i][jj]*Pflow.etay[i][jj];
    Pflow.A3[i][jj][1] = etax_imeio*etax_imeio + etay_imeio*etay_imeio;
    //Pflow.A3[i][jj][2] = etax_jmeio*etax_jmeio + etay_jmeio*etay_jmeio;
}

    Pflow.qsix[IMAX-1][jj] = Pflow.qsix[0][jj];
    Pflow.etax[IMAX-1][jj] = Pflow.etax[0][jj];

    Pflow.qsiy[IMAX-1][jj] = Pflow.qsiy[0][jj];
    Pflow.etay[IMAX-1][jj] = Pflow.etay[0][jj];

    // ATRIBUINDO VALORES PARA O JACOBIANO NA ESTEIRA EM i=0 E i=imax-1 PARA a parede do perfil, ou seja, JMAX-1
    Pflow.J[IMAX-1][JMAX-1][0] = Pflow.J[0][JMAX-1][0];  // i,j
    Pflow.J[IMAX-1][JMAX-1][1] = Pflow.J[0][JMAX-1][1]; // i+1/2,j
    //Pflow.J[IMAX-1][JMAX-1][2] = Pflow.J[0][JMAX-1][2]; // i,j+1/2

    Pflow.A1[IMAX-1][JMAX-1][0] = Pflow.A1[0][JMAX-1][0];
    Pflow.A1[IMAX-1][JMAX-1][1] = Pflow.A1[0][JMAX-1][1];
    //Pflow.A1[IMAX-1][JMAX-1][2] = Pflow.A1[0][JMAX-1][2];

    Pflow.A2[IMAX-1][JMAX-1][0] = Pflow.A2[0][JMAX-1][0];
    Pflow.A2[IMAX-1][JMAX-1][1] = Pflow.A2[0][JMAX-1][1];
    //Pflow.A2[IMAX-1][JMAX-1][2] = Pflow.A2[0][JMAX-1][2];

    Pflow.A3[IMAX-1][JMAX-1][0] = Pflow.A3[0][JMAX-1][0];
    Pflow.A3[IMAX-1][JMAX-1][1] = Pflow.A3[0][JMAX-1][1];
    //Pflow.A3[IMAX-1][JMAX-1][2] = Pflow.A3[0][JMAX-1][2];

return Pflow;

}

PotFlow02 VelContra(PotFlow02 Pflow,int i, int jj, int IMAX, int JMAX, double omega, double alfa, double xqsi, double yqsi, double xeta, double yeta, double gamma, double gammadiv, double ainf, double c){

    double phiqsi, phiqsi_imeio, phiqsi_jmeio, phieta, phieta_imeio, phieta_jmeio;

for (jj = 0; jj<=JMAX-2; jj++){


    for (i = 0; i<=IMAX-2; i++){

    if (i==0){
    Pflow.Phiqsi[i][jj][0] = 0.5*(Pflow.Phi[i+1][jj] - Pflow.Phi[IMAX-2][jj]);
    Pflow.Phiqsi[i][jj][1] = Pflow.Phi[i+1][jj] - Pflow.Phi[i][jj];
    Pflow.Phiqsi[i][jj][2] = 0.25*(Pflow.Phi[i+1][jj+1] + Pflow.Phi[i+1][jj] - Pflow.Phi[IMAX-2][jj+1] - Pflow.Phi[IMAX-2][jj]);
    }else{
    Pflow.Phiqsi[i][jj][0] = 0.5*(Pflow.Phi[i+1][jj] - Pflow.Phi[i-1][jj]);
    Pflow.Phiqsi[i][jj][1] = Pflow.Phi[i+1][jj] - Pflow.Phi[i][jj];
    Pflow.Phiqsi[i][jj][2] = 0.25*(Pflow.Phi[i+1][jj+1] + Pflow.Phi[i+1][jj] - Pflow.Phi[i-1][jj+1] - Pflow.Phi[i-1][jj]);
    }

    if (jj==0){

    Pflow.Phieta[i][jj][0] = Pflow.Phi[i][0];
    Pflow.Phieta[i][jj][1] = Pflow.Phi[i][0];
    Pflow.Phieta[i][jj][2] = Pflow.Phi[i][jj+1] - Pflow.Phi[i][0];
    //Pflow.Phieta[i][jj][0] = Pflow.Phi[i][jj+1] - Pflow.Phi[i][jj];
    //Pflow.Phieta[i][jj][1] = 0.5*((Pflow.Phi[i][jj+1] - Pflow.Phi[i][jj]) +(Pflow.Phi[i+1][jj+1] - Pflow.Phi[i+1][jj]));
    // Apenas para as derivadas Phiqsi em j=0

    }else{

    Pflow.Phieta[i][jj][0] = 0.5*(Pflow.Phi[i][jj+1] - Pflow.Phi[i][jj-1]);
    //cout << Pflow.Phieta[i][jj][0] << endl;
    Pflow.Phieta[i][jj][1] = 0.25*(Pflow.Phi[i+1][jj+1] + Pflow.Phi[i][jj+1] - Pflow.Phi[i+1][jj-1] - Pflow.Phi[i][jj-1]);
    Pflow.Phieta[i][jj][2] = Pflow.Phi[i][jj+1] - Pflow.Phi[i][jj];

    }

    // CALCULO NO PONTO i,j
    Pflow.U[i][jj][0] = Pflow.A1[i][jj][0]*Pflow.Phiqsi[i][jj][0] + Pflow.A2[i][jj][0]*Pflow.Phieta[i][jj][0];
    Pflow.U[i][jj][1] = Pflow.A1[i][jj][1]*Pflow.Phiqsi[i][jj][1] + Pflow.A2[i][jj][1]*Pflow.Phieta[i][jj][1];
    Pflow.U[i][jj][2] = Pflow.A1[i][jj][2]*Pflow.Phiqsi[i][jj][2] + Pflow.A2[i][jj][2]*Pflow.Phieta[i][jj][2];

    Pflow.V[i][jj][0] = Pflow.A2[i][jj][0]*Pflow.Phiqsi[i][jj][0] + Pflow.A3[i][jj][0]*Pflow.Phieta[i][jj][0];
    //cout << Pflow.V[i][jj][0] << endl;
    Pflow.V[i][jj][1] = Pflow.A2[i][jj][1]*Pflow.Phiqsi[i][jj][1] + Pflow.A3[i][jj][1]*Pflow.Phieta[i][jj][1];
    Pflow.V[i][jj][2] = Pflow.A2[i][jj][2]*Pflow.Phiqsi[i][jj][2] + Pflow.A3[i][jj][2]*Pflow.Phieta[i][jj][2];
}

    Pflow.Phiqsi[IMAX-1][jj][0] = Pflow.Phiqsi[0][jj][0];
    Pflow.Phieta[IMAX-1][jj][0] = Pflow.Phieta[0][jj][0];

    // CALCULO NO PONTO i,j
    Pflow.U[IMAX-1][jj][0] = Pflow.A1[0][jj][0]*Pflow.Phiqsi[0][jj][0] + Pflow.A2[0][jj][0]*Pflow.Phieta[0][jj][0];
    Pflow.V[IMAX-1][jj][0] = Pflow.A2[0][jj][0]*Pflow.Phiqsi[0][jj][0] + Pflow.A3[0][jj][0]*Pflow.Phieta[0][jj][0];
    // CALCULO NO PONTO i+1/2,j
    Pflow.U[IMAX-1][jj][1] = Pflow.A1[0][jj][1]*Pflow.Phiqsi[0][jj][1] + Pflow.A2[0][jj][1]*Pflow.Phieta[0][jj][1];
    Pflow.V[IMAX-1][jj][1] = Pflow.A2[0][jj][1]*Pflow.Phiqsi[0][jj][1] + Pflow.A3[0][jj][1]*Pflow.Phieta[0][jj][1];
    // CALCULO NO PONTO i,j+1/2
    Pflow.U[IMAX-1][jj][2] = Pflow.A1[0][jj][2]*Pflow.Phiqsi[0][jj][2] + Pflow.A2[0][jj][2]*Pflow.Phieta[0][jj][2];
    Pflow.V[IMAX-1][jj][2] = Pflow.A2[0][jj][2]*Pflow.Phiqsi[0][jj][2] + Pflow.A3[0][jj][2]*Pflow.Phieta[0][jj][2];

}

jj = JMAX-1; // PARA PAREDE

for (i = 0; i<=IMAX-2; i++){

    if (i==0){

    Pflow.Phiqsi[i][jj][0] = 0.5*(Pflow.Phi[i+1][jj] - Pflow.Phi[IMAX-2][jj]);

    }else{

    Pflow.Phiqsi[i][jj][0] = 0.5*(Pflow.Phi[i+1][jj] - Pflow.Phi[i-1][jj]);

   }

    // Acho que não precisa calcular
    Pflow.Phiqsi[i][jj][1] = Pflow.Phi[i+1][jj] - Pflow.Phi[i][jj];
    Pflow.Phiqsi[i][jj][2] = Pflow.Phi[i][jj] - Pflow.Phi[i][jj-1];  // i,j+1/2 Utilizar um operador one-sided ou backward na parede

    Pflow.Phieta[i][jj][0] = -(Pflow.A2[i][jj][0]/Pflow.A3[i][jj][0])*Pflow.Phiqsi[i][jj][0]; // Phieta em i,JMAX  (Vi,jmax) = ( A2phiqsi + A3phieta) = 0;
    // Acho que não precisa calcular
    Pflow.Phieta[i][jj][1] = -(Pflow.A2[i][jj][1]/Pflow.A3[i][jj][1])*Pflow.Phiqsi[i][jj][1]; // Phieta em i,JMAX  (Vi,jmax) = ( A2phiqsi + A3phieta) = 0;
    Pflow.Phieta[i][jj][2] = Pflow.Phi[i][jj] - Pflow.Phi[i][jj-1];

    Pflow.V[i][jj][0] = 0.0; // i,j
    Pflow.V[i][jj][1] = 0.0; // i+1/2,j
    Pflow.V[i][jj][2] = - Pflow.V[i][jj-1][1]; // V[i][JMAX-2][1] = V[i][JMAX-1][-1/2] //i,j+1/2

    //cout << Pflow.V[i][jj][0] << " " << Pflow.V[i][jj][1] << " " << Pflow.V[i][jj][2] << endl;

    Pflow.U[i][jj][0] = (Pflow.A1[i][jj][0] - (Pflow.A2[i][jj][0]*Pflow.A2[i][jj][0]/Pflow.A3[i][jj][0]))*Pflow.Phiqsi[i][jj][0];
    Pflow.U[i][jj][1] = (Pflow.A1[i][jj][1] - (Pflow.A2[i][jj][1]*Pflow.A2[i][jj][1]/Pflow.A3[i][jj][1]))*Pflow.Phiqsi[i][jj][1];
    Pflow.U[i][jj][2] = Pflow.U[i][jj-1][2];
}

    Pflow.Phiqsi[IMAX-1][jj][0] = Pflow.Phiqsi[0][jj][0];
    Pflow.Phieta[IMAX-1][jj][0] = Pflow.Phieta[0][jj][0];

   // ATUALIZAÇÃO DA ESTEIRA EM i = 0 E I = IMAX-1
    Pflow.V[IMAX-1][jj][0] = Pflow.V[0][jj][0]; // i,j
    Pflow.V[IMAX-1][jj][1] = Pflow.V[0][jj][1]; // i+1/2,j
    Pflow.V[IMAX-1][jj][2] = - Pflow.V[0][jj-1][1]; // V[i][JMAX-2][1] = V[i][JMAX-1][-1/2] //i,j+1/2

    Pflow.U[IMAX-1][jj][0] = (Pflow.A1[0][jj][0] - (Pflow.A2[0][jj][0]*Pflow.A2[0][jj][0]/Pflow.A3[0][jj][0]));
    Pflow.U[IMAX-1][jj][1] = (Pflow.A1[0][jj][1] - (Pflow.A2[0][jj][1]*Pflow.A2[0][jj][1]/Pflow.A3[0][jj][1]));
    Pflow.U[IMAX-1][jj][2] = - Pflow.U[0][jj-1][1];

return Pflow;

}

PotFlow02 Densidades (PotFlow02 Pflow,int i, int j, int IMAX, int JMAX, double omega, double alfa, double xqsi, double yqsi, double xeta, double yeta, double gamma, double gammadiv, double ainf, double c, double rhoinf){

    int jj,r,s;
    //gammadiv = 0.166667 ; (1.0/(gamma-1.0)) = 2.5
// CALCULO PARA AS DENSIDADES NO CAMPO INTERNO DO ESCOAMENTO
for (jj = 0; jj<=JMAX-1; jj++){

    for (i = 0; i<=IMAX-2; i++){

    if (jj==0){

    Pflow.rho[i][jj][0] = rhoinf;
    Pflow.rho[i][jj][1] = rhoinf;
    }else{

    Pflow.rho[i][jj][0] = pow((1.0 - 0.166667*(Pflow.U[i][jj][0]*Pflow.Phiqsi[i][jj][0] + Pflow.V[i][jj][0]*Pflow.Phieta[i][jj][0])),2.5); // i,j
    Pflow.rho[i][jj][1] = pow((1.0 - 0.166667*(Pflow.U[i][jj][1]*Pflow.Phiqsi[i][jj][1] + Pflow.V[i][jj][1]*Pflow.Phieta[i][jj][1])),2.5); // i+1/2,j
    }

    if (jj ==JMAX-1){
    Pflow.rho[i][jj][2] = Pflow.rho[i][jj-1][2];
    }else{
    Pflow.rho[i][jj][2] = pow((1.0 - 0.166667*(Pflow.U[i][jj][2]*Pflow.Phiqsi[i][jj][2] + Pflow.V[i][jj][2]*Pflow.Phieta[i][jj][2])),2.5); // i,j+1/2
    }

    }
    Pflow.rho[IMAX-1][jj][0] = Pflow.rho[0][jj][0] ; // i,j
    Pflow.rho[IMAX-1][jj][1] = Pflow.rho[0][jj][1] ; // i+1/2,j
    Pflow.rho[IMAX-1][jj][2] = Pflow.rho[0][jj][2]; // i,j+1/2

}
return Pflow;
}

PotFlow02 Densidadetilbarra (PotFlow02 Pflow,int i, int j, int IMAX, int JMAX, double omega, double alfa, double xqsi, double yqsi, double xeta, double yeta, double gamma, double gammadiv, double ainf, double c, double rhoinf){

 int jj,r,s;

// CALCULO DAS DENSIDADES BARRAS E TIL
for (jj = 0; jj<=JMAX-1; jj++){

    for(i=0;i<=IMAX-2;i++){

    int r = Rparametro(Pflow,i,jj); // Direção preferencialmente Upwind
    int s = Sparametro(Pflow,i,jj); // Direção preferencialmente Upwind

    double vi = vimaismeio(Pflow,i,jj,c);
    double vj = vjmaismeio(Pflow,i,jj,c,JMAX);

    if (i==0 && r==-1){
    Pflow.rhotil[i][jj][1] = (1.0 - vimaismeio(Pflow,i,jj,c))*Pflow.rho[i][jj][1] + vimaismeio(Pflow,i,jj,c)*Pflow.rho[IMAX-2][jj][1]; // if i==0 e r==-1 tenho que pegar a posição IMAX-1
    }else {
    Pflow.rhotil[i][jj][1] = (1.0 - vimaismeio(Pflow,i,jj,c) )*Pflow.rho[i][jj][1] + vimaismeio(Pflow,i,jj,c)*Pflow.rho[i+r][jj][1];
    }

    if (jj==0 && s==-1){
    Pflow.rhobarra[i][jj][2] = (1.0 - vjmaismeio(Pflow,i,jj,c,JMAX))*Pflow.rho[i][jj][2] + vjmaismeio(Pflow,i,jj,c,JMAX)*rhoinf;
    }else{
        if (jj==JMAX-1 && s==1){
        Pflow.rhobarra[i][jj][2] = (1.0 - vjmaismeio(Pflow,i,jj,c,JMAX))*Pflow.rho[i][jj][2] + vjmaismeio(Pflow,i,jj,c,JMAX)*Pflow.rho[i][jj+s-1][2];
        //cout << Pflow.rhobarra[i][jj][2] << endl;
        }else{
        Pflow.rhobarra[i][jj][2] = (1.0 - vjmaismeio(Pflow,i,jj,c,JMAX))*Pflow.rho[i][jj][2] + vjmaismeio(Pflow,i,jj,c,JMAX)*Pflow.rho[i][jj+s][2];
        }
    }
    }
    Pflow.rhotil[IMAX-1][jj][1] = Pflow.rhotil[0][jj][1] ;
    Pflow.rhobarra[IMAX-1][jj][2] = Pflow.rhobarra[0][jj][2];
}

for (jj = 1 ; jj<=JMAX-1; jj++){

    for (i = 0; i<=IMAX-2; i++){

    if (i==0){
    Pflow.Btil[i][jj][1] = (Pflow.rhotil[IMAX-2][jj][1]*Pflow.A1[IMAX-2][jj][1])/Pflow.J[IMAX-2][jj][1] ;// CORRESPONDE AO VALOR DE B_TIL EM (i-1/2,j);
    }else{
    Pflow.Btil[i][jj][1] = (Pflow.rhotil[i-1][jj][1]*Pflow.A1[i-1][jj][1])/Pflow.J[i-1][jj][1] ;// CORRESPONDE AO VALOR DE B_TIL EM (i-1/2,j);
    }

    Pflow.Bbarra[i][jj][2] = (Pflow.rhobarra[i][jj-1][2]*Pflow.A3[i][jj-1][2])/Pflow.J[i][jj-1][2];// CORRESPONDE AO VALOR DE B_BARRA EM (i,j-1/2);
    }
    Pflow.Btil[IMAX-1][jj][1] = Pflow.Btil[0][jj][1];
    Pflow.Bbarra[IMAX-1][jj][2] = Pflow.Bbarra[0][jj][2];
}
return Pflow;
}

PotFlow02 Residuo (PotFlow02 Pflow,int i, int j, int IMAX, int JMAX, double omega, double alfa, double xqsi, double yqsi, double xeta, double yeta, double gamma, double gammadiv, double ainf){

// STEP 1 PODE SER RESOLVIDO POR ELIMINAÇÃO DE GAUSS MAS VOU TENTAR USAR O SOLVER DE TRIDIAGONAL SIMPLES
for (int i = 0; i<=IMAX-2; i++){

    for (int jj = 1; jj <=JMAX-1; jj++){

        if(jj== JMAX-1){ // CALCULO DO RESIDUO NA PAREDE É DIFERENTE
            if(i==0){
            Pflow.LPhi[i][jj] = (Pflow.rhotil[i][jj][1]*Pflow.U[i][jj][1]/Pflow.J[i][jj][1]) -(Pflow.rhotil[IMAX-2][jj][1]*Pflow.U[IMAX-2][jj][1]/Pflow.J[IMAX-2][jj][1])- 2.0*(Pflow.rhobarra[i][jj-1][2]*Pflow.V[i][jj-1][2]/Pflow.J[i][jj-1][2]);
            }else{
            Pflow.LPhi[i][jj] = (Pflow.rhotil[i][jj][1]*Pflow.U[i][jj][1]/Pflow.J[i][jj][1]) -(Pflow.rhotil[i-1][jj][1]*Pflow.U[i-1][jj][1]/Pflow.J[i-1][jj][1])- 2.0*(Pflow.rhobarra[i][jj-1][2]*Pflow.V[i][jj-1][2]/Pflow.J[i][jj-1][2]);
            }
        }else{ // CALCULO DO RESIDUO FORA DA PAREDE
            if (i==0){
            Pflow.LPhi[i][jj] = (Pflow.rhotil[i][jj][1]*Pflow.U[i][jj][1]/Pflow.J[i][jj][1]) -(Pflow.rhotil[IMAX-2][jj][1]*Pflow.U[IMAX-2][jj][1]/Pflow.J[IMAX-2][jj][1]) + (Pflow.rhobarra[i][jj][2]*Pflow.V[i][jj][2]/Pflow.J[i][jj][2]) - (Pflow.rhobarra[i][jj-1][2]*Pflow.V[i][jj-1][2]/Pflow.J[i][jj-1][2]);
            }else{
            Pflow.LPhi[i][jj] = (Pflow.rhotil[i][jj][1]*Pflow.U[i][jj][1]/Pflow.J[i][jj][1]) -(Pflow.rhotil[i-1][jj][1]*Pflow.U[i-1][jj][1]/Pflow.J[i-1][jj][1]) + (Pflow.rhobarra[i][jj][2]*Pflow.V[i][jj][2]/Pflow.J[i][jj][2]) - (Pflow.rhobarra[i][jj-1][2]*Pflow.V[i][jj-1][2]/Pflow.J[i][jj-1][2]);
            }
        }
            Pflow.LPhi[IMAX-1][jj] = Pflow.LPhi[0][jj];
    }
}

return Pflow;

}

PotFlow02 AF2 (PotFlow02 Pflow,int i, int j, int IMAX, int JMAX, double omega, double alfa, double xqsi, double yqsi, double xeta, double yeta, double gamma, double gammadiv, double ainf){

    double *x;
    x = new double [IMAX];
    double *y;
    y = new double [IMAX];

    int M = IMAX-2; //
    int n = JMAX-1; //

// STEP 1 PODE SER RESOLVIDO POR backward DE GAUSS MAS VOU TENTAR USAR O
for (int i = 0; i<=IMAX-2; i++){

    for (int jj = 1; jj <=JMAX-1; jj++){
        Pflow.AT [jj-1] = 0.0 ; // LOW DIAGONAL
        Pflow.BT [jj-1] = alfa + Pflow.Bbarra[i][jj][2] ; // MAIN DIAGONAL

       if(jj==JMAX-1){
        Pflow.CT [jj-1] = 0.0;
        }else{
        Pflow.CT [jj-1] = -Pflow.Bbarra[i][jj+1][2] ; // UPPER DIAGONAL
        }

       Pflow.VT[jj-1] = alfa*omega*Pflow.LPhi[i][jj] ; // LADO DIREITO DA EQUAÇÃO
        if ( jj == 1){
         Pflow.AT[0] = 0.0;
        }
        if (jj == n){
        Pflow.CT[n] = 0.0;
       }
    } // END i

    TDMAPHI(Pflow.AT, Pflow.BT, Pflow.CT, Pflow.VT, x, n);

    for ( int jj = 1; jj <=JMAX-1; jj++){
        Pflow.f[i][jj] = x[jj-1];
        Pflow.f[IMAX-1][jj] = Pflow.f[0][jj];
    }

}

// STEP 2 - TRIDIAGONAL PERIODICA
for (int jj = 1; jj<=JMAX-1; jj++){

    for (i = 0; i <=IMAX-2; i++){

        if (i==0){
        Pflow.AT [M] = -Pflow.Btil[i][jj][1]; // ultimo elemento da tridiagonal
        }else{
        Pflow.AT [i-1] = -Pflow.Btil[i][jj][1]; // [1] significa que estou avaliando em i-1/2,j // LOW DIAGONAL
        }

        Pflow.BT [i] = Pflow.Btil[i][jj][1] + Pflow.Btil[i+1][jj][1] + alfa; // MAIN DIAGONAL
        Pflow.CT [i+1] = -Pflow.Btil[i+1][jj][1] ; // UPPER DIAGONAL
        Pflow.VT [i] = Pflow.f[i][jj] + alfa*Pflow.Cij[i][jj-1]; // RHS

    }

        Pflow.CT [0] = -Pflow.Btil[M][jj][1] ;

        TridiagPhi(Pflow.AT, Pflow.BT, Pflow.CT, Pflow.VT, x, M);

        for (i =0; i<=M; i++){
        Pflow.Cij[i][jj] = x[i];
        }
        for ( int i = 0; i <= IMAX-2; i++){
        Pflow.Phi[i][jj] = Pflow.Cij[i][jj] + Pflow.Phi[i][jj];
        }
       // correção para a esteira
        Pflow.Phi[IMAX-1][jj] =  Pflow.Phi[0][jj];

}

return Pflow;

}

PotFlow02 AF2alfa (PotFlow02 Pflow,int i, int j, int IMAX, int JMAX, double omega,double xqsi, double yqsi, double xeta, double yeta, double gamma, double gammadiv, double ainf, double alfa, int Melements,double alfaLow,double alfaHigh){

    double *x;
    x = new double [IMAX];
    double *y;
    y = new double [IMAX];

    int M = IMAX-2; // VOU DEIXAR IMAX-2 POIS IRA RESOLVER ATÉ O PENULTIMO PONTO DA PERIODICA
    int n = JMAX-1; // ANTES ERA JMAX-2

// STEP 1 PODE SER RESOLVIDO POR ELIMINAÇÃO DE GAUSS MAS VOU TENTAR USAR O SOLVER DE TRIDIAGONAL SIMPLES
for (int i = 0; i<=IMAX-2; i++){

    for (int jj = 1; jj <=JMAX-1; jj++){
        Pflow.AT [jj-1] = 0.0 ; // LOW DIAGONAL
        Pflow.BT [jj-1] = alfa + Pflow.Bbarra[i][jj][2] ; // MAIN DIAGONAL
        if(jj==JMAX-1){
        Pflow.CT [jj-1] = 0.0;
        }else{
        Pflow.CT [jj-1] = -Pflow.Bbarra[i][jj+1][2] ; // UPPER DIAGONAL
        }
        Pflow.VT[jj-1] = alfa*omega*Pflow.LPhi[i][jj] ; // LADO DIREITO DA EQUAÇÃO
        if ( jj == 1){
         Pflow.AT[0] = 0.0;
        }
        if (jj == n){
        Pflow.CT[n] = 0.0;
       }
    } // END i

    TDMAPHI(Pflow.AT, Pflow.BT, Pflow.CT, Pflow.VT, x, n);

    for ( int jj = 1; jj <=JMAX-1; jj++){
        Pflow.f[i][jj] = x[jj-1];
        Pflow.f[IMAX-1][jj] = Pflow.f[0][jj];
   }
}

// STEP 2 - TRIDIAGONAL PERIODICA
for (int jj = 1; jj<=JMAX-1; jj++){

    for (i = 0; i <=IMAX-2; i++){

        if (i==0){
        Pflow.AT [M] = -Pflow.Btil[i][jj][1]; // ultimo elemento da tridiagonal
        }else{
        Pflow.AT [i-1] = -Pflow.Btil[i][jj][1]; // [1] significa que estou avaliando em i-1/2,j // LOW DIAGONAL
        }
        Pflow.BT [i] = Pflow.Btil[i][jj][1] + Pflow.Btil[i+1][jj][1] + alfa; // MAIN DIAGONAL
        Pflow.CT [i+1] = -Pflow.Btil[i+1][jj][1] ; // UPPER DIAGONAL
        Pflow.VT [i] = Pflow.f[i][jj] + alfa*Pflow.Cij[i][jj-1]; // RHS
    }
        Pflow.CT [0] = -Pflow.Btil[M][jj][1] ;
        TridiagPhi(Pflow.AT, Pflow.BT, Pflow.CT, Pflow.VT, x, M);
        for (i =0; i<=M; i++){
        Pflow.Cij[i][jj] = x[i];
        }

        for ( int i = 0; i <= IMAX-2; i++){
        Pflow.Phi[i][jj] = Pflow.Cij[i][jj] + Pflow.Phi[i][jj];
        }
        Pflow.Phi[IMAX-1][jj] =  Pflow.Phi[0][jj];
}

return Pflow;

}

PotFlow02 PropriedadesEscoamento (PotFlow02 Pflow,int i, int j, int IMAX, int JMAX, double Uinf, double ainf){

    for (j=1;j<=JMAX-1;j++){

        for (i=0;i<=IMAX-2;i++){
        Pflow.u[i][j] = Pflow.qsix[i][j]*Pflow.Phiqsi[i][j][0] + Pflow.etax[i][j]*Pflow.Phieta[i][j][0];
        Pflow.v[i][j] = Pflow.qsiy[i][j]*Pflow.Phiqsi[i][j][0] + Pflow.etay[i][j]*Pflow.Phieta[i][j][0];
        Pflow.Ucampo[i][j] = sqrt(Pflow.u[i][j]*Pflow.u[i][j] +Pflow.v[i][j]*Pflow.v[i][j] );
        Pflow.rhocampo[i][j] = Pflow.rho[i][j][0];
        Pflow.MACH[i][j] = pow(((2.0*Pflow.Ucampo[i][j]*Pflow.Ucampo[i][j])/(2.4 + Pflow.Ucampo[i][j]*Pflow.Ucampo[i][j]*(0.4))),0.5);
        }
        Pflow.u[IMAX-1][j] = Pflow.u[0][j];
        Pflow.v[IMAX-1][j] = Pflow.v[0][j];
        Pflow.Ucampo[IMAX-1][j] = Pflow.Ucampo[0][j];
        Pflow.rhocampo[IMAX-1][j] = Pflow.rhocampo[0][j];
        Pflow.MACH[IMAX-1][j] = Pflow.MACH[0][j];
    }
    for (j=1;j<=JMAX-1;j++){
        for (i=0;i<=IMAX-2;i++){
        Pflow.cp[i][j] = 1.0 -(Pflow.u[i][j]*Pflow.u[i][j] + Pflow.v[i][j]*Pflow.v[i][j])/(Uinf*Uinf);
        Pflow.cp[IMAX-1][j] = Pflow.cp[0][j];
    }
    }

        for(i=0;i<=IMAX-2;i++){
        Pflow.cplinha[i] = 1.0 - (Pflow.u[i][JMAX-1]*Pflow.u[i][JMAX-1] + Pflow.v[i][JMAX-1]*Pflow.v[i][JMAX-1])/(Uinf*Uinf);
        }
        Pflow.cplinha[IMAX-1] = Pflow.cplinha[0];

return Pflow;

}

void TridiagPhi(double* a, double* b, double* c, double* dv, double *x, int M){

    int i;

    double *l1,*l2,*l3,*u1,*u2,*u3,*f,*g;
    double Ta, Tb;

    l1 = new double [M];
    l2 = new double [M];
    l3 = new double [M];
    u1 = new double [M];
    u2 = new double [M];
    u3 = new double [M];
    f = new double [M];
    g = new double [M];

    for (i =0; i<=M;i++){ // Uma vez que foi alocado dinamicamente (via pointer) é conveniente zerar as posições alocadas, evitando pegar lixo de máquina

        l1[i] = 0.0;
        l2[i] = 0.0;
        l3[i] = 0.0;
        u1[i] = 0.0;
        u2[i] = 0.0;
        u3[i] = 0.0;
        f[i]  = 0.0;
        g[i]  = dv[i];
    }

    l1[0] = b[0]; // L1 = B1 ZERO no c++!!!
    l2[0] = a[0]; // M1 = A1
    l3[0] = c[0]; // N1 = Cn+1
    f[0]  = g[0]/l1[0]; // x_barra1 = L1^-1*y1
    u2[1] = a[M]/l1[0];  // V1 = L1^-1*A0

    Ta = l3[0]*u2[1]; // Ta = N1V1
    Tb = l3[0]*f[0];  //

         /////////////////////////
        // Forward Sweep //////
        ////////////////////////

    for (i = 1; i <=(M-2); i++){

        u1[i] = c[i]/l1[i-1];
        l1[i] = b[i] - u1[i]*l2[i-1];
        l2[i] = a[i];
        l3[i] = -l3[i-1]*u1[i];
        f[i]  = (g[i]-l2[i-1]*f[i-1])/l1[i];
        u2[i+1] = -l2[i-1]*u2[i]/l1[i];
        Ta = Ta + l3[i]*u2[i+1];
        Tb = Tb + l3[i]*f[i];
    }

        u1[M-1] = c[M-1]/l1[M-2];
        l1[M-1] = b[M-1] - l2[M-2]*u1[M-1];
        l2[M-1] = a[M-1] - l3[M-2]*u1[M-1];
         f[M-1] = (g[M-1]-l2[M-2]*f[M-2])/l1[M-1];

         u1[M] = (c[M]-l2[M-2]*u2[M-1])/l1[M-1];
         Ta = Ta + l2[M-1]*u1[M];
         l1[M] = b[M]-Ta;
         f[M] = (g[M] -l2[M-1]*f[M-1]-Tb)/l1[M];

        /////////////////////////
        // Backward Sweep //////
        ////////////////////////

        x[M] = f[M];

        x[M-1] = f[M-1] - u1[M]*x[M];

        for (i = (M-2); i>=0; i--){
            x[i] = f[i] - u1[i+1]*x[i+1] - u2[i+1]*x[M];
        }

 return;

}

void TDMAPHI(double* a, double* b, double* c, double* d, double *x, int n){

    int k;
    double m;

    for ( k = 1 ; k <= n-1; k++){
    m = a[k]/b[k-1];
    b[k] = b[k] -m*c[k-1];
    d[k] = d[k] -m*d[k-1];
    }
    x[n-1] = d[n-1]/b[n-1];
    for (k = n-1; k >= 0; k--){
        x[k] = (d[k] - c[k]*x[k+1])/b[k];
    }

return;

}

double MaxResiduoPhi (PotFlow02 Pflow,int IMAX,int JMAX,int kiter, double resmaxX){

      for (int i = 0; i <= IMAX-2; i++){
            for (int j = 1; j <= JMAX-1; j++ ){
                if (fabs(Pflow.LPhi[i][j])> resmaxX ){
                    resmaxX = fabs(Pflow.LPhi[i][j]);
            }
      }
      }

return resmaxX;

}

double valorAbsoluto(double x){
      if (x >= 0.0) return x;
      else return -x;
    }

double Rparametro(PotFlow02 Pflow,int i, int jj){

      if (Pflow.U[i][jj][1] >= 0.0){
      return -1;
      }else{
      return 1;
     }

}

double Sparametro(PotFlow02 Pflow,int i, int jj){

      if (Pflow.V[i][jj][2] >= 0.0){
     return -1;
      }else{
      return 1;
      }
}

// !!!! WARNING !!! FAZER ESSE CALCULO USANDO COMO REFERÊNCIA A DENSIDADE E NÃO O NUMERO DE MACH !!!!
double vimaismeio(PotFlow02 Pflow,int i, int jj, double c ){

      double c1 = 0.633938, c2 = 4.9325; // c = 1.2 Segundo recomendação do HOLST.
      double A, B, C = 0.0;

      if (Pflow.U[i][jj][1] >= 0.0){
        A = max(0.0,(c1 - Pflow.rho[i][jj][0])*c2*c);
           return A;
      }else{ // SE A VELOCIDADE FOR MENOR QUE ZERO
       B = max(0.0,(c1 - Pflow.rho[i+1][jj][0])*c2*c);
            return B;
       }
}

double vjmaismeio(PotFlow02 Pflow,int i, int jj, double c, int JMAX ){

      double c1 = 0.6339, c2 = 4.9325;
      double A, B, C = 0.0;

     if (Pflow.V[i][jj][2] >= 0.0){

        A = max(0.0,(c1 - Pflow.rho[i][jj][0])*c2*c) ;
            return A;
       }else{
        if (jj ==JMAX-1){
        B = max(0.0,(c1 - Pflow.rho[i][jj][0])*c2*c);
        }else{
        B = max(0.0,(c1 - Pflow.rho[i][jj+1][0])*c2*c);
        }
            return B;
       }
}

void TransMatrixVetor(double** Mx, double** My, int IMAX, int JMAX, double* xVet, double* yVet){

    int t1 = 0, t2 = 0;

    for(int j = 0; j < JMAX; j++){
        for(int i = 0; i < IMAX ; i++){
            xVet[t1++] = Mx[i][j];
            yVet[t2++] = My[i][j];
        }
    }
}

void TransMatrixVetorTriplo(double** Mx, double** My,double** Mz, int IMAX, int JMAX, double* xVet, double* yVet, double* zVet){

    int t1 = 0, t2 = 0, t3 = 0;

    for(int j = 0; j < JMAX; j++){
        for(int i = 0; i < IMAX ; i++){
            xVet[t1++] = Mx[i][j];
            yVet[t2++] = My[i][j];
            zVet[t3++] = Mz[i][j];
            cout << xVet[i] << " " << yVet[i] << " " << zVet[i] << endl;
        }
    }

}

void AlfaSequence (int Melements, double alfaLow, double alfaHigh,double *alfaa){

 double expoente;
 int k, kaux;

    for (k = 0; k<Melements;k++){
        kaux = k+1;
        if (k == Melements){
            kaux = k;
        }
        expoente = (kaux -1.0)/(Melements-1.0);
        alfaa[k] = alfaHigh*pow((alfaLow/alfaHigh),expoente);
    }
return;
}
