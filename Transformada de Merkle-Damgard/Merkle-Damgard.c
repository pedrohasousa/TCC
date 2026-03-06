// Implementacao da trasformada de Merkle-Damgard, sem funcao de compressao

#include <stdio.h>
#include <stdlib.h>

// Usaremos n,m como multiplos de 64 para simplificar. e serao usados
// todos os bits de X, incluindo os zeros a esquerda no ultimo long. Assim, seu comprmiento sera um multiplo de 64


// FUNCAO DE COMPRESSAO ILUSTRATIVA
// Entrada X de comprimento n+m e saída Y de comprimento n. m >=n
// Chave c
void h(unsigned long *X, unsigned long *Y, unsigned long n, unsigned long m, unsigned long c){
  unsigned long r = n/64, t = m/64;

  for(int i = 0; i<r; i++){
    Y[i] = c^X[i]^X[i+t];
  }
}

// Entrada de comprimento qualquer e saída de comprimento n
// n+m e o comprimento da entrada da funcao de compressao usada aqui, m>=n
// c e uma chave para a funcao de compressao
// Escolher l de forma que o texto X tenha comprimento menor que 2^l
// L e o comprimento de X em long ints
// O resultado sera guardado em z, que tem tamanho n
void Merkle_Damgard(unsigned long *X, unsigned long n, unsigned long m, unsigned long c, unsigned long l, unsigned long L, unsigned long *z){

  unsigned long um = 1, r = m>>6, s = L/r, t = L + 2*r - L%r, *Y;

  // Completando o X
  X[0] |= L<<6;
  X[t-L-1] |= um<<63;

  // Geracao do Hash
  t = n>>6;
  Y = malloc((n+m)*sizeof(long));

  for(int i = 0; i<s; i++){
    for(int j = 0; j<t; j++){
      Y[j] = z[j];
    }

    for(int j = 0; j<r; j++){
      Y[t+j] = X[r*i + j];
    }

    h(Y, z, n, m, c);
  }
  free(Y);
}



void main(){
  
  unsigned long r, n = 256, m = 256, k = m>>6, t = n>>6, L = 5, *X, *z0;
  r = L + 2*k - L%k;
  X = malloc(r*sizeof(long)); // X ja com o tamanho necessario para o append
  z0 = malloc(t*sizeof(long));
  
  // Define X usando seus ultimos L longs de tras para frente. Os outros sao apenas para o completamento
  for(int i = 0; i<L; i++){
    X[r-i-1] = 1;
  }

  for(int i = 0; i<t; i++){
    z0[i] = 0;
  }

  Merkle_Damgard(X, n, m, 0, 10, L, z0);
  // h(X, z0, n, m, 0);

  for(int i = 0; i<t; i++){
    printf("%8.8lX ", z0[i]);
  }
  printf("\n");
  
  printf("X: ");
  for(int i = 0; i<r; i++){
    printf("%lX ", X[i]);
  }
  printf("\n");

}