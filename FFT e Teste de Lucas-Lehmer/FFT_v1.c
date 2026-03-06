// Implementacao do algoritmo FFT para multiplicar inteiros
// compilar com biblioteca_Zp.s

#include <stdio.h>
#include <stdlib.h>

unsigned long summod(unsigned long a, unsigned long b, unsigned long m);
unsigned long submod(unsigned long a, unsigned long b, unsigned long m);
unsigned long multmod(unsigned long a, unsigned long b, unsigned long m);
unsigned long powermod(unsigned long a, unsigned long b, unsigned long m);
unsigned long invmod(unsigned long a, unsigned long m);


//P deve ter 2^n coeficientes
//p e o primo do Zp usado
//Entrar com i = 1
void FFT(unsigned long *P, unsigned long *R, unsigned long w, unsigned long n, unsigned long p, unsigned long i){
  if(n == 0){
    R[0] = P[0];
    return;
  }

  else{
    unsigned long l;
    l = 2*i;
    FFT(P, R, powermod(w,2,p), n-1, p, l);
    FFT(P+i, R+i, powermod(w,2,p), n-1, p, l);
    
    unsigned long k, m, *aux;
    m = powermod(2, n-1, p);
    aux = malloc(2*m*sizeof(unsigned long));

    for(int j = 0; j < m; j++){
      k = powermod(w, j, p);
      k = multmod(k, R[2*j*i+i], p);

      //printf("R[%i] + k[%i] = %lu + %lu\n", j, j, R[j], k);
      // printf("j e j+im  %i %i\n", j, j + i*m);

      aux[j+m] = submod(R[2*j*i], k, p);
      aux[j] = summod(R[2*j*i], k, p);
    }

    for(int j = 0; j < m; j++){
      R[j*i] = aux[j];
      R[j*i+m*i] = aux[j+m];
    }
    free(aux);

    // printf("n = %lu\n", n);
    // for(int j = 0; j < 2*m; j++){
    //   printf("%lx ", R[j*i]);
    // }
    // printf("\n");
  }
}





//Faz S = P*Q
//P e Q devem ter ate 2^n coeficientes. Devem ter 2^(n+1) espacos de memoria alocada, com zeros nos espacos apos seus coeficientes
//R e S devem ser nulos com 2^(n+1) espacos alocados.
void Mult_FFT(unsigned long *P, unsigned long *Q, unsigned long *R, unsigned long *S, unsigned long w, unsigned long n, unsigned long p){
  unsigned long a, b;
  a = powermod(2, n, p);
  b = invmod(a, p);
  printf("a = %lu\nb = %lu\n", a, b);

  FFT(P, R, w, n, p, 1);
  FFT(Q, S, w, n, p, 1);
  

  for(int i = 0; i < a; i++){
    R[i] = multmod(R[i], S[i], p);
  }

  FFT(R, S, w, n, p, 1);
  
  for(int i = 0; i < a; i++){
    S[i] = multmod(S[i], b, p);
  }

  //Inversao do vetor
  unsigned long aux = 0;
  for(int i = 1; i < a/2; i++){
    aux = S[a-i];
    S[a-i] = S[i];
    S[i] = aux;
  }

  //"vai 1" em base 2^d
  unsigned long d = 20;
  d = powermod(2, d, p);
  for(int i = 0; i < a-1; i++){
    S[i+1] += S[i]/d;
    S[i] = S[i]%d;
  }
}


void main(){
  long n = 4;
  unsigned long *P;
  unsigned long *Q;
  unsigned long *R;
  unsigned long *S;
  unsigned long p = 18446744069414584321;  //2^64 - 2^32 + 1 ou 0xFFFFFFFF00000001
  unsigned long w = 1753635133440165772;   //raiz da unidade de ordem 2^32 em Zp
  
  long a;
  a = powermod(2, n, p);

  //--------------------------------Aplicacao da FFT em um polinomio----------------------------------
  // P = malloc(a*sizeof(long)); //polinomios com 2^n coeficientes
  // Q = malloc(a*sizeof(long));
  // R = malloc(a*sizeof(long));
  // S = malloc(a*sizeof(long));

  // ajuste em w, para torna-lo uma raiz da unidade de ordem 2^n

  // w = powermod(w, powermod(2,32-n,p), p);
  //printf("teste %lu\n", w);


  // coeficientes de P e zerando R
  // for(int i = 0; i < a; i++){
  //   P[i] = 0;
  //   R[i] = 0;
  // }
  // P[0] = 1;
  

  // FFT(P, R, w, n, p, 1);

  // for(int i = 0; i < a; i++){
  //   printf("R[%i] = %lu\n", i, R[i]);
  // }



  //-------------------------Multiplicacao de inteiros positivos---------------------------------
  // os numeros tem a = 2^n digitos, logo precisamos de 2*a de espaco para guardar a multiplicacao
  P = malloc(2*a*sizeof(long)); //polinomios com 2^n coeficientes, mas 2^(n+1) espacos alocados
  Q = malloc(2*a*sizeof(long));
  R = malloc(2*a*sizeof(long));
  S = malloc(2*a*sizeof(long));

  w = powermod(w, powermod(2,32-(n+1),p), p);
  // printf("w %lu\n", w);

  for(int i = a; i < 2*a; i++){
    P[i] = Q[i] = R[i] = S[i] = 0;
  }

  for(int i = 0; i < a; i++){
    P[i] = Q[i] = 0;
    R[i] = S[i] = 0;
  }
  
  //Ver na funcao qual dase esta sendo usada
  // P[0] = 0x116;
  // P[1] = 0x1;
  // P[2] = 0x132;
  // P[3] = 0x765;
  P[0] = 0b10001000001000001001;
  P[1] = 0b00001001010001000100;

  Q[0] = 0b10001000001000001001;
  Q[1] = 0b00001001010001000100;
  
  // Q[0] = 0x2;
  // Q[1] = 0xefd;
  // Q[2] = 0x12c;
  // Q[3] = 0x0;

  //entrar com n+1, pois essa e a profundidade da FFT usada para multiplicar inteiros com 2^n digitos
  Mult_FFT(P, Q, R, S, w, n+1, p);

  for(int i = 0; i < 2*a; i++){
    printf("S[%i] : 0%lx;\n", i, S[i]);
  }
  printf("\n");


  // printf("Resultado\n");
  // for(int i = 0; i < a; i++){
  //   printf("R[%i] = %lu ", i, R[i]);
  // }
}
