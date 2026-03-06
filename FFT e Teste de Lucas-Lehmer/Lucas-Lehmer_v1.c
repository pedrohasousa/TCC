// Implementacao do teste de primalidade de Lucas-Lehmer
// compilar junto com biblioteca_Zp.s e aux_Lucas_Lehmer.s

#include <stdio.h>
#include <stdlib.h>

// Funcoes da biblioteca_Zp
unsigned long summod(unsigned long a, unsigned long b, unsigned long m);
unsigned long submod(unsigned long a, unsigned long b, unsigned long m);
unsigned long multmod(unsigned long a, unsigned long b, unsigned long m);
unsigned long powermod(unsigned long a, unsigned long b, unsigned long m);
unsigned long invmod(unsigned long a, unsigned long m);

// Profundidade da FFT: Faremos a conta modulo Mp = 2^p -1. Esse numero ocupa p bits (Mp = 0b1...1, p vezes).
// Entao, basta fazer com que o vetor S armazene um numero com o dobro do comprimento de Mp, pois 
// se Sn tem comprimento p e nao foi reduzido, S(n+1) tera comprimento ate 2p, e sera reduzido.

// Logo, precisamos de 2p bits utilizaveis. Em base 2^20, usamos 2p/20 = p/10 long ints. Como p deve ser uma potencia
// de 2, tomamos o menor 2^n >= p/10 arredondao para cima.

// Usada para calcular a profundidade da transformada de Fourier usada no teste de Lucas-Lehmer sobre Mp com base 2^d
// Dennotando-o por m, a rofundidade da FFT sera 2^m, isto e, 2^m e a menor potencia de 2 >= teto(2p/d).
unsigned long profundidade(unsigned long p, unsigned long d);

// Faz S (mod Mp), onde S e um numero vindo da FFT com base 2^d
unsigned long modMp(unsigned long *S, unsigned long p, unsigned long d);

unsigned long sub2(unsigned long *S, unsigned long d);

// O primo usado aqui e 0xffffffff00000001. Podemos ter profundidades de ate 2^32 coeficientes.
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

// Faz S = P*Q
// S deve ter ate 2^n coeficientes. Deve ter 2^(n+1) espacos de memoria alocada, com zeros nos espacos apos seus coeficientes
// R deve ser nulo com 2^(n+1) espacos alocados.
void Mult_FFT(unsigned long *S, unsigned long *R, unsigned long w, unsigned long n, unsigned long p, unsigned long d){
  unsigned long a, b;
  a = powermod(2, n, p);
  b = invmod(a, p);

  FFT(S, R, w, n, p, 1);
  

  for(int i = 0; i < a; i++){
    R[i] = powermod(R[i], 2, p);
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
  for(int i = 0; i < a-1; i++){
    S[i+1] += S[i]>>d;
    S[i] = (S[i]<<(64-d))>>(64-d);
  }
}

// Retorna 1 se 2^p - 1 for primo, e 0 caso contrario
// p e um long,entao p <= 2^64 - 1  
int Lucas_Lehmer(unsigned long p, unsigned long d){
  unsigned long a, m, q, w, *R, *S;
  m = profundidade(p, d);
  q = 0xffffffff00000001;
  w = 1753635133440165772;   //raiz da unidade de ordem 2^32 em Zq
  
  a = powermod(2, m, q); //numero de coeficientes de R e S
  w = powermod(w, powermod(2,32-m,q), q); //agora w tem ordem 2^m
  
  R = malloc(a*sizeof(unsigned long));
  S = malloc(a*sizeof(unsigned long));
  for(int i = 0; i < a; i++){
    R[i] = S[i] = 0;
  }

  S[0] = 4;  //inicializacao da sequencia do teste

  for(int i = 0; i < p-2; i++){
    Mult_FFT(S, R, w, m, q, d);
    modMp(S, p, d);
    sub2(S, d);

    //TESTE
    // if(i == p-3){
    // printf("i = %i ", i);
    // for(int j = 0; j<a; j++){  
    //   printf("%lx ", S[j]);
    // }
    // printf("\n");
    // }

  }

  for(int i = 0; i < a; i++){
    if(S[i] != 0){
      free(S);
      return 0;
    }
  }
  free(S);
  return 1;
}





void main(){
  unsigned long p = 44497;  //44497
  unsigned long d = 20;

  printf("%i\n", Lucas_Lehmer(p, d));
}
