/*
Usa as funções do arquivo biblioteca_Fp.s para fazer contas com polinomios de Fp (p = 2^n).
*/

#include <stdio.h>
#include <stdlib.h>

// print de números em base 2
void print_bin(unsigned long a){
  for(int i = 63; i; i--){
    printf("%li",(a>>i)&1);
  }
  printf("%li\n", a&1);
}


unsigned long div_resto(unsigned long a, unsigned long b);

unsigned long resto(unsigned long a, unsigned long b);

// m deve ser de grau 64, e sera representado sem o primeiro bit para caber em 64 bits
unsigned long mult_mod_gr64(unsigned long a, unsigned long b, unsigned long m);

// m de grau ate 63 (ate 64 bits)
unsigned long mult_mod(unsigned long a, unsigned long b, unsigned long m);

unsigned long inv_mod(unsigned long a, unsigned long m);

// b e um numero natural e os outros sao polinomios de Z_2[x]
unsigned long power_mod(unsigned long a, unsigned long b, unsigned long m);

// b e um numero natural e os outros sao polinomios de Z_2[x]
unsigned long power_mod_gr64(unsigned long a, unsigned long b, unsigned long m);

void main(){
  unsigned long a = 0b0000000000000000000000000000000000000000000000000000000000001101;
  unsigned long b = 0b0000000000000000000000000000000000000000000000000000000000000111;
  unsigned long m = 0b0000000000000000000000000000000000000000000000000000000000011111;

 print_bin(mult_mod(0b110, 0b110, m));
 

  // teste se m e irredutivel: fazer power_mod(p, k, m) com varios ps diferentes, e k = ordem do grupo multiplicativo gerado por m
  // (teorema de Lagrange)
  // print_bin(power_mod(0b1111,  15, m));
}
