/*
Usa as funções do arquivo biblioteca_Z_p.s
*/

#include <stdio.h>
#include <stdlib.h>

unsigned long multmod(unsigned long a, unsigned long b, unsigned long m);

unsigned long invmod(unsigned long a, unsigned long m);

unsigned long pwm(unsigned long a, unsigned long b, unsigned long m);


void main(){
  unsigned long a = 123; //0xc00500a000900001; 
  unsigned long b = 2; //123123123;
  unsigned long m = 5; //0x80001000d00f0002;

  printf("%lu\n", pwm(a, b, m));
}