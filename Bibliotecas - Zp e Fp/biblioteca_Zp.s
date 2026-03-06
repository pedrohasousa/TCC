/*
soma em Z_p
*/
.globl summod
.type summod, @function
/*%rdi + %rsi mod(%rdx)*/

summod:
  mov %rdx, %r8
  mov %rsi, %rax
  add %rdi, %rax
  jc sum1
  mov $0, %rdx
  jmp sum2

  sum1:
    mov $1, %rdx

  sum2:
    div %r8
    mov %rdx, %rax
    ret


/*
subtracao em Z_p
*/
.globl submod
.type submod, @function
/*%rdi - %rsi mod(%rdx)*/

submod:
  sub %rsi, %rdi
  jnc fim_sub
  sub1:
  add %rdx, %rdi
  jnc sub1

  fim_sub:
    mov %rdi, %rax
    ret




/*
multiplicacao em Z_p
*/
.globl multmod
.type multmod, @function

/*
Faz %rdi * %rsi mod(%rdx)
- %r8 guardara o valor inicial de %rdx
*/
multmod:
  mov %rdx, %r8
  mov %rdi, %rax
  mul %rsi
  div %r8
  mov %rdx, %rax
  ret













/*
inverso multiplicativo em Z_p
*/
.globl invmod
.type invmod, @function

/*
(%rdi)^(-1) (mod %rsi)
%rax guardara o quociente e %rdx guardara o resto das divisoes feitas no processo
%r11 e %r12 guardarao os coeficientes do algoritmo de euclides extendido referentes ao numero em %rdi
%r13 guardara o valor inicial de %rsi
*/

invmod:
  mov $1, %r11
  mov $0, %r12
  mov %rsi, %r13
  
  invmod_loop:
  mov $0, %rdx
  mov %rdi, %rax
  div %rsi
  
  mov %rsi, %rdi
  mov %rdx, %rsi


  push %r12
  
  push %rdi
  push %rsi
  mov %rax, %rdi
  mov %r12, %rsi
  mov %r13, %rdx
  call multmod  /* %rax*%r12 (mod %r13) */
  pop %rsi
  pop %rdi

  push %rdi
  push %rsi
  mov %r11, %rdi
  mov %rax, %rsi
  mov %r13, %rdx
  call submod  /*%r11 - %rax (mod %r13)*/
  pop %rsi
  pop %rdi
  mov %rax, %r12

  pop %r11

  cmp $0, %rsi  /*se o resto for zero, finaliza*/
  jne invmod_loop
  mov %r11, %rax
  ret    



 

/*
potenciacao em Z_p
*/
.globl   powermod
.type    powermod, @function

powermod:
  mov %rsi, %r13 /*%r13 vai guardar o valor inicial de %rsi*/
  bsr %r13, %r11 /*%r11 guarda o numero de bits de %rsi*/
  mov %rdx, %r12 /*%r12 guarda o valor inicial de %rdx*/

  mov $1, %rax /*retorna 1 se %rsi = 0*/
  cmp $0, %rsi
  je fim_power

  mov %rdi, %rax
  mov $0, %rdx   /*reduzir %rdi mod(%rdx)*/
  div %r12
  mov %rdx, %rdi

  mov $1, %rax /*recolocando q no rax, por causa da divisao acima*/
  mov $0, %rcx
  power_loop:
  bt %rcx, %r13 /*verifica se o bit da iteracao atual de %r13 e 1 ou 0*/
  jnc continuar3
  
  mov %r12, %rdx
  mov %rax, %rsi
  push %rdi
  call multmod
  pop %rdi

  continuar3:
  mov %r12, %rdx
  push %rax
  mov %rdi, %rsi
  call multmod
  mov %rax, %rdi
  pop %rax

  cmp %rcx, %r11 /*verfifica se deve comecar outra iteracao*/
  je fim_power
  inc %rcx
  jmp power_loop

  fim_power:
    ret
