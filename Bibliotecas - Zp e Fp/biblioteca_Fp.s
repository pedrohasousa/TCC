fim:
  ret

epilog:
  mov %rbp, %rsp
  pop %rbp
  ret





/*
Divide dois polinômios de Z_2[x], guarda o quociente em %rax e o resto em %rdx
*/
.globl div_resto
.type div_resto, @function

/*
%rdi/%rsi
guarda o quociente em %rax e o resto em %rdx
usa os registros %rdi, %rsi, %rax, %rdx, %r8, %r9 e %rcx

(do jeito que esta, se for para um codigo separado em C, so retorna o quociente. Usar
endereco para retornar o resto tambem)
*/

div_resto:
  mov %rdi, %rdx /*inicializacao do resto e do quociente, respectivamente*/
  mov $0, %rax

  bsr %rdi, %r8 /*contar bits de %rdi e guardar em %r8*/
  jz fim        /*se %rdi = 0, finaliza*/
  bsr %rsi, %r9 /*contar bits de %rsi e guardar em %r9*/
  
  divisao:
  mov %r8, %rcx  /*finaliza se deg(%rsi) > deg(resto atual)*/
  sub %r9, %rcx
  jl fim

  push %rsi
  shl %rcx, %rsi /*incrementando o resto*/
  xor %rsi, %rdx
  
  mov $1, %rsi
  shl %rcx, %rsi /*incrementando o quociente*/
  xor %rsi, %rax
  pop %rsi

  bsr %rdx, %r8 /*contar bits de %rdx e guardar em %r8*/
  jz fim

  jmp divisao











/*
copia de div_resto que guarda o resto em %rax, para retornar apenas o resto para codigo separado em C
*/
.globl resto
.type resto, @function

/*
%rdi/%rsi
guarda o quociente em %rdx e o resto em %rax
usa os registros %rdi, %rsi, %rax, %rdx, %r8, %r9 e %rcx

*/

resto:
  mov %rdi, %rax /*inicializacao do quociente e do resto*/
  mov $0, %rdx

  bsr %rdi, %r8 /*contar bits de %rdi e guardar em %r8*/
  jz fim        /*se %rdi = 0, finaliza*/
  bsr %rsi, %r9 /*contar bits de %rsi e guardar em %r9*/
  
  divisao2:
  mov %r8, %rcx
  sub %r9, %rcx
  jl fim

  push %rsi
  shl %rcx, %rsi /*incrementando o resto*/
  xor %rsi, %rax
  
  mov $1, %rsi
  shl %rcx, %rsi /*incrementando o quociente*/
  xor %rsi, %rdx
  pop %rsi

  bsr %rax, %r8 /*contar bits de %rax e guardar em %r8*/
  jz fim

  jmp divisao2





/*
multiplica dois polinomios de Z_2[x]/m, polinomios de ate 64 bits e guarda o resultado no %rax
m deve ser de grau 64, e dele serao guardados apenas os ultimos 64 bits pois ja sabemos que o primeiro e 1.
*/
.globl mult_mod_gr64
.type mult_mod_gr64, @function

/*
%rdi * %rsi (mod %rdx)
*/

mult_mod_gr64:
  mov $0, %rax  /*retorna 0 caso algum dos polinomios seja nulo*/
  cmp $0, %rdi
  je fim
  cmp $0, %rsi
  je fim

  mov %rdi, %rax
  bsr %rsi, %rcx /*guarda o bit mais significativo de %rsi*/
  mult_loop64:
  cmp $0, %rcx
  je fim

  shl %rax
  jnc continuar_64
  xor %rdx, %rax
  
  continuar_64:
  dec %rcx
  bt %rcx, %rsi
  jnc mult_loop64
  xor %rdi, %rax
  jmp mult_loop64



/*
multiplica dois polinomios de Z_2[x]/m, todos polinomios de ate 64 bits (grau ate 63).
guarda o resultado no %rax
*/
.globl mult_mod
.type mult_mod, @function

/*
%rdi * %rsi (mod %rdx)
*/

mult_mod:
  mov $0, %rax  /*retorna 0 caso algum dos polinomios seja nulo*/
  cmp $0, %rdi
  je fim
  cmp $0, %rsi
  je fim
  
  push %rsi      /*%rdi (mod %rdx)*/
  mov %rdx, %rsi
  push %rdx
  call div_resto
  mov %rdx, %rdi
  pop %rdx
  pop %rsi
  
  bsr %rdx, %r9 /*contar os bits de %rdx e guardar em %r9*/

  mov %rdi, %rax
  bsr %rsi, %rcx /*guarda o bit mais significativo de %rsi*/
  mult_loop:
  bsr %rax, %r8 /*conta os bits de %rax para essa iteracao*/
  cmp $0, %rcx
  je fim

  shl %rax    /*multiplicar %rax por x e reduzir mod %rdx se preciso*/
  inc %r8
  cmp %r9, %r8  /*comparacao usada para ver se precisa reduzir %rax*/
  jl continuar
  xor %rdx, %rax
  
  continuar:
  dec %rcx   /*somar %rdi se o coeficiente da proxima potencia de x em %rsi for 1*/
  bt %rcx, %rsi
  jnc mult_loop
  xor %rdi, %rax
  jmp mult_loop




/*
Calcula o inverso de um polinomio de Z_2[x]/m, polinomios de ate 64 bits (grau ate 63)
Deve ser coprimo com m. (retorna um valor errado se nao for o caso)
*/
.globl inv_mod
.type inv_mod, @function

/*
(%rdi)^(-1) (mod %rsi)

%rax guardara o quociente e %rdx guardara o resto das divisoes feitas no processo
%r11 e %r12 guardarao os coeficientes do algoritmo de euclides extendido referentes ao numero em %rdi
%r13 guardara o valor inicial de %rsi
*/

inv_mod:
  mov $1, %r11
  mov $0, %r12
  mov %rsi, %r13
  
  inv_mod_loop:
  call div_resto
  
  mov %rsi, %rdi
  mov %rdx, %rsi

  push %r12
  
  push %rdi
  push %rsi
  mov %rax, %rdi
  mov %r12, %rsi
  mov %r13, %rdx
  call mult_mod  /* %rax*%r12 (mod %r13) */
  pop %rsi
  pop %rdi

  mov %r11, %r12
  xor %rax, %r12

  pop %r11

  cmp $0, %rsi  /*se o resto for zero, finaliza*/
  jne inv_mod_loop
  mov %r11, %rax
  ret




/*
Faz a^b (mod m), onde a,m sao polinomios de Z_2[x] e b e um numero natural.
Guarda o resultado no rax
*/
.globl power_mod
.type power_mod, @function

/*
%rdi^%rsi (mod %rdx)

registros usados: rdi, rsi, rdx, r10, r11, rcx
*/

power_mod:
  mov %rsi, %r10 /*%r10 vai guardar o valor inicial de %rsi*/
  bsr %r10, %r11 /*%r11 guarda o numero de bits de %rsi*/

  mov $1, %rax /*retorna 1 se %rsi = 0*/
  cmp $0, %rsi
  je fim
  
  push %rdx
  push %rcx
  mov %rdx, %rsi
  call div_resto /*reduzir %rdi (mod %rdx)*/
  mov %rdx, %rdi
  pop %rcx
  pop %rdx

  mov $1, %rax /*recolocando q no rax, por cauda da divisao acima*/
  mov $0, %rcx
  power_loop:
  bt %rcx, %r10 /*verifica se o bit da iteracao atual de %r10 e 1 ou 0*/
  jnc continuar3
  
  push %rdx
  push %rcx
  mov %rax, %rsi
  call mult_mod   /*faz %rdi*%rax se o bit atual de %rsi for 1*/
  pop %rcx
  pop %rdx

  continuar3:
  push %rdx
  push %rax
  push %rcx
  mov %rdi, %rsi
  call mult_mod /*calcula a proxima potencia de %rdi e guarda no proprio %rdi*/
  pop %rcx
  mov %rax, %rdi
  pop %rax
  pop %rdx

  cmp %rcx, %r11 /*verfifica se deve comecar outra iteracao*/
  je fim
  inc %rcx
  jmp power_loop






  /*
Faz a^b (mod m), onde a,m sao polinomios de Z_2[x] e b e um numero natural.
m deve ter grau 64 e a grau ate 63
Guarda o resultado no rax.
*/
.globl power_mod_gr64
.type power_mod_gr64, @function

/*
%rdi^%rsi (mod %rdx)

registros usados: rdi, rsi, rdx, r8
*/

power_mod_gr64:
  mov $1, %rax
  cmp $0, %rsi /*retorna 1 se %rsi = 0*/
  je fim

  mov $0, %rcx
  bsr %rsi, %r8 /*%r8 guarda o numero de bits de %rsi*/
  mov %rsi, %r9 /*%r9 vai guardar o valor inicial de %rsi*/
  power_loop64:
  bt %rcx, %r9 /*verifica se o bit da iteracao atual de %r9 e 1 ou 0*/
  jnc continuar2_64
  
  push %rcx
  mov %rax, %rsi
  call mult_mod_gr64 /*faz %rdi*%rax se o bit atual de %rsi for 1*/
  pop %rcx
  
  continuar2_64:
  push %rcx
  push %rax
  mov %rdi, %rsi
  call mult_mod_gr64 /*calcula a proxima potencia de %rdi e guarda no proprio %rdi*/
  mov %rax, %rdi
  pop %rax
  pop %rcx

  cmp %rcx, %r8 /*verfifica se deve comecar outra iteracao*/
  je fim
  inc %rcx
  jmp power_loop64










/*
Calcula o mdc entre dois polinomios de grau ate 63
*/
.globl mdc
.type mdc, @function
/*
mdc(%rdi, %rsi)
*/
mdc: 
  call resto  /*resto de %rdi/%rsi */ 
  cmp $0, %rax
  je fim_mdc
  mov %rsi, %rdi
  mov %rax, %rsi
  jmp mdc

fim_mdc: 
  mov %rsi, %rax
  ret










/*
Calcula a derivada de um polinomio de grau ate 63
*/
.globl derivada
.type derivada, @function
/*
deriva %rdi
*/

derivada:
  shr %rdi
  mov $0x5555555555555555, %r15
  and %r15, %rdi
  mov %rdi, %rax
  ret
