
/*Funcoes auxiliares para o teste de Lucas-Lehmer*/


/*
Usada para calcular a profundidade da transformada de Fourier usada no teste de Lucas-Lehmer sobre Mp com base 2^d.
Retorna a potencia de 2 da profundidade, que sera 2^retorno.
*/
.globl profundidade
.type profundidade, @function

/*
Mp = 2^(%rdi) - 1
Base = 2^(%rsi)
*/

profundidade:
  shl %rdi /*2p*/

  mov $0, %rdx
  mov %rdi, %rax  /*%rdi \ rsi%*/
  div %rsi

  cmp $0, %rdx
  je continuar
  add $1, %rax

  continuar:
  bsr %rax, %rcx
  mov $1, %r8
  shl %rcx, %r8
  cmp %r8, %rax
  je fim_prof  /*se %rax nao for uma potencia de 2, adicionamos 1 a %rcx antes de retornar*/

  add $1, %rcx

  fim_prof:
  mov %rcx, %rax
  ret


/*
Recebe um ponteiro para a lista de longs que representa um numero n na forma da FFT, e dois numero p, d.
O numero esta em base 2^d, entao apenas os d primeiros bits de cada long podem estar ocupados (cada um
e um coeficiente do numero na base 2^d).

Faz n (mod Mp), onde Mp = 2^p - 1. n tem no maximo 2p bits, e para fazer o modulo basta tomar os p bits
mas significativos, multiplica-los por 2^(-p) e somar com os p bits menos significativos. Se estourar
algum bit, repetir o processo.

De fato, 2^p == 1 (mod Mp), logo n = n1+2^p + n0 == n1 + n0 (mod Mp).

OBS: Chamamos 2^m a profundidade da FFT (funcao acima) para o numero Mp com base 2^d.
*/
.globl modMp
.type modMp, @function

/*
 OBS: (registro) = o endereco indicado pelo numero no registro
      n(registro) = o endereco n + (registro), onde n esta em bytes 

%rdi: ponteiro para o numero n. Sera n = (%rdi) + 8(%rdi)*2^d + ... + ((2^m - 1)*8)(%rdi)*2^((2^m - 1)*d)
      o numero 2^m e pequeno o suficiente para caber em um long, entao nao ha problema nesses indices
%rsi: p
%rdx: d
*/

modMp:
  mov %rdx, %r8  /*salvar os valores*/

  mov $0, %rdx
  mov %rsi, %rax
  div %r8  /*p/d em %rax e p%d em %rdx. O p-esimo bit de n (primeiro de n1, indice p-1) sera o de indice %rdx em (%rdi, %rax, 8)*/
  
  mov $0, %rcx /*contador para as somas*/
  mov %rax, %r9 /*quociente de p/d guardado*/
  mov $0, %r11

  somas: /*faz todas as somas, exceto a ultima que sera feita separadamente*/
    cmp %r9, %rcx
    je ultima_soma

    push %rcx
  
    mov (%rdi, %rax, 8), %r10
    mov %rdx, %rcx   /*primeira parte do %r10. Ultimos d - %rdx bits do coeficiente atual*/
    shr %rcx, %r10

    mov $1, %r12
    shl %rcx, %r12
    sub $1, %r12  /*%r12 = 1...1, com %rdx 1's.*/
    and %r12, (%rdi, %rax, 8) /*zera os ultimos bits do coeficiente atual, que ja estao em %r10*/

    inc %rax /*atualiza o endereco que sera acessado para montar o %r10 para o fim dessa it e inicio da proxima*/

    mov %r8, %rcx   /*%rcx = d - %rdx. sera usado abaixo*/
    sub %rdx, %rcx
  
    mov (%rdi, %rax, 8), %r13 /*segunda parte do %r10. Primeiros %rdx bits do proximo coeficiente*/
    and %r12, %r13
    shl %rcx, %r13
    add %r13, %r10
  
    mov %r8, %rcx
    mov $1, %r13
    shl %rcx, %r13
    push %r13 /*esse numero sera usado de novo no tratamento de overflow*/
    sub $1, %r13   /*%r13 = 1...1 com d 1's.*/
    sub %r12, %r13 /*%r13 = 1...10...0, com %rdx 0's e d-%rdx 1's.*/
    and %r13, (%rdi, %rax, 8) /*zera os primeiros bits do proximo coeficiente, que ja estao em %r10*/
    
    pop %r13
    pop %rcx

    add %r10, (%rdi, %r11, 8) /*somando %r10 no respectivo coeficiente do numero*/

    push %r11
    call somar_overflow
    pop %r11
    
    inc %rcx /*incrementos para a proxima iteracao*/
    inc %r11

    jmp somas
  

  ultima_soma:
    push %rcx
  
    mov (%rdi, %rax, 8), %r10
    mov %rdx, %rcx   /*primeira parte do %r10. Ultimos d - %rdx bits do coeficiente atual*/
    shr %rcx, %r10

    mov $1, %r12
    shl %rcx, %r12 /*%r12 = 1...1, com %rdx 1's.*/
    sub $1, %r12
    and %r12, (%rdi, %rax, 8) /*zera os ultimos bits do coeficiente atual, que ja estao em %r10*/

    inc %rax /*atualiza o endereco que sera acessado para montar o %r10 para o fim dessa it e inicio da proxima*/

    mov %r8, %rcx   /*%rcx = d - %rdx. sera usado abaixo*/
    sub %rdx, %rcx
  
    mov (%rdi, %rax, 8), %r13 /*segunda parte do %r10. Primeiros %rdx bits do proximo coeficiente*/
    and %r12, %r13
    shl %rcx, %r13
    add %r13, %r10
  
    mov %r8, %rcx
    mov $1, %r13
    shl %rcx, %r13
    push %r13 /*esse numero sera usado de novo no tratamento de overflow*/
    sub $1, %r13   /*%r13 = 1...1 com d 1's.*/
    sub %r12, %r13 /*%r13 = 1...10...0, com %rdx 0's e d-%rdx 1's.*/
    and %r13, (%rdi, %rax, 8) /*zera os primeiros bits do proximo coeficiente, que ja estao em %r10*/
    
    pop %r13
    pop %rcx

    add %r10, (%rdi, %r11, 8) /*somando %r10 no respectivo coeficiente do numero*/

    push %r11
    call somar_overflow
    pop %r11

    ret

/* 
 Verifica se houve overflow na ultima soma feita. Se sim soma 1 ao proximo long int do resultado e chama a si mesma novamente,
 para verificar se houve overflow ao somar esse 1. Se nao, nao faz nada.

 Aqui, o overflow e verificado pelo bit de indice d, ja que cada coeficiente esta em base 2^d. Apenas o ultimo e diferente.
 */
somar_overflow:
  push %rbp
  mov %rsp, %rbp
  
  overflow_loop:
  cmp %r11, %r9
  je ultimo_overflow /*ultimo overflow separado*/

  
  bsr (%rdi, %r11, 8), %r12
  cmp %r12, %r8
  jne sair
  sub %r13, (%rdi, %r11, 8)


  inc %r11
  addq $1, (%rdi, %r11, 8) /*soma 1 ao proximo long int*/
  
  jmp overflow_loop  /*volta para verificar se o incremento acima gerou overflow*/

  sair:
  mov %rbp, %rsp
  pop %rbp
  ret


/*
Se for %r11 = %r9, estamos no coeficiente que contem o p-esimo bit de n. Como estamos fazendo n (mod Mp),
o p-esimo bit deve continuar nulo, logo e nele que a verificacao sera feita. Se for 1, o zeramos, somamos
1 ao primeiro coeficiente e verificamos de novo se houve overflow. Se for 0, nao fazemos nada.
*/
ultimo_overflow:
  bsr (%rdi, %r11, 8), %r12
  jz sair /*necessario para se o coeficiente for zero, pois ano houve overflow e o bsr so vai ativar a zero flag*/
  cmp %r12, %rdx
  jne sair
  
  push %r13
  mov $1, %r13
  push %rcx
  mov %rdx, %rcx
  shl %rcx, %r13  /* %r13 = 10...0 onde 1 esta no bit de indice %rcx*/
  sub %r13, (%rdi, %r11, 8)
  addq $1, (%rdi) 
  
  pop %rcx
  pop %r13

  mov $0, %r11
  call somar_overflow

  jmp sair









/*
Funcao para subtrair 2 de um termo da sequencia do teste de Lucas-Lehmer. Pode ocorrer de seu
primeiro coeficiente ser 1 ou zero, entao teremos que fazer o "vai 1".
*/
.globl sub2
.type sub2, @function

/*
(%rdi) = S[0]
%rsi = d
*/

sub2:
  mov $0, %r8
  cmpq $2, (%rdi)
  jl cont_sub
  
  subq $2, (%rdi)
  ret
  
  cont_sub:
  inc %r8
  cmpq $0, (%rdi, %r8, 8)
  je cont_sub   /*%r8 = indice do primeiro coeficiente nao nulo, exceto possivelmente a_0*/

  subq $1, (%rdi, %r8, 8)
  
  mov %rsi, %rcx
  mov $1, %r9
  shl %rcx, %r9
  dec %r9  /*%r9 = 1...1 com d 1's*/

  dec %r8
  sum_loop: /*soma %r9 a todos os coeficientes de S antes do primeiro nao nulo, exceto ao primeiro coeficiente*/
    cmp $0, %r8
    je fim_sub
    addq %r9, (%rdi, %r8, 8)
    sub $1, %r8
    jmp sum_loop
  
  fim_sub:
  add %r9, (%rdi)
  subq $1, (%rdi)
  ret
