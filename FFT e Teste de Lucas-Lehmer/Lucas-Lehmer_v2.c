// Implementacao melhorada do teste de primalidade de Lucas-Lehmer
// Eleva números ao quadrado mod Mp usando apenas p bits
// Compilar com biblioteca_Z_p.s e aux_Lucas_Lehmer.s

#include <stdio.h>
#include <stdlib.h>

// Funcoes da biblioteca_Z_p
unsigned long summod(unsigned long a, unsigned long b, unsigned long m);
unsigned long submod(unsigned long a, unsigned long b, unsigned long m);
unsigned long multmod(unsigned long a, unsigned long b, unsigned long m);
unsigned long powermod(unsigned long a, unsigned long b, unsigned long m);
unsigned long invmod(unsigned long a, unsigned long m);

// Funcao de aux_Lucas_Lehmer.s
unsigned long profundidade(unsigned long p, unsigned long d);

unsigned long profund(unsigned long p){
    unsigned long m = 1, d = 32;
    d = 32;  //maior d possivel
    m = profundidade(p, d);

    while(d > 32 - m/2){
        d = 32 - m/2;
        m = profundidade(p, d);
    }

    return m == 0 ? 1 : m;
}


// Por causa da lista de comprimentos, a funcao de subtrair 2 de aux_Lucas_Lehmer.s nao pode ser usada aqui
void sub2novo(unsigned long *S, short int *comprimentos, unsigned long a){
    if(S[0] >= 2){
        S[0] -= 2;
        return;
    }

    long i;
    for(i = 1; i < a; i++){
        if(S[i] > 0){
            break;
        }
    }
    S[i]--;
    i--;
    for(i; i >= 0; i--){
        S[i] += (1<<comprimentos[i]) - 1;
    }
    S[0] -= 1;
}

// O primo usado aqui e 0xffffffff00000001. Podemos ter profundidades de ate 2^32 coeficientes.
void FFT(unsigned long *P, unsigned long* raizes, unsigned long n, unsigned long p){
    unsigned long aux, a = 1, m = powermod(2, n-1, p), j = 0;
    
    for(int l = 0; l < m; l++){
        aux = summod(P[l], P[l + m], p);
        P[l + m] = submod(P[l], P[l + m], p);
        P[l] = aux;
    }

    for (int i = 1; i < n; i++){
        a <<= 1;
        m >>= 1;
        j = 0;
        for(int l = 0; l < a; l++){
            for(int k = 0; k < m; k++){
                aux = summod(P[j + k], multmod(raizes[l], P[j + k + m], p), p);
                P[j + k + m] = submod(P[j + k], multmod(raizes[l], P[j + k + m], p) , p);
                P[j + k] = aux;
            }
            j+=(m<<1);
        }
    }
}



// lista de razies utilizada na FFT de profundidade n, tem apenas 2^(n-1) elementos
// raizes deve ter 2^(n-1) espacos alocados
void criarListaRaizes(unsigned long *raizes, unsigned long n, unsigned long w, unsigned long p){
    if(n == 1){
        raizes[0] = 1;
        return;
    }
    unsigned long m = powermod(2, n-2, p), j = 1<<(n-2);

    raizes[0] = 1;
    for(int i = 1; i <= m; i<<=1){
        raizes[i] = powermod(w, j, p);
        j >>= 1;

        for(int k = i; k < i<<1; k++){
            raizes[k] = multmod(raizes[i], raizes[k-i], p);
        }
    }
}

// inverte as raizes para usar na FFT inversa
void inverterListaRaizes(unsigned long *inversas, unsigned long *raizes, unsigned long n, unsigned long w, unsigned long p){
    unsigned long m = powermod(2, n-1, p);

    for(int i = 0; i < m; i++){
        inversas[i] = invmod(raizes[i], p);
    }
}

void FFTinversa(unsigned long *P, unsigned long* inversas, unsigned long n, unsigned long p){
    unsigned long aux, a = powermod(2, n-1, p), m = 1, j = 0;
    
    for(int l = 0; l < a; l++){
        aux = summod(P[j], P[j + m], p);
        P[j + m] = multmod(inversas[l], submod(P[j], P[j + m], p), p);
        P[j] = aux;
        j+=(m<<1);
    }

    
    for (int i = 1; i < n; i++){
        a >>= 1;
        m <<= 1;
        j = 0;
        for(int l = 0; l < a; l++){
            for(int k = 0; k < m; k++){
                aux = summod(P[j + k], P[j + k + m], p);
                P[j + k + m] = multmod(inversas[l], submod(P[j + k], P[j + k + m], p), p);
                P[j + k] = aux;
            }
            j+=(m<<1);
        }
    }

    // Em vez de multiplicar por 2^(-n) aqui, isso é feito na função voltar_base_sigmap
    return;
}

// p = expoente do Mp
void criarListaComprimentos(short int *comprimentos, unsigned long p, unsigned long a){
    unsigned long d, soma = 0;

    for(unsigned long i = 1; i < a; i++){
        d = (i*p)/a + 1 - soma;
        comprimentos[i-1] = d;
        soma += d;
    }
    comprimentos[a-1] = p - soma;

}


void criarListasSigma(unsigned long *sigmas, unsigned long *inversos, unsigned long p, unsigned long q, unsigned long sigma, unsigned long a){
    unsigned long r, potencia = 1, m = invmod(a, q);
    unsigned long *aux = (unsigned long*)malloc(a*sizeof(unsigned long));

    for(unsigned long i = 0; i < a; i++){
        aux[i] = potencia;
        potencia = multmod(potencia, sigma, q);
    }
    
    sigmas[0] = 1;
    inversos[0] = m;
    for(unsigned long i = 1; i < a; i++){
        r = a - (i*p)%a;
        sigmas[i] = aux[r];
        inversos[i] = multmod(invmod(aux[r], q), m, q);
    }

    free(aux);
}


// m = profundidade da FFT
// q = 2^64 - 2^32 + 1
// p = expoente do Mp
// comprimentos[i] eh o comprimento máximo que P[i] deve ter ao final de todas as operações
void transformar_base_sigmap(unsigned long *P, unsigned long *sigmas, unsigned long q, unsigned long p, unsigned long a){
    for(unsigned long i = 1; i < a; i++){
        P[i] = multmod(P[i], sigmas[i], q);
    }
}


void voltar_base_sigmap(unsigned long *P, unsigned long *inversos, unsigned long q, unsigned long p, unsigned long a){
    for(unsigned long i = 0; i < a; i++){
        P[i] = multmod(P[i], inversos[i], q);
    }
}


void aoQuadrado(unsigned long *P, unsigned long *raizes, unsigned long *inversas, unsigned long n, unsigned long p){
    FFT(P, raizes, n, p);

    unsigned long a = powermod(2, n, p);
    for(int i = 0; i< a; i++){
        P[i] = powermod(P[i], 2, p);
    }

    FFTinversa(P, inversas, n, p);
}


// ajusta os tamanhos dos coeficientes de acordo com a lista de comprimentos
void vai_1(unsigned long *P, short int *comprimentos, unsigned int comprimento1, unsigned int comprimento2, unsigned int mask1, unsigned int mask2, unsigned long a){
    long continuar = 1;
    
    while(continuar){
        for(int i = 0; i < a-1; i++){
            if(comprimentos[i] == comprimento1){
                P[i+1] += P[i]>>comprimento1;
                P[i] = P[i]&mask1;
            } else {
                P[i+1] += P[i]>>comprimento2;
                P[i] = P[i]&mask2;
            }
        }

        if(P[a-1]>>comprimentos[a-1]){
            if(comprimentos[a-1] == comprimento1){
                P[0] += P[a-1]>>comprimento1;
                P[a-1] = P[a-1]&mask1;
            } else {
                P[0] += P[a-1]>>comprimento2;
                P[a-1] = P[a-1]&mask2;
            }
        } else {
            continuar = 0;
        }
    }
}


// p e um long,entao p <= 2^64 - 1  
int Lucas_Lehmer(unsigned long p, unsigned long m){
    unsigned long a, q, w, *S, *raizes, *inversas, *sigmas, *sigmasInversos;
    short int *comprimentos;
    q = 0xffffffff00000001;
    w = 0X0E2AF70635744FD26;   //raiz da unidade de ordem 2^32 em Zq

    unsigned long sigma = 16817914963793830057; // sigma^(2^26) = 2 em Zp
    sigma = powermod(sigma, 1<<(26-m), q); // sigma^(2^m) = 2 em Zp
    
    a = powermod(2, m, q); //numero de coeficientes de S
    w = powermod(w, powermod(2, 32-m, q), q); //agora w tem ordem 2^m
    
    S = (unsigned long*)calloc(a,sizeof(unsigned long));
    comprimentos = (short int*)malloc(a*sizeof(short int));
    raizes = (unsigned long*)malloc((a>>1)*sizeof(long));
    inversas = (unsigned long*)malloc((a>>1)*sizeof(long));
    sigmas = (unsigned long*)malloc((a)*sizeof(long));
    sigmasInversos = (unsigned long*)malloc((a)*sizeof(long));
    
    criarListaRaizes(raizes, m, w, q);
    inverterListaRaizes(inversas, raizes, m, w, q);
    criarListaComprimentos(comprimentos, p, a);
    criarListasSigma(sigmas, sigmasInversos, p, q, sigma, a);

    unsigned int comprimento1 = p/a;
    unsigned int comprimento2 = comprimento1 + 1;
    unsigned int mask1 = (1<<comprimento1) -1, mask2 = (1<<comprimento2) - 1;
    
    S[0] = 4;  //inicializacao da sequencia do teste
    
    for(int i = 0; i < p-2; i++){ //aplicação do teste
        transformar_base_sigmap(S, sigmas, q, p, a);
        aoQuadrado(S, raizes, inversas, m, q);
        voltar_base_sigmap(S, sigmasInversos, q, p, a);
        sub2novo(S, comprimentos, a);
        vai_1(S, comprimentos, comprimento1, comprimento2, mask1, mask2, a);
    }

    short int continuar = 0;
    for(int i = 0; i < a; i++){
        if(S[i] != (1<<comprimentos[i]) - 1){
            continuar++;
            break;
        }
    }
    if(continuar == 0){
        for(int i = 0; i < a; i++)
            S[i] = 0;
    }

    for(int i = 0; i < a; i++){
        if(S[i] != 0){
            free(S);
            free(raizes);
            free(inversas);
            free(sigmas);
            free(sigmasInversos);
            free(comprimentos);
            return 0;
        }
    }
    free(S);
    free(raizes);
    free(inversas);
    free(sigmas);
    free(sigmasInversos);
    free(comprimentos);
    return 1;
}




void main(){
    unsigned long p = 44497;  //44497 859433
    unsigned long m;

    m = profund(p) <= 1 ? 1 : profund(p) - 1;

    printf("m = %lu\n",m);

    printf("%i\n", Lucas_Lehmer(p, m));
}