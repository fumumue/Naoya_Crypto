#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//符号のパラーメータの指定。通常[N,K,T]として、
//Nは符号の長さ、Kが符号の次元、Tは訂正エラー数
//を表す。ここではDは符号長にしている。
#define N 257 // set small prime ex. p=2053
#define M N// puncture code
#define K (120) // degree of polynomial
#define E (9)   // bit size of prime
#define DEG 520 // set (K * E) < N
#define T (K / 2) // weight of error vector
#define Q K*2

unsigned int mat[N][K*E] = {0};
unsigned int ma[N][K*E] = {0};
