
/* -*- mode: C; coding:utf-8 -*- */

//monomial
typedef struct
{
   unsigned short n; //単項式の次数
   unsigned short a; //単項式の係数
} oterm;

//polynomial
typedef struct
{
  oterm t[DEG]; //単項式の配列として多項式を表現する
} OP;

typedef struct 
{
   unsigned short x[DEG]; //配列の添字を次数に、配列の値を係数に持つ多項式の表現
} vec;

typedef struct {
  vec q;
  vec r;
} rem;



#define I8T char
#define U8C(v) (v##U)

#define U8V(v) (( char)(v)&U8C(0xFF))
#define ROTL8(v, n) \
  (U8V((v) << (n)) | ((v) >> (8 - (n))))

#define R(x, n) (((x) << (n)) | ((x) >> (32 - (n))))


typedef struct {
   short x[N][N];
   vec f;
   int row; //行
   int col; //列
} MTX;

