#include <stdio.h>
#include <stdlib.h>

/* $BAH9g$;$N@8@.(B: $B%l%8%9%?$r;H$C$?9bB.HG!#(B */

#define NN 23
#define KK 3


/* [23,12,7]-�g���S�[���C���� 
unsigned int ht[11]={
0b10000000000111110010010,
0b01000000000011111001001,
0b00100000000110001110110,
0b00010000000011000111011,
0b00001000000110010001111,
0b00000100000100111010101,
0b00000010000101101111000,
0b00000001000010110111100,
0b00000000100001011011110,
0b00000000010000101101111,
0b00000000001111100100101
};
*/

unsigned int count,count2;
typedef unsigned int seti;
#define first(n) ((seti) ((1U << (n)) - 1U))
typedef struct {
  int x[4096];
} tri;

unsigned int u[32],ss[256];

unsigned int model;

unsigned sindy[2048]={0};

static unsigned short h2[8]={32884,16562,8424,4235,2135,1054,717,303};
/*
0b1000000001110100,
0b0100000010110010,
0b0010000011101000,
0b0001000010001011,
0b0000100001010111,
0b0000010000011110,
0b0000001011001101,
0b0000000100101111
};
*/


unsigned int enc(unsigned int y, unsigned int z)
{
  unsigned int c;
  
  c=0;
  while(y!=0){
    if(y&1) c ^=z;
    z<<=1; y>>=1;
  }
  return c;
  
}

int sind(int c){
  int i;
  int ss=0;
      for(i=0;i<11;i++){
      ss<<=1;
    ss^=__builtin_popcount(c&ht[i])%2;
    printf("%d,",__builtin_popcount(c&ht[i])%2);
    }
return ss;
}

/*
(23,12,7)-Golay Code generator function is
	G(x)=x^11+x^9+x^7+x^6+x^5+x+1.  
  101011100011
*/
int codec(int mm){
  int c=0b101011100011;

  return enc(c,mm);
}

int vx[2048]={0};
seti nextset(seti x)
{
  seti smallest,ripple,n_small,ones;
  
  smallest=x& -x;
  ripple=x+smallest;
  n_small=ripple & -ripple;
  ones=((n_small/smallest)>>1)-1;
  return ripple|ones;
}

tri ple={0};
int uu=0;
tri printest(seti s,int t1)
{
  int i,l,k=0;
  unsigned int b=0,j;
  char t[24];
  tri tt={0};
  unsigned short tmp=0;
  unsigned vv=0;

  for(i=0;i<NN;i++){
    if(k==t1)
      k=0;
    if(s&1){
      printf("%d ",i);
      //tmp^=gt[i];
      ple.x[uu]^=(1<<i);
      //tmp^=ht[i];
      //vv^=(1<<(i));
    }
    s >>= 1;
  }
  printf("vv=%b\n",vv);
printf("\n");

  return ple;
}


void gappa(){
int i,j,k=0;

for(j=0;j<23;j++){
  k=0;
  for(i=0;i<11;i++){
  printf("%d",(ht[i]%2));
  k<<=1;
  k^=ht[i]&(1);
  ht[i]=ht[i]>>1;
  }
printf("%d=%d %b\n",j,k,k);
printf("\n");
}
}

int fugo(void){
  int i,j;
  seti x;
  int t1;


  count=0;count2=0,i=1;
for(t1=1;t1<4;t1++)
{
 x=first(t1);
  while(! (x & ~first(N))){
    printf("%4d:",i); 
    tri vx=printest(x,t1);
    x=nextset(x); i++;
    uu++;
  }
printf("uu=%d\n",uu);
}
printf("%b\n",2047);

//exit(1);
int e=0;
//exit(1);

    int uu=0;
    int ss=0;
    int table[2048]={0};
    for(int j=0;j<2048;j++)
    {
        uu=0;
        ss=0;
    for(i=0;i<11;i++){
      ss<<=1;
    ss^=__builtin_popcount(ple.x[j]&ht[i])%2;
    printf("%d,",__builtin_popcount(ple.x[j]&ht[i])%2);
    }
    if(sindy[ss]==0){
    sindy[ss]=ple.x[j];
    }else{
      printf("baka\n");
      exit(1);
    }
    printf("\n");
    }
  
    for(i=0;i<2048;i++)
    printf("%d=%b\n",i,sindy[i]);
    //exit(1);


int mm[4096]={0};
for(i=0;i<4096;i++){
  mm[i]=codec(i);
  printf("%d=%b\n",i,mm[i]);
}
int c=codec(517)^0b10000000000000100000001;
//ss=0;
  ss=sind(c);
    printf("%b\n",sindy[ss]);
    c^=sindy[ss];
 for(i=0;i<4096;i++){
 if(mm[i]==c){
  printf("m=%d %d\n",i,c);
  break;
 }
}
//gappa();
exit(1);


return 0;
}


/*
int itob(int n,char s[])
{
  int i,j,k=0;
  
  for(i=(N-1),j=0;i>=0;i--,j++){
    s[j]=((n>>i) & 0x0001) + '0';
  }
  
  for(j=0;j<N;j++){
    if(s[j]-48==1)
      k++;
  }
  return k;
}
*/
/*
  model==strtoul("1001001001001000",(char **)NULL,2);
  scanf("%u",&N);
  scanf("%u",&K);
  for(i=0;i<N;i++)
    u[i]=(1<<i);
  for(i=0;i<256;i++)
    ss[i]=0;
*/
 /* printf("%u %s\n",model,itoa(model,2)); */
