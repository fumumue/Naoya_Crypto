#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>


typedef struct {
   unsigned x[23];
} vec;

vec conv(vec a,vec b,int n){
    int i,l=0,j=0;
    vec c={0};

    for(int k=0;k<n;k++){
        l=0;
        for(i=0;i<n;i++){
            for(j=0;j<n;j++){
            if((i+j)==k && a.x[i]<2 && b.x[j]<2){
            l+=(a.x[i]^b.x[j]);
            printf("i=%d,j=%d, %d %d\n",i,j,a.x[i],b.x[j]);
            }
            }
        }
        printf("\n");
        c.x[k]=l%2;
    }

    return c;
}

vec and(vec a,vec b){
    int i;
    vec c={0};

    for(i=0;i<23;i++)
    c.x[i]=a.x[i]&b.x[i];

    return c;
}

vec or(vec a,vec b){
    int i;
    vec c={0};

    for(i=0;i<23;i++)
    c.x[i]=(a.x[i]|b.x[i]);

return c;
}

vec xor(vec a,vec b){
    int i=0;
    vec c={0};

    for(i=0;i<23;i++)
    c.x[i]=a.x[i]^b.x[i];

    return c;
}

int wt(vec a){
    int i=0,k=0;

    for(i=0;i<23;i++){
    if(a.x[i]>0)
    k++;
    }

    return k;
}


void main(void){
    int k=11,i,count[4]={0},l=0;
    vec x={0},h={1},r1={0},r2={0},y={0};
    vec a={1,2,3};
    vec b={4,5,6};

    vec v=conv(a,b,3);
    for(i=0;i<3;i++)
    printf("%d,",v.x[i]);
    printf("\n");
    //exit(1);


    srand(clock());
    
    while(1){
    memset(x.x,0,sizeof(vec));
    memset(h.x,1,sizeof(23));
    memset(r1.x,0,sizeof(23));
    memset(r2.x,0,sizeof(23));
    memset(y.x,0,sizeof(23));

    count[2]++;
    x.x[rand()%23]=1;
    r1.x[rand()%23]=1;
    vec e={0};
    
    e.x[rand()%23]=1;

    while(count[0]<k){
        l=rand()%23;
        if(y.x[l]==0){
        y.x[l]=1;
        count[0]++;
        }
        printf("a");
    }
    while(count[1]<k){
        l=rand()%23;
        if(r2.x[l]==0){
        r2.x[l]=1;
        count[1]++;
        }
        printf("b");
    }
    printf("\n");

    vec s={0},rr=xor(r1,or(r2,h));
    s=xor(x,or(h,y));
    vec t={0};
    t=xor(e,xor(or(s,r2),or(rr,y)));
    vec u=xor(r1,or(h,r2)),vv=xor(or(s,r1),e);

    int n=wt(t),m=0,o=wt(or(or(h,y),r2));

    if(n<3 && o>3 && wt(u)>3 && wt(vv)>3)
    {
        printf("n=%d o=%d u=%d,v=%d\n",n,o,wt(u),wt(vv));
        printf("%d %b %d\n",n,m,count[2]);
        for(i=0;i<23;i++){
            m=(m<<1);
            m^=t.x[i];
        }
        exit(1);
    }
        n=0;
        m=0;
    }


    return;
}
