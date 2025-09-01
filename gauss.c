# include <stdio.h>
#include <assert.h>
#include <time.h>
#include <stdlib.h>


//# define nn  3   //the dimnnsion of equation

//#include "hqc_golay.c"
void vec_diff(int a[N], int b[N],int nn){
    /* Calcurate the differnnce of two vectors. Be caution that b[N] changes.*/
    for (int i = 0; i < nn; i++){
        b[i] = (N+b[i]-a[i])%N;
    }
}


int sankaku(MTX m,int nn){
    int mm[3][4] = {{5,1,1,1}, {2,1,3,2},{1,1,1,3}};    // The matrix
    int b[3] = {0,5,6};
    int i,j;

    
    printf("The coefficinnt matrix is : \n");
    for (int i = 0; i < nn; i++){
        for (int j = 0; j < nn+1; j++){
            printf("%d ", m.x[i][j]);
            if (j == nn){
                printf("\n");
            }
        }
    }

    printf("\nUse Gauss method to solve equations : \n");
    //上三角
    for (int i = 0; i < nn; i++){
        for (int j = i+1; j < nn; j++){
            int coef = m.x[j][i] * inv(m.x[i][i],N);
            int del[nn];

            for (int k = 0; k < nn; k++){
                del[k] = m.x[i][k] * coef%N;
            }
            vec_diff(del, m.x[j],nn);
            b[j] = (N+b[j]-b[i] * coef)%N;
        }
    }
    
    for (int i = 0; i <nn; i++){
        int x =  inv(m.x[i][i],N);
        m.x[i][i] = (m.x[i][i]*x)%N;
        b[i] = b[i]*x%N;
        for(j=i+1;j<nn+1;j++)
        m.x[i][j]=m.x[i][j]*x%N;
        /*
        //下三角
        for (int j = 0; j < i; j++){
            b[i] -= b[j]*m[j][i];
            m[j][i] = 0;
        }
            */
    }
    
    

    for (int i = 0; i < nn; i++){
        for (int j = 0; j < nn+1; j++){
            printf("%d ", m.x[i][j]);
            if (j == nn){
                printf("\n");
            }
        }
    }
    printf("\n");


    return 0;
}

void vdiff(float a[3], float b[3],int nn){
    /* Calcurate the difference of two vectors. Be caution that b[N] changes.*/
    for (int i = 0; i < nn; i++){
        b[i] -= a[i];
    }
}
int gauss(int nn){
    float m[3][3] = {{5,-1,-1}, {2,1,-3},{1,1,1}};    // The matrix
    float b[3] = {0,-5,6};


    printf("The coefficient matrix is : \n");
    for (int i = 0; i < nn; i++){
        for (int j = 0; j < nn; j++){
            printf("%1.f ", m[i][j]);
            if (j == nn-1){
                printf("\n");
            }
        }
    }

    printf("\nUse Gauss method to solve equations : \n");
    for (int i = 0; i < nn; i++){
        for (int j = i+1; j < nn; j++){
            float coef = m[j][i] / m[i][i];
            float del[N];

            for (int k = 0; k < nn; k++){
                del[k] = m[i][k] * coef;
            }
            vdiff(del, m[j],nn);
            b[j] -= b[i] * coef;
        }
    }

    for (int i = nn -1; i >= 0; i--){
        float x = 1. / m[i][i];
        m[i][i] *= x;
        b[i] *= x;

        for (int j = 0; j < i; j++){
            b[j] -= b[i]*m[j][i];
            m[j][i] = 0;
        }
    }

    for (int i = 0; i < nn; i++){
        for (int j = 0; j < nn; j++){
            printf("%1.f ", m[i][j]);
            if (j == nn - 1){
                printf("\n");
            }
        }
    }

    for (int i = 0; i < nn; i++){
        printf("%f ", b[i]);
    }

    return 0;
}
