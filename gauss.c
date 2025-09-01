# include <stdio.h>


# define EN  3   //the dimension of equation

void vec_diff(int a[EN], int b[EN]){
    /* Calcurate the difference of two vectors. Be caution that b[N] changes.*/
    for (int i = 0; i < EN; i++){
        b[i] = (N+(b[i]-a[i])%N)%N;
    }
}

int gauss(){
    int m[EN][EN] = {{5,N-1,N-1}, {2,1,N-3},{1,1,1}};    // The matrix
    int b[EN] = {0,N-5,6};


    printf("The coefficient matrix is : \n");
    for (int i = 0; i < EN; i++){
        for (int j = 0; j < EN; j++){
            printf("%d ", m[i][j]);
            if (j == EN-1){
                printf("\n");
            }
        }
    }

    printf("\nUse Gauss method to solve equations : \n");
    for (int i = 0; i < EN; i++){
        for (int j = i+1; j < EN; j++){
            int coef = m[j][i] * inv(m[i][i],N)%N;
            int del[EN];

            for (int k = 0; k < EN; k++){
                del[k] = m[i][k] * coef%N;
            }
            vec_diff(del, m[j]);
            b[j] = (N-b[j]%N-b[i] * coef%N)%N;
        }
    }

    for (int i = EN -1; i >= 0; i--){
        int x = 1*inv(m[i][i],N)%N;
        m[i][i] = m[i][i]*x%N;
        b[i] = b[i]*x%N;

        for (int j = 0; j < i; j++){
            b[j] = (N+b[j]%N-b[i]*m[j][i]%N)%N;
            m[j][i] = 0;
        }
    }

    for (int i = 0; i < EN; i++){
        for (int j = 0; j < EN; j++){
            printf("%d ", (N+m[i][j]%257)%N);
            if (j == EN - 1){
                printf("\n");
            }
        }
    }

    for (int i = 0; i < EN; i++){
        printf("%d ", (N+b[i]%257)%N);
    }
    printf("%d\n",N);

    return 0;
}
