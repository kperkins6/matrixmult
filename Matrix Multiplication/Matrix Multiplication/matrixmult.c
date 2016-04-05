//
//  matrixmult.c
//  Matrix Multiplication
//
//  Created by Kevin Perkins on 3/26/16.
//  Copyright © 2016 Kevin Perkins. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <string.h>

void Check_for_error(int local_ok, char fname[], char message[],
                     MPI_Comm comm);
void Get_dims( int* n_p, int* local_n_p,
              int my_rank, int comm_sz, MPI_Comm comm);
void Allocate_arrays(int** local_A, int** local_B,
                     int** local_C, int   int n, int local_n,
                     MPI_Comm comm);
void Read_matrix(char prompt[], int local_A[],
                 int n, int my_rank, MPI_Comm comm);
void Generate_matrix(int local_A[], int   int n);
void Print_matrix(char title[], int local_A[], int m, int
                  int n, int my_rank, MPI_Comm comm);
// void ijk_mult(int local_A[], int local_B[], int local_C[], int n,
              // int local_n, MPI_Comm comm);
// void ikj_mult(int local_A[], int local_B[], int local_C[], int n,
              // int local_n, MPI_Comm comm);
// void kij_mult(int local_A[], int local_B[], int local_C[], int n,
              // int local_n, MPI_Comm comm);

/*-------------------------------------------------------------------*/
int main(void) {
    int* local_A;
    int* local_B;
    int* local_C;
    int* C, *B, *A;
    char* flag = malloc(sizeof(char));
    int* form = malloc(sizeof(int));
    int n, local_n, a_part;
    int my_rank, comm_sz;
    MPI_Comm comm;
    int start, finish, loc_elapsed, elapsed;
    MPI_Init(NULL, NULL);
    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &comm_sz);
    MPI_Comm_rank(comm, &my_rank);
    srandom(time(NULL));

    Get_dims(&m, &n, &local_n, my_rank, comm_sz, comm);
    Allocate_arrays(&local_A, &local_B, &local_C, &C, n, local_n, comm);
    if (*flag == "I" ) {
        printf("Enter the A matrix <A>\n");
        Read_matrix("A", local_A, n, local_n, n, my_rank, comm);
        Print_matrix("A", local_A, n, local_n, n, my_rank, comm);
        printf("Enter the B matrix <B>\n");
        Read_matrix("A", local_A, n, local_n, n, my_rank, comm);
        Print_matrix("A", local_A, n, local_n, n, my_rank, comm);
    }
    else {
        Generate_matrix(local_A, n);
        Generate_matrix(local_B, n);

    }

    MPI_Barrier(comm);
//     start = MPI_Wtime();
//     if ( !strcmp(form, "ijk"))
//         ijk_mult(local_A, local_B, local_C,   n, local_n, comm);
//     else if ( !strcmp(form, "ikj"))
//         ikj_mult(local_A, local_B, local_C,   n, local_n, comm);
//     else if ( !strcmp(form, "kij"))
//         kij_mult(local_A, local_B, local_C,   n, local_n, comm);
//     else
//         printf("Invalid <form>\n");
//     finish = MPI_Wtime();
//     loc_elapsed = finish-start;
//     MPI_Reduce(&loc_elapsed, &elapsed, 1, MPI_INT, MPI_MAX, 0, comm);
//
// #  ifdef DEBUG
//     Print_vector("C", local_C, m,   my_rank, comm);
// #  endif
//
//     if (my_rank == 0)
//         printf("Elapsed time = %e\n", elapsed);
//
//     free(local_A);
//     free(local_B);
//     free(local_C);
    MPI_Finalize();
    return 0;
}  /* main */


/*-------------------------------------------------------------------*/
void Check_for_error(
                     int       local_ok   /* in */,
                     char      fname[]    /* in */,
                     char      message[]  /* in */,
                     MPI_Comm  comm       /* in */) {
    int ok;

    MPI_Allreduce(&local_ok, &ok, 1, MPI_INT, MPI_MIN, comm);
    if (ok == 0) {
        int my_rank;
        MPI_Comm_rank(comm, &my_rank);
        if (my_rank == 0) {
            fprintf(stderr, "Proc %d > In %s, %s\n", my_rank, fname,
                    message);
            fflush(stderr);
        }
        MPI_Finalize();
        exit(-1);
    }
}  /* Check_for_error */


/*-------------------------------------------------------------------*/
void Get_dims(
              int*      m_p        /* out */,
              int*      n_p        /* out */,
              int*      local_n_p  /* out */,
              char*     flag       /* out */,
              char*     form       /* out */,
              int       my_rank    /* in  */,
              int       comm_sz    /* in  */,
              MPI_Comm  comm       /* in  */) {
    int local_ok = 1;

    if (my_rank == 0) {
        printf("Enter the <form>\n");
        scanf("%s", form);
        printf("Enter the <flag>\n");
        scanf("%s", flag);
        printf("Enter the dimentions <n>\n");
        scanf("%i", &m_p);
    }
    MPI_Bcast(m_p, 1, MPI_INT, 0, comm);
    MPI_Bcast(n_p, 1, MPI_INT, 0, comm);
    if (*m_p <= 0 || *n_p <= 0 || *m_p % comm_sz != 0
        || *n_p % comm_sz != 0) local_ok = 0;
    Check_for_error(local_ok, "Get_dims",
                    "m and n must be positive and evenly divisible by comm_sz",
                    comm);

    *local_n_p = *n_p/comm_sz;
}  /* Get_dims */

/*-------------------------------------------------------------------*/
void Allocate_arrays(
                     int**  local_A  /* out */,
                     int**  local_B  /* out */,
                     int**  local_C  /* out */,
                     int       n           /* in  */,
                     int       local_n     /* in  */,
                     MPI_Comm  comm        /* in  */) {

    int local_ok = 1;

    *local_A = malloc(n*n*sizeof(int));
    *local_B = malloc(n*n*sizeof(int));
    *local_C = malloc(n*n*sizeof(int));

    if (*local_A == NULL || local_B == NULL ||
        local_C == NULL) local_ok = 0;
    Check_for_error(local_ok, "Allocate_arrays",
                    "Can't allocate local arrays", comm);
}  /* Allocate_arrays */

/*-------------------------------------------------------------------*/
void Read_matrix(
                 char      prompt[]   /* in  */,
                 int    local_A[]  /* out */,
                 int       local_n    /* in  */,
                 int       n          /* in  */,
                 int       my_rank    /* in  */,
                 MPI_Comm  comm       /* in  */) {
    int* A = NULL;
    int local_ok = 1;
    int i, j;

    if (my_rank == 0) {
        A = malloc(n*n*sizeof(int));
        if (A == NULL) local_ok = 0;
        Check_for_error(local_ok, "Read_matrix",
                        "Can't allocate temporary matrix", comm);
        printf("Enter the matrix %s\n", prompt);
        for (i = 0; i < m; i++)
            for (j = 0; j < n; j++)
                scanf   ("%lf", &A[i*n+j]);
        MPI_Scatter(A, local_n*n, MPI_INT,
                    local_A, local_n*n, MPI_INT, 0, comm);
        free(A);
    } else {
        Check_for_error(local_ok, "Read_matrix",
                        "Can't allocate temporary matrix", comm);
        MPI_Scatter(A, local_n*n, MPI_INT,
                    local_A, local_n*n, MPI_INT, 0, comm);
    }
}  /* Read_matrix */

/*-------------------------------------------------------------------*/
void Read_vector(
                 char      prompt[]     /* in  */,
                 int    local_vec[]  /* out */,
                 int       n            /* in  */,
                 int       local_n      /* in  */,
                 int       my_rank      /* in  */,
                 MPI_Comm  comm         /* in  */) {
    int* vec = NULL;
    int i, local_ok = 1;

    if (my_rank == 0) {
        vec = malloc(n*sizeof(int));
        if (vec == NULL) local_ok = 0;
        Check_for_error(local_ok, "Read_vector",
                        "Can't allocate temporary vector", comm);
        printf("Enter the vector %s\n", prompt);
        for (i = 0; i < n; i++)
            scanf("%lf", &vec[i]);
        MPI_Scatter(vec, local_n, MPI_INT,
                    local_vec, local_n, MPI_INT, 0, comm);
        free(vec);
    } else {
        Check_for_error(local_ok, "Read_vector",
                        "Can't allocate temporary vector", comm);
        MPI_Scatter(vec, local_n, MPI_INT,
                    local_vec, local_n, MPI_INT, 0, comm);
    }
}  /* Read_vector */

/*-------------------------------------------------------------------*/
void Generate_matrix(
                     int local_A[]  /* out */,
                     int local_n,   /* in */,
                     int    n          /* in  */){

    int i, j, q;
    printf("Generating Matrix, local_n=%i, n= %i\n", local_n, n);
    for (i = 0; i < local_n; i++) {
        for (j = 0; j < n; j++) {
            q = (( int) random()%42);
            local_A[i*n + j] = q;
            printf("q = %i ", q);
        }
        printf("\n");
    }
}  /* Generate_matrix */

/*-------------------------------------------------------------------*/
void Generate_vector(
                     int local_B[] /* out */,
                     int    local_n   /* in  */) {
    int i;

    for (i = 0; i < local_n; i++)
        local_B[i] = ((int) random())%42);
}  /* Generate_vector */

/*-------------------------------------------------------------------*/
void Print_matrix(
                  char      title[]    /* in */,
                  int    local_A[]  /* in */,
                  int       local_n    /* in */,
                  int       n          /* in */,
                  int       my_rank    /* in */,
                  MPI_Comm  comm       /* in */) {
    int* A = NULL;
    int i, j, local_ok = 1;

    if (my_rank == 0) {
        A = malloc(m*n*sizeof(int));
        if (A == NULL) local_ok = 0;
        Check_for_error(local_ok, "Print_matrix",
                        "Can't allocate temporary matrix", comm);
        MPI_Gather(local_A, local_n*n, MPI_INT,
                   A, local_n*n, MPI_INT, 0, comm);
        printf("\nThe matrix %s\n", title);
        for (i = 0; i < m; i++) {
            for (j = 0; j < n; j++)
                printf("%f ", A[i*n+j]);
            printf("\n");
        }
        printf("\n");
        free(A);
    } else {
        Check_for_error(local_ok, "Print_matrix",
                        "Can't allocate temporary matrix", comm);
        MPI_Gather(local_A, local_n*n, MPI_INT,
                   A, local_n*n, MPI_INT, 0, comm);
    }
}  /* Print_matrix */

/*-------------------------------------------------------------------*/
void Print_vector(
                  char      title[]     /* in */,
                  int    local_vec[] /* in */,
                  int       n           /* in */,
                  int       local_n     /* in */,
                  int       my_rank     /* in */,
                  MPI_Comm  comm        /* in */) {
    int* vec = NULL;
    int i, local_ok = 1;

    if (my_rank == 0) {
        vec = malloc(n*sizeof(int));
        if (vec == NULL) local_ok = 0;
        Check_for_error(local_ok, "Print_vector",
                        "Can't allocate temporary vector", comm);
        MPI_Gather(local_vec, local_n, MPI_INT,
                   vec, local_n, MPI_INT, 0, comm);
        printf("\nThe vector %s\n", title);
        for (i = 0; i < n; i++)
            printf("%f ", vec[i]);
        printf("\n");
        free(vec);
    }  else {
        Check_for_error(local_ok, "Print_vector",
                        "Can't allocate temporary vector", comm);
        MPI_Gather(local_vec, local_n, MPI_INT,
                   vec, local_n, MPI_INT, 0, comm);
    }
}  /* Print_vector */

/*-------------------------------------------------------------------*/
//void Mat_vect_mult(
//                   int    local_A[]  /* in  */,
//                   int    local_B[]  /* in  */,
//                   int    local_C[]  /* out */,
//                   int       n          /* in  */,
//                   int       local_n    /* in  */,
//                   MPI_Comm  comm       /* in  */) {
//    int* B;
//    int local_i, j;
//    int local_ok = 1;
//
//    B = malloc(n*sizeof(int));
//    if (B == NULL) local_ok = 0;
//    Check_for_error(local_ok, "Mat_vect_mult",
//                    "Can't allocate temporary vector", comm);
//    MPI_Allgather(local_B, local_n, MPI_INT,
//                  B, local_n, MPI_INT, comm);
//
//    for (local_i = 0; local_i < local_n; local_i++) {
//        local_C[local_i] = 0.0;
//        for (j = 0; j < n; j++)
//            local_C[local_i] += local_A[local_i*n+j]*B[j];
//    }
//    free(B);
//}  /* Mat_vect_mult */

/*


 x = malloc(n*sizeof(int));
 if (x == NULL) local_ok = 0;
    Check_for_error(local_ok, "Mat_vect_mult",
    "Can't allocate temporary vector", comm);
 MPI_Allgather(local_x, local_n, MPI_INT,
 x, local_n, MPI_INT, comm);

 for (local_i = 0; local_i < local_n; local_i++) {
    local_y[local_i] = 0.0;
    for (j = 0; j < n; j++)
        local_y[local_i] += local_A[local_i*n+j]*x[j];
}
 free(x);


*/


void ijk_mult(
              int    local_A[]  /* in  */,
              int    local_B[]  /* in  */,
              int    local_C[]  /* out */,
              int       n          /* in  */,
              int       local_n    /* in  */,
              MPI_Comm  comm       /* in  */) {

    int i, k, j;
    for (i = 0; i < *local_n; i++) {
        for(j=0; j < *n; j++) {
//            local_C[i*(*n)+j] = 0;
            for (k = 0; k < *n; k++) {
                local_C[i*(*n)+j] += local_A[i*(*n)+k] * local_B[k*(*n)+j];
            }
        }
    }
}


void ikj_mult(
              int    local_A[]  /* in  */,
              int    local_B[]  /* in  */,
              int    local_C[]  /* out */,
              int       n          /* in  */,
              int       local_n    /* in  */,
              MPI_Comm  comm       /* in  */) {

//    int* B;
//    int local_i, j;
//    int local_ok = 1;
//
//    B = malloc(n*sizeof(int));
//    if (B == NULL) local_ok = 0;
//    Check_for_error(local_ok, "Mat_vect_mult",
//                    "Can't allocate temporary vector", comm);
//    MPI_Allgather(local_B, local_n, MPI_INT,
//                  B, local_n, MPI_INT, comm);




    int i, k, j;
    for (i = 0; i < *local_n; i++) {
        for (k = 0; k < *n; k++) {
//            local_C[i*(*n)+j] = 0;
            for(j=0; j < *n; j++) {
                local_C[i*(*n)+j] += localA[i*(*n)+k] * local_B[k*(*n)+j];
            }
        }
    }
//    free(x);
}

void kij_mult(
              int    local_A[]  /* in  */,
              int    local_B[]  /* in  */,
              int    local_C[]  /* out */,
              int       n          /* in  */,
              int       local_n    /* in  */,
              MPI_Comm  comm       /* in  */) {

//    int* B;
//    int local_i, j;
//    int local_ok = 1;
//
//    B = malloc(n*sizeof(int));
//    if (B == NULL) local_ok = 0;
//    Check_for_error(local_ok, "Mat_vect_mult",
//                    "Can't allocate temporary vector", comm);
//    MPI_Allgather(local_B, local_n, MPI_INT,
//                  B, local_n, MPI_INT, comm);




    int i, k, j;
    for (k = 0; k < *n; k++) {
        for (i = 0; i < *local_n; i++) {
//            local_C[i*(*n)+j] = 0;
            for(j=0; j < *n; j++) {
                local_C[i*(*n)+j] += local_A[i*(*n)+k] * local_B[k*(*n)+j];
            }
        }
    }
//    free(x);
}
