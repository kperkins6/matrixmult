//
//  matrixmult2.c
//  Matrix Multiplication
//
//  Created by Kevin Perkins on 3/28/16.
//  Copyright © 2016 Kevin Perkins. All rights reserved.
//  * Compile:  mpicc -g -Wall -o mpi_mat_vect_time mpi_mat_vect_time.c
//  * Run:      mpiexec -n <number of processes> ./mpi_mat_vect_time

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <string.h>

void Get_dims(int* n, int* local_n, char* flag, int* form, int my_rank, int comm_sz, MPI_Comm comm);
void Allocate_arrays(int** local_A_pp,  int** local_B_pp, int** local_C_pp, int my_rank, int n, int local_n, MPI_Comm comm);
void Construct_form( int local_A[],  int local_B[], int n, int my_rank, char* flag, int a_part );
void Generate_matrix(int local_A[],  int local_n, int n);
void Print_matrix(char title[],  int local_A[], int local_n, int n, int my_rank, MPI_Comm comm);
void Check_for_error(int local_ok, char fname[], char message[], MPI_Comm comm);
void ijk_mult(double local_A[], double local_B[], double local_C[], int   int n, int local_n, MPI_Comm comm);
void ikj_mult(double local_A[], double local_B[], double local_C[], int   int n, int local_n, MPI_Comm comm);
void kij_mult(double local_A[], double local_B[], double local_C[], int   int n, int local_n, MPI_Comm comm);

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

    //Don’t worry about P > n.
    Get_dims(&n, &local_n, flag,  form, my_rank, comm_sz, comm);
    //In the case of n not evenly divisible by P.
    //Give the last process all of the extra from N%P.
    Allocate_arrays(&local_A, &local_B, &local_C, my_rank, n, local_n, comm);
    a_part = n*n/comm_sz;
    srandom(time(NULL));
    C = malloc(n*n*sizeof(int));
    B = malloc(n*n*sizeof(int));
    A = malloc(n*n*sizeof(int));
    memset(A, 0, n*n*sizeof(int));
    memset(B, 0, n*n*sizeof(int));
    memset(C, 0, n*n*sizeof(int));
    MPI_Barrier(MPI_COMM_WORLD);

  if (my_rank == 0) {
    if (!strcmp(flag, "I")) {
      Construct_form(A, B, n, my_rank, flag, a_part);
      // Print_matrix("A", A, local_n, n, my_rank, comm);
      // Print_matrix("B", B, local_n, n, my_rank, comm);
    }
    else {
      printf("Generating Local(A and B)\n");
      printf("Generating A \n");
      Generate_matrix(local_A, local_n, n);
      printf("Generating B\n");
      Generate_matrix(local_B, n, n);
      // Print_matrix("local_A", local_A, local_n, n, my_rank, comm);
      // Print_matrix("local_B", local_B, local_n, n, my_rank, comm);
    }
  }
    MPI_Scatter(A, a_part, MPI_INT, local_A, a_part, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(B, n*n, MPI_INT, 0, MPI_COMM_WORLD);
    Print_matrix("A", local_A, local_n, n, my_rank, comm);
    Print_matrix("B", local_B, local_n, n, my_rank, comm);
    // free(local_A);
    // free(local_B);
    // free(local_C);
    // Read_matrix("A", local_A, m,   n, my_rank, comm);

//    MPI_Barrier(comm);
//    start = MPI_Wtime();
//    if ( !strcmp(form, "ijk"))
//        ijk_mult(local_A, local_B, local_C,   n, local_n, comm);
//    else if ( !strcmp(form, "ikj"))
//        ikj_mult(local_A, local_B, local_C,   n, local_n, comm);
//    else if ( !strcmp(form, "kij"))
//        kij_mult(local_A, local_B, local_C,   n, local_n, comm);
//    else
//        printf("Invalid <form>\n");
//    finish = MPI_Wtime();
//    loc_elapsed = finish-start;
//    MPI_Reduce(&loc_elapsed, &elapsed, 1, MPI_INT, MPI_MAX, 0, comm);
//
//#  ifdef DEBUG
//    Print_vector("C", local_C, m,   my_rank, comm);
//#  endif
//
//    if (my_rank == 0)
//        printf("Elapsed time = %e\n", elapsed);
// hello
    // free(local_A);
    // free(local_B);
    // free(local_C);
    MPI_Finalize();
    return 0;
}  /* main */

/*-------------------------------------------------------------------*/
void Get_dims(
              int*      n        /* out */,
              int*      local_n  /* out */,
              char*     flag       /* out */,
              int*      form       /* out */,
              int       my_rank    /* in  */,
              int       comm_sz    /* in  */,
              MPI_Comm  comm       /* in  */) {
    int local_ok = 1;
    char* temp = malloc(4*sizeof(char));
    if (my_rank == 0) {
        printf("Enter the <form>\n");
        scanf("%s", temp);
        printf("Enter the <flag> ('I' or 'R')\n");
        scanf("%s", flag);
        if (strcmp(temp, "ijk"))
            *form = 1;
        else if (strcmp(temp, "ikj"))
            *form = 2;
        else *form = 3;
        printf("Enter the dimentions <n>\n");
        scanf("%i", n);

    }
    MPI_Bcast(n, 1, MPI_INT, 0, comm);
    MPI_Bcast(flag, 1, MPI_CHAR, 0, comm);
    MPI_Bcast(form, 4*sizeof(char), MPI_CHAR, 0, MPI_COMM_WORLD);

    if (*n <= 0 || *n % comm_sz != 0) local_ok = 0;
    Check_for_error(local_ok, "Get_dims",
                    "n must be positive and evenly divisible by comm_sz",
                    comm);
    *local_n = *n/comm_sz;
    printf("n = %i, local_n= %i, comm_sz = %i, my_rank= %i.\n", *n, *local_n, comm_sz, my_rank);
}  /* Get_dims */

void Allocate_arrays(
                     int**  local_A_pp  /* out */,
                     int**  local_B_pp  /* out */,
                     int**  local_C_pp  /* out */,
                     int       my_rank    /* in  */,
                     int       n           /* in  */,
                     int       local_n     /* in  */,
                     MPI_Comm  comm        /* in  */) {

    int local_ok = 1;
    // if (my_rank == 0) {
        *local_A_pp = malloc(local_n*n*sizeof( int));
        *local_B_pp = malloc(n*n*sizeof( int));
        *local_C_pp = malloc(local_n*n*sizeof( int));
         memset(*local_A_pp, 0, local_n*n*sizeof( int));
         memset(*local_B_pp, 0, n*n*sizeof( int));
         memset(*local_C_pp, 0, local_n*n*sizeof( int));
    // } //
    printf("n = %i, local_n= %i, my_rank= %i.\n", n, local_n, my_rank);

    if (*local_A_pp == NULL || local_B_pp == NULL ||
        local_C_pp == NULL) local_ok = 0;
    Check_for_error(local_ok, "Allocate_arrays",
                    "Can't allocate local arrays", comm);
}  /* Allocate_arrays */

void Construct_form(
                    int local_A[],
                    int local_B[],
                    int n,
                    int my_rank,
                    char* flag,
                    int a_part ) {

        int *A = malloc(n*n*sizeof(int));
        if (!strcmp(flag, "I") ) {
            if (my_rank == 0) {
              printf("Enter the A matrix <A>, flag = %s \n", flag);
              int m_ai;
              int* m_x = malloc(sizeof( int));
              for (m_ai = 0; m_ai < (n*n); m_ai++) {
                  scanf("%i", m_x);
                  A[m_ai] = *m_x;
                }
                printf("Enter the B matrix <B>\n");
                for (m_ai = 0; m_ai < (n*n); m_ai++) {
                  scanf("%i", m_x);
                  local_B[m_ai] = *m_x;
                }
                // MPI_Scatter(A, a_part, MPI_INT, local_A, a_part, MPI_INT, 0, MPI_COMM_WORLD);
                // MPI_Bcast(local_B, n*n, MPI_INT, 0, MPI_COMM_WORLD);
                // free(m_x);
            }
        }
        else {
            // srandom(time(NULL));
            // if (my_rank == 0) {
            //   printf("Generating A and B\n");
            //   printf("Generating A \n");
            //   Generate_matrix(A, n);
            //   printf("Scattering A\n");
            //   MPI_Scatter(A, a_part, MPI_INT, local_A, a_part, MPI_INT, 0, MPI_COMM_WORLD);
            //   printf("Generating B\n");
            //   Generate_matrix(local_B, n);
            //   printf("Scattering B\n");
            //   MPI_Bcast(local_B, n*n, MPI_INT, 0, MPI_COMM_WORLD);
            // }
        }
}

void Generate_matrix(
                     int local_A[]  /* out */,
                     int local_n,   /* in */
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

void Print_matrix(
                  char      title[]    /* in */,
                  int    local_A[]     /* in */,
                  int       local_n    /* in */,
                  int       n          /* in */,
                  int       my_rank    /* in */,
                  MPI_Comm  comm       /* in */) {
    int* A = NULL;
    int i, j, local_ok = 1;
    printf("In Print\n");
    if (my_rank == 0) {
        A = malloc(n*n*sizeof( int));
        if (A == NULL) local_ok = 0;
        Check_for_error(local_ok, "Print_matrix",
                        "Can't allocate temporary matrix", comm);
        MPI_Gather(local_A, local_n*n, MPI_INT,
                   A, local_n*n, MPI_INT, 0, comm);
        printf("\nThe matrix %s\n", title);
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++)
                printf("%i ", A[i*n+j]);
            printf("\n");
        }
        printf("\n");
        // free(A);
    } else {
        Check_for_error(local_ok, "Print_matrix",
                        "Can't allocate temporary matrix", comm);

        printf("Gathering from %i\n", my_rank);
        MPI_Gather(local_A, local_n*n, MPI_INT,
                   A, local_n*n, MPI_INT, 0, comm);
    }
}  /* Print_matrix */

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

void ijk_mult(
              double    local_A[]  /* in  */,
              double    local_B[]  /* in  */,
              double    local_C[]  /* out */,
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
              double    local_A[]  /* in  */,
              double    local_B[]  /* in  */,
              double    local_C[]  /* out */,
              int       n          /* in  */,
              int       local_n    /* in  */,
              MPI_Comm  comm       /* in  */) {

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
              double    local_A[]  /* in  */,
              double    local_B[]  /* in  */,
              double    local_C[]  /* out */,
              int       n          /* in  */,
              int       local_n    /* in  */,
              MPI_Comm  comm       /* in  */) {

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
