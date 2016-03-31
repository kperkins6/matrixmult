//
//  matrixmult2.c
//  Matrix Multiplication
//
//  Created by Kevin Perkins on 3/28/16.
//  Copyright © 2016 Kevin Perkins. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <string.h>

void Get_dims( int* n_p, int* local_n_p, char* flag[], char* form,
                     int my_rank, int comm_sz, MPI_Comm comm);
void Allocate_arrays( int** local_A_pp,  int** local_B_pp,
                      int** local_C_pp, int my_rank, int n, int local_n,
                     MPI_Comm comm);
void Construct_form(  int local_A[],  int local_B[], int n,
                     char* flag );
void Generate_matrix( int local_A[], int n);
void Print_matrix(char title[],  int local_A[], int local_n,
                  int n, int my_rank, MPI_Comm comm);
void Check_for_error(int local_ok, char fname[], char message[],
                     MPI_Comm comm);


int main(void) {
     int* local_A;
     int* local_B;
     int* local_C;
     int* C;
    char* flag;
    char* form[3];
    int n, local_n;
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
    Construct_form(local_A, local_B, n, flag);
    Print_matrix("A", local_A, local_n, local_n, n, my_rank, comm);
    Print_matrix("B", local_B, local_n,  local_n, n, my_rank, comm);
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
//    
    free(local_A);
    free(local_B);
    free(local_C);
    MPI_Finalize();
    return 0;
}  /* main */

/*-------------------------------------------------------------------*/
void Get_dims(
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
        scanf("%i", n_p);
    }
    MPI_Bcast(n_p, 1, MPI_INT, 0, comm);
    if (*n_p <= 0 || *n_p % comm_sz != 0) local_ok = 0;
    Check_for_error(local_ok, "Get_dims",
                    "n must be positive and evenly divisible by comm_sz",
                    comm);
    
    *local_n_p = *n_p/comm_sz;
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
    if (my_rank == 0) {
        *local_A_pp = malloc(n*n*sizeof( int));
        *local_B_pp = malloc(n*n*sizeof( int));
        *local_C_pp = malloc(n*n*sizeof( int));
        memset(local_A_pp, 0, n*n*sizeof( int));
        memset(local_B_pp, 0, n*n*sizeof( int));
        memset(local_C_pp, 0, n*n*sizeof( int));
    }

    if (*local_A_pp == NULL || local_B_pp == NULL ||
        local_C_pp == NULL) local_ok = 0;
    Check_for_error(local_ok, "Allocate_arrays",
                    "Can't allocate local arrays", comm);
}  /* Allocate_arrays */

void Construct_form(
                     int local_A[],
                     int local_B[],
                    int n,
                    char* flag ) {
    if (my_rank == 0) {
        if (strcmp(flag, "I") ) {
            printf("Enter the A matrix <A>\n");
            int m_ai, *m_x;
            for (m_ai = 0; m_ai < (n*n); m_ai++) {
                scanf("%i", m_x);
                local_A[m_ai] = *m_x;
            }
            printf("Enter the B matrix <B>\n");
            for (m_ai = 0; m_ai < (n*n); m_ai++) {
                scanf("%i", m_x);
                local_B[m_ai] = *m_x;
            }
        }
        else {
            srandom(time(NULL));
            printf("Generating A and B\n");
            Generate_matrix(local_A, n);
            Generate_matrix(local_B, n);
            
        }
    }
}

void Generate_matrix(
                      int local_A[]  /* out */,
                     int    n          /* in  */){
    
    int i, j;
    
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            local_A[i*n + j] = (( int) random())/(( int) RAND_MAX);
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
    
    if (my_rank == 0) {
        A = malloc(n*n*sizeof( int));
        if (A == NULL) local_ok = 0;
        Check_for_error(local_ok, "Print_matrix",
                        "Can't allocate temporary matrix", comm);
        MPI_Gather(local_A, local_n*n, MPI_INT,
                   A, local_n*n, MPI_INT, 0, comm);
        printf("\nThe matrix %s\n", title);
        for (i = 0; i < m; i++) {
            for (j = 0; j < n; j++)
                printf("%i ", A[i*n+j]);
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