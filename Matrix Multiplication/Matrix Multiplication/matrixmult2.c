//
//  matrixmult2.c
//  Matrix Multiplication
//
//  Created by Kevin Perkins on 3/28/16.
//  Copyright © 2016 Kevin Perkins. All rights reserved.
//  * Compile:  mpicc -g -Wall -o matrixmult matrixmult.c
//  * Run:      mpirun -n <number of processes> ./matrixmult

// Notes:
//
// Issues:
//  The case of n/comm_sz != 0 is not handeled. the residual size will NOT be
//  acounted for.

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <string.h>

void get_format(int* n, int* local_n, char* flag, int* form, int my_rank, int comm_sz, MPI_Comm comm);
void Allocate_memory(int** local_A,  int** local_B, int** local_C, int my_rank, int n, int local_n, MPI_Comm comm);
void Input_matrices(int local_A[],  int local_B[], int local_n, int n, int my_rank, char* flag, int a_part );
void Generate_matrix(int local_A[],  int local_n, int n);
void Print_matrix(char title[],  int matrix[], int local_n, int n);
void Check_for_error(int local_ok, char fname[], char message[], MPI_Comm comm);
void ijk_mult(int local_A[], int local_B[], int local_C[], int n, int local_n, MPI_Comm comm);
void ikj_mult(int local_A[], int local_B[], int local_C[], int n, int local_n, MPI_Comm comm);
void kij_mult(int local_A[], int local_B[], int local_C[], int n, int local_n, MPI_Comm comm);

int main(void) {
    int* local_A;
    int* local_B;
    int* local_C;
    int* C;
    int* B;
    int* A;
    char* flag = malloc(sizeof(char));
    int* form = malloc(sizeof(int));
    int n, local_n, a_part;
    int my_rank, comm_sz;
    MPI_Comm comm;
    double start, finish, loc_elapsed, elapsed;

    MPI_Init(NULL, NULL);
    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &comm_sz);
    MPI_Comm_rank(comm, &my_rank);

    //Don’t worry about P > n.
    get_format(&n, &local_n, flag,  form, my_rank, comm_sz, comm);
    //In the case of n not evenly divisible by P.
    //Give the last process all of the extra from N%P.
    Allocate_memory(&local_A, &local_B, &local_C, my_rank, n, local_n, comm);
    a_part = (n*n)/comm_sz;
    // printf("a_part = %i\n", a_part);
    srandom(time(NULL));
    B = malloc(n*n*sizeof(int));
    memset(B, 0, n*n*sizeof(int));
    MPI_Barrier(MPI_COMM_WORLD);

  if (my_rank == 0) {
    A = malloc(n*n*sizeof(int));
    C = malloc(n*n*sizeof(int));
    memset(C, 0, n*n*sizeof(int));
    memset(A, 0, n*n*sizeof(int));
    if (!strcmp(flag, "I")) {
      Input_matrices(A, B, local_n, n, my_rank, flag, a_part);
      // Print_matrix("A", A, local_n, n, my_rank, comm);
      // Print_matrix("B", B, n, n);
    }
    else {
      // printf("Generating Local(A and B)\n");
      // printf("Generating A \n");
      Generate_matrix(A, n, n);
      // printf("Generating B\n");
      Generate_matrix(B, n, n);
      // Print_matrix("local_A", local_A, local_n, n, my_rank, comm);
      // Print_matrix("local_B", local_B, local_n, n, my_rank, comm);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(B, n*n, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Scatter(A, a_part, MPI_INT, local_A, a_part, MPI_INT, 0, MPI_COMM_WORLD);

   MPI_Barrier(comm);
   start = MPI_Wtime();
  //  printf( "Form = %i\n", *form);
   if ( *form == 1) {
       ijk_mult(local_A, B, local_C, n, local_n, comm);
      //  printf("IJK Form\n");
    }
   else if ( *form == 2) {
       ikj_mult(local_A, B, local_C, n, local_n, comm);
      //  printf("IKJ Form\n");
    }
   else if ( *form == 3) {
       kij_mult(local_A, B, local_C, n, local_n, comm);
      //  printf("KIJ Form\n");
    }
   else
       printf("Invalid <form>\n");
   MPI_Gather(local_C, local_n*n, MPI_INT, C, local_n*n, MPI_INT, 0, comm);
   finish = MPI_Wtime();
   loc_elapsed = finish-start;
   MPI_Reduce(&loc_elapsed, &elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, comm);

  //  if (my_rank == 0) {
  //    Print_matrix("matrix A", A, n, n);
  //    Print_matrix("matrix B", B, n, n);
  //  }
  //  Print_matrix("matrix local_A", local_A, local_n, n);
  //  Print_matrix("matrix local_C", local_C, local_n, n);
   MPI_Barrier(MPI_COMM_WORLD);

   if (my_rank == 0) {
      // printf("\n\n--------RESULT MATRIX-------\n");
      // Print_matrix("matrix C", C, n, n);
      printf("\n\n--------RESULT TIMING-------\n");
      if (comm_sz > 1)
        printf("Running on %i processors\n", comm_sz);
      else  printf("Running on 1 processor\n");
      printf("Elapsed Time = %f\n", elapsed);
      printf("\n\n\n\n");
    }
    // free(local_A);
    // free(local_B);
    // free(local_C);
    MPI_Finalize();
    return 0;
}
/* ---------------------------
  Function: Get Format
  Arguments:
    int*      n          output,
    int*      local_n    output,
    char*     flag       output,
    int*      form       output,
    int       my_rank    input,
    int       comm_sz    input,
    MPI_Comm  comm       input
  Functionality:
    Accepts setup parameters for
    the program, including flag,
    form, and dimensions
*/
void get_format(int*  n, int*  local_n, char* flag, int*  form, int my_rank,
                int comm_sz, MPI_Comm comm) {

    int local_ok = 1;
    char* temp = malloc(4*sizeof(char));
    if (my_rank == 0) {
        printf("Enter the <form>\n");
        scanf("%s", temp);
        printf("Enter the <flag> ('I' or 'R')\n");
        scanf("%s", flag);
        if (!strcmp(temp, "ijk"))
            *form = 1;
        else if (!strcmp(temp, "ikj"))
            *form = 2;
        else *form = 3;
        printf("Enter the dimentions <n>\n");
        scanf("%i", n);

    }
    MPI_Bcast(n, 1, MPI_INT, 0, comm);
    MPI_Bcast(flag, 1, MPI_CHAR, 0, comm);
    MPI_Bcast(form, 4*sizeof(char), MPI_CHAR, 0, MPI_COMM_WORLD);

    if (*n <= 0 || *n % comm_sz != 0) local_ok = 0;
    Check_for_error(local_ok, "get_format",
                    "n must be positive and evenly divisible by comm_sz",
                    comm);
    *local_n = *n/comm_sz;
    // printf("n = %i, local_n= %i, comm_sz = %i, my_rank= %i.\n", *n, *local_n, comm_sz, my_rank);
}

/* ---------------------------
  Function: Allocate Memory
  Arguments:
    int**     local_A    out,
    int**     local_B    out,
    int**     local_C    out,
    int       my_rank    in,
    int       n          in,
    int       local_n    in,
    MPI_Comm  comm       in
  Functionality:
    Allocates the necessary
    memory for each matrix.
    An error is raised if
    memory cannot be allocated
*/
void Allocate_memory( int** local_A, int** local_B, int** local_C, int my_rank,
                      int n, int local_n, MPI_Comm comm) {

    int local_ok = 1;
        *local_A = malloc(local_n*n*sizeof( int));
        *local_B = malloc(n*n*sizeof( int));
        *local_C = malloc(local_n*n*sizeof( int));
         memset(*local_A, 0, local_n*n*sizeof( int));
         memset(*local_B, 0, n*n*sizeof( int));
         memset(*local_C, 0, local_n*n*sizeof( int));
    // printf("n = %i, local_n= %i, my_rank= %i.\n", n, local_n, my_rank);

    if (*local_A == NULL || local_B == NULL ||
        local_C == NULL) local_ok = 0;
    Check_for_error(local_ok, "Allocate_memory",
                    "Can't allocate local arrays", comm);
}

/* ---------------------------
  Function: Input Matrices
  Arguments:
    int local_A[]   output,
    int local_B[]   output,
    int local_n     input,
    int n           input,
    int my_rank     input,
    char* flag      input,
    int a_part      input
  Functionality:
    Reads in the numeric
    values for the A and B
    matrices. Only works on
    thread rank 0.
*/
void Input_matrices(
                    int local_A[],
                    int local_B[],
                    int local_n,
                    int n,
                    int my_rank,
                    char* flag,
                    int a_part ) {

        if (!strcmp(flag, "I") ) {
            if (my_rank == 0) {
              printf("Enter the A matrix <A>\n");
              int m_ai;
              int* m_x = malloc(sizeof( int));
              for (m_ai = 0; m_ai < (n*n); m_ai++) {
                  scanf("%i", m_x);
                  local_A[m_ai] = *m_x;
                }
                printf("Enter the B matrix <B>");
                for (m_ai = 0; m_ai < (n*n); m_ai++) {
                  scanf("%i", m_x);
                  local_B[m_ai] = *m_x;
                }
            }
        }
        // else {
        //     // srandom(time(NULL));
        //     // if (my_rank == 0) {
        //     //   printf("Generating A and B\n");
        //     //   printf("Generating A \n");
        //     //   Generate_matrix(A, n);
        //     //   printf("Scattering A\n");
        //     //   MPI_Scatter(A, a_part, MPI_INT, local_A, a_part, MPI_INT, 0, MPI_COMM_WORLD);
        //     //   printf("Generating B\n");
        //     //   Generate_matrix(local_B, n);
        //     //   printf("Scattering B\n");
        //     //   MPI_Bcast(local_B, n*n, MPI_INT, 0, MPI_COMM_WORLD);
        //     // }
        // }
}

/* ---------------------------
  Function: Generate Matrix
  Arguments:
    int       local_A[]  out,
    int       local_n    in,
    int       n          in
  Functionality:
    Generates a random number
    matrix based on the size
    parameters local_n, n
    It is filled with random
    integers.
*/
void Generate_matrix( int local_A[], int local_n, int n){

    int i, j, q;
    // fprintf(stderr, "Generating the Matrix\n");
    for (i = 0; i < local_n; i++) {
        for (j = 0; j < n; j++) {
            q = (( int) random()%42);
            local_A[i*n + j] = q;
            // fprintf(stderr, "q = %i ", q);
        }
        // fprintf(stderr, "\n");
    }
}

void Print_matrix( char title[], int matrix[], int local_n, int n) {
  int i, j;

  fprintf(stderr, "\nPrinting %s\n", title);
  for (i = 0; i < local_n; i++) {
     for (j = 0; j < n; j++)
        fprintf(stderr, "%i ", matrix[i*n+j]);
     fprintf(stderr, "\n");
  }
}

/*-------------------------------------------------------------------*/
/*
  THIS ERROR FUNCTION WAS RETRIEVED FROM PACHECO'S CODE.
  SINCE IT'S PURELY FOR TESTING, AND ADDS NO FUNCTIONALITY
  TO THE PROGRAM, I COPIED IT DIRECTLY
*/
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

/* -------------------------------------------
Form 'ijk'
  int       local_A[]   input
  int       local_B[]   input
  int       local_C[]   output
  int       n           input
  int       local_n     input
  MPI_Comm  comm        input
Functionality:
  This function performs matrix multiplication
  on a section of a (Called local_A). It
  multiplies this section with the appropriate
  section of the entirety of matrix B, and puts
  the result into local_C
*/
void ijk_mult(int local_A[], int local_B[], int local_C[], int n, int local_n,
              MPI_Comm comm) {

    int i, k, j, q;
    // printf("In IJK\n");
    for (i = 0; i < local_n; i++) {
      // printf("\ni= %i", i);
        for(j=0; j < n; j++) {
//            local_C[i*(*n)+j] = 0;
            // printf("\n");
            for (k = 0; k < n; k++) {
                q = local_A[i*(n)+k] * local_B[k*(n)+j];
                // printf("[%i*%i=]%i",local_A[i*(n)+k], local_B[k*(n)+j], q);
                local_C[i*(n)+j] += q;
                // if (k==(n-1))
                //   printf("=%i\n", local_C[i*n+j]);
                // else printf("+");
            }
        }
    }
}

/* -------------------------------------------
Form 'ijk'
  int       local_A[]   input
  int       local_B[]   input
  int       local_C[]   output
  int       n           input
  int       local_n     input
  MPI_Comm  comm        input
Functionality:
  This function performs matrix multiplication
  on a section of a (Called local_A). It
  multiplies this section with the appropriate
  section of the entirety of matrix B, and puts
  the result into local_C
*/
void ikj_mult(int local_A[], int local_B[], int local_C[], int n, int local_n,
              MPI_Comm comm) {

    int i, k, j;
    for (i = 0; i < local_n; i++) {
        for (k = 0; k < n; k++) {
//            local_C[i*(*n)+j] = 0;
            for(j=0; j < n; j++) {
                local_C[i*(n)+j] += local_A[i*(n)+k] * local_B[k*(n)+j];
            }
        }
    }
//    free(x);
}

/* -------------------------------------------
Form 'kij'
  int       local_A[]   input
  int       local_B[]   input
  int       local_C[]   output
  int       n           input
  int       local_n     input
  MPI_Comm  comm        input
Functionality:
  This function performs matrix multiplication
  on a section of a (Called local_A). It
  multiplies this section with the appropriate
  section of the entirety of matrix B, and puts
  the result into local_C
*/
void kij_mult(int local_A[], int local_B[], int local_C[], int n, int local_n,
              MPI_Comm comm) {

    int i, k, j;
    for (k = 0; k < n; k++) {
        for (i = 0; i < local_n; i++) {
//            local_C[i*(*n)+j] = 0;
            for(j=0; j < n; j++) {
                local_C[i*(n)+j] += local_A[i*(n)+k] * local_B[k*(n)+j];
            }
        }
    }
//    free(x);
}
