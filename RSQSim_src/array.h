#ifndef ARRAY_H
#define ARRAY_H

char **carray2(int m, int n);
void cfree2(char **A, int m, int n);

double ****array4(int k, int l, int m, int n);
double ***array3(int l, int m, int n);
double **darray2(int m, int n);
double ***darray3(int l, int m, int n);
void dfree2(double **A, int l, int m);
void dfree3(double **A, int l, int m, int n);

float **farray2(int m, int n);
float ***farray3(int l, int m, int n);
float ****farray4(int k, int l, int m, int n);
float *****farray5(int l, int m, int n, int o, int p);
void ffree2(float **A, int m, int n);
void ffree3(float ***A, int l, int m, int n);
void ffree4(float ****A, int k, int l, int m, int n);

long **larray2(int m, int n);
long ***larray3(int l, int m, int n);

short **sarray2(int m, int n);

unsigned char **ucarray2(int m, int n);

int **iarray2(int m, int n);
int **i_ragged_array2(int l, int *m);
void i_ragged_free2(int **A, int l, int *m);
void i_ragged_free_nonblock2(int **A, int l, int *m);
int ***iarray3(int l, int m, int n);
void ifree2(int **A, int m, int n);

char **struct_array2(int m, int n, int size);
char **struct_ragged_array2(int l, int *m, int size);
char **read_struct_ragged_array2(char *fname, int *l, int **m, int size);
void write_struct_ragged_array2(char *fname, int l, int *m, int size, char *ptr);
char ***struct_array3(int l, int m, int n, int size);
char ***struct_ragged_array3(int l, int m, int **n, int size);
char ****struct_array4(int k, int l, int m, int n, int size);
void struct_free2(char **A, int m, int n);

void **varray2(int m, int n, int size);
void vfree2(void **A, int m, int n);
void *****varray5(int l, int m, int n, int o, int p, int size);

int ****iarray4(int k, int l, int m, int n);
void ifree4(int ****A, int k, int l, int m, int n);

double ***gsharray3(int lmax);

float **farray2SLT(int n);
void ffree2SLT(float **A, int n);
double **darray2SLT(int n);
void dfree2SLT(double **A, int n);

#endif /* ndef ARRAY_H */
