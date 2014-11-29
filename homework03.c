#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* ----------------------------------------------------------------------------
 * MATRIX struct and functions
 * --------------------------------------------------------------------------*/
typedef struct Matrix
{
  int _rows;
  int _cols;
  double **_data;
} MATRIX;

MATRIX *createMatrix(int rows, int cols)
{
  int i;
  MATRIX *m = (MATRIX *)malloc(sizeof(MATRIX));
  m->_rows = rows;
  m->_cols = cols;
  m->_data = (double **)malloc(rows * sizeof(double *));

  for (i = 0; i < rows; ++i)
    m->_data[i] = (double *)malloc(cols * sizeof(double));

  return m;
}

void destroyMatrix(MATRIX *m)
{
  int i;

  for (i = 0; i < m->_rows; ++i)
    free(m->_data[i]);

  free(m->_data);
  free(m);
  m = NULL;
}

void swapMatrixRows(MATRIX *m, int a, int b)
{
  double *temp = m->_data[a];
  m->_data[a] = m->_data[b];
  m->_data[b] = temp;
}

MATRIX *copyMatrix(MATRIX *m)
{
  MATRIX *copy = createMatrix(m->_rows, m->_cols);

  int i, j;

  for (i = 0; i < copy->_rows; ++i)
  {
    for (j = 0; j < copy->_cols; ++j)
      copy->_data[i][j] = m->_data[i][j];
  }

  return copy;
}

/* ----------------------------------------------------------------------------
 * VECTOR struct and functions
 * --------------------------------------------------------------------------*/
typedef struct Vector
{
  int _items;
  double *_data;
} VECTOR;

VECTOR *createVector(int items)
{
  VECTOR *v = (VECTOR *)malloc(sizeof(VECTOR));
  v->_items = items;
  v->_data = (double *)malloc(items * sizeof(double));
  return v;
}

void destroyVector(VECTOR *v)
{
  free(v->_data);
  free(v);
  v = NULL;
}

void swapVectorItems(VECTOR *v, int a, int b)
{
  double temp = v->_data[a];
  v->_data[a] = v->_data[b];
  v->_data[b] = temp;
}

void printVector(VECTOR *v)
{
  int i;

  for (i = 0; i < v->_items; ++i)
    printf("%.8lf ", v->_data[i]);
}

double innerProduct(VECTOR *a, VECTOR *b)
{
  double result = 0;
  int i;

  for (i = 0; i < a->_items; ++i)
  {
    result += a->_data[i] * b->_data[i];
  }

  return result;
}

VECTOR *multiplyMatrixVector(MATRIX *A, VECTOR *B)
{
  VECTOR *result = createVector(A->_rows);

  int i, j;

  for (i = 0; i < A->_rows; ++i)
  {
    result->_data[i] = 0;

    for (j = 0; j < A->_cols; ++j)
      result->_data[i] += A->_data[i][j] * B->_data[j];
  }

  return result;
}

double infiniteNorm(VECTOR *v)
{
  double infnorm = fabs(v->_data[0]);
  int i;

  for (i = 1; i < v->_items; ++i)
  {
    double value = fabs(v->_data[i]);
    infnorm = value > infnorm ? value : infnorm;
  }

  return infnorm;
}

/* ----------------------------------------------------------------------------
 * SYSTEM struct and functions
 * --------------------------------------------------------------------------*/
typedef double(*equation_ptr)(VECTOR *);

typedef struct System
{
  int _items;
  equation_ptr *_data;
} SYSTEM;

SYSTEM *createSystem(int items)
{
  SYSTEM *s = (SYSTEM *)malloc(sizeof(SYSTEM));
  s->_items = items;
  s->_data = (equation_ptr *)malloc(items * sizeof(equation_ptr));
  return s;
}

void destroySystem(SYSTEM *s)
{
  free(s->_data);
  free(s);
  s = NULL;
}

/* ----------------------------------------------------------------------------
 * Utility functions
 * --------------------------------------------------------------------------*/
double min(double a, double b)
{
  return a < b ? a : b;
}

double func1(VECTOR *v)
{
  double x1 = v->_data[0];
  double x3 = v->_data[2];
  return -(x1 * x1) + x3 + 3;
}

double func2(VECTOR *v)
{
  double x1 = v->_data[0];
  double x2 = v->_data[1];
  double x3 = v->_data[2];
  return -x1 + (2 * x2 * x2) - (x3 * x3) - 3;
}

double func3(VECTOR *v)
{
  double x2 = v->_data[1];
  double x3 = v->_data[2];
  return x2 - (3 * x3 * x3) + 2;
}

VECTOR *calcSystem(SYSTEM *s, VECTOR *v)
{
  VECTOR *result = createVector(v->_items);

  int i;

  for (i = 0; i < result->_items; ++i)
    result->_data[i] = s->_data[i](v);

  return result;
}

int main()
{
  int i;
  SYSTEM *system1 = createSystem(3);
  system1->_data[0] = func1;
  system1->_data[1] = func2;
  system1->_data[2] = func3;

  VECTOR *vec = createVector(3);

  for (i = 0; i < 3; ++i)
    scanf("%lf", &vec->_data[i]);

  VECTOR *res = calcSystem(system1, vec);

  printVector(res);
  printf("\n");

  destroyVector(res);
  destroyVector(vec);
  destroySystem(system1);

  return 0;
}
