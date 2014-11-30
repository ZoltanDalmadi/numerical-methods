#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* ----------------------------------------------------------------------------
 * MATRIX struct and functions
 * ------------------------------------------------------------------------- */
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

void printMatrix(MATRIX *m)
{
  int i, j;

  for (i = 0; i < m->_rows; ++i)
  {
    printf("|");

    for (j = 0; j < m->_cols; ++j)
    {
      printf("%lf", m->_data[i][j]);

      if (j == m->_cols - 1)
        printf("|\n");
      else
        printf(" ");
    }
  }

  printf("\n");
}

/* ----------------------------------------------------------------------------
 * VECTOR struct and functions
 * ------------------------------------------------------------------------- */
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
 * ------------------------------------------------------------------------- */
typedef double(*equation_ptr)(VECTOR *);
typedef MATRIX *(*jacobifunc_ptr)(VECTOR *);

typedef struct System
{
  int _items;
  equation_ptr *_data;
  jacobifunc_ptr _jacobi;
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
  s->_data = NULL;
  free(s);
  s = NULL;
}

/* ----------------------------------------------------------------------------
 * Utility functions
 * ------------------------------------------------------------------------- */
double min(double a, double b)
{
  return a < b ? a : b;
}

/* ----------------------------------------------------------------------------
 * Equation system #1
 * ------------------------------------------------------------------------- */
double sys1_func1(VECTOR *v)
{
  double x1 = v->_data[0];
  double x3 = v->_data[2];
  return -(x1 * x1) + x3 + 3;
}

double sys1_func2(VECTOR *v)
{
  double x1 = v->_data[0];
  double x2 = v->_data[1];
  double x3 = v->_data[2];
  return -x1 + (2 * x2 * x2) - (x3 * x3) - 3;
}

double sys1_func3(VECTOR *v)
{
  double x2 = v->_data[1];
  double x3 = v->_data[2];
  return x2 - (3 * x3 * x3) + 2;
}

MATRIX *sys1_jacobi(VECTOR *v)
{
  MATRIX *result = createMatrix(v->_items, v->_items);

  double x1 = v->_data[0];
  double x2 = v->_data[1];
  double x3 = v->_data[2];

  result->_data[0][0] = -2 * x1;
  result->_data[0][1] = 0;
  result->_data[0][2] = 1;
  result->_data[1][0] = -1;
  result->_data[1][1] = 4 * x2;
  result->_data[1][2] = -2 * x3;
  result->_data[2][0] = 0;
  result->_data[2][1] = 1;
  result->_data[2][2] = -6 * x3;

  return result;
}

/* ----------------------------------------------------------------------------
 * Equation system #2
 * ------------------------------------------------------------------------- */
double sys2_func1(VECTOR *v)
{
  double x1 = v->_data[0];
  double x2 = v->_data[1];
  return (2 * x1 * x1) - x2 - 1;
}

double sys2_func2(VECTOR *v)
{
  double x1 = v->_data[0];
  double x2 = v->_data[1];
  return -x1 + (2 * x2 * x2) - 1;
}

MATRIX *sys2_jacobi(VECTOR *v)
{
  MATRIX *result = createMatrix(v->_items, v->_items);

  double x1 = v->_data[0];
  double x2 = v->_data[1];

  result->_data[0][0] = 4 * x1;
  result->_data[0][1] = -1;
  result->_data[1][0] = -1;
  result->_data[1][1] = 4 * x2;

  return result;
}

/* ----------------------------------------------------------------------------
 * Equation system #3
 * ------------------------------------------------------------------------- */
double sys3_func1(VECTOR *v)
{
  double x1 = v->_data[0];
  double x2 = v->_data[1];
  return -4 * x1 + cos(2 * x1 - x2) - 3;
}

double sys3_func2(VECTOR *v)
{
  double x1 = v->_data[0];
  double x2 = v->_data[1];
  return sin(x1) - 3 * x2 - 2;
}

MATRIX *sys3_jacobi(VECTOR *v)
{
  MATRIX *result = createMatrix(v->_items, v->_items);

  double x1 = v->_data[0];
  double x2 = v->_data[1];

  result->_data[0][0] = -4 - sin(2 * x1 - x2) * 2;
  result->_data[0][1] = sin(2 * x1 - x2);
  result->_data[1][0] = cos(x1);
  result->_data[1][1] = -3;

  return result;
}

/* ----------------------------------------------------------------------------
 * Equation system #4
 * ------------------------------------------------------------------------- */
double sys4_func1(VECTOR *v)
{
  double x1 = v->_data[0];
  double x2 = v->_data[1];
  return (x1 * x2 * x2) - (4 * x1 * x2) + (4 * x1) - 1;
}

double sys4_func2(VECTOR *v)
{
  double x1 = v->_data[0];
  double x2 = v->_data[1];
  return exp(x1 - 1) - x2 + 1;
}

MATRIX *sys4_jacobi(VECTOR *v)
{
  MATRIX *result = createMatrix(v->_items, v->_items);

  double x1 = v->_data[0];
  double x2 = v->_data[1];

  result->_data[0][0] = (x2 * x2) - 4 * (x2 + 1);
  result->_data[0][1] = 2 * x1 * (x2 - 2);
  result->_data[1][0] = exp(x1 - 1);
  result->_data[1][1] = -1;

  return result;
}

/* ----------------------------------------------------------------------------
 * Equation systems miscellaneous
 * ------------------------------------------------------------------------- */
VECTOR *calcSystem(SYSTEM *s, VECTOR *v)
{
  VECTOR *result = createVector(v->_items);

  int i;

  for (i = 0; i < result->_items; ++i)
    result->_data[i] = s->_data[i](v);

  return result;
}

/* ----------------------------------------------------------------------------
 * Main function
 * ------------------------------------------------------------------------- */
int main()
{
  /* variables */
  int i, n, N, k, task, maxit;
  double epsilon, t;
  SYSTEM *activeSystem;

  /* construct equation systems */
  SYSTEM *system1 = createSystem(3);
  system1->_data[0] = sys1_func1;
  system1->_data[1] = sys1_func2;
  system1->_data[2] = sys1_func3;
  system1->_jacobi = sys1_jacobi;

  SYSTEM *system2 = createSystem(2);
  system2->_data[0] = sys2_func1;
  system2->_data[1] = sys2_func2;
  system2->_jacobi = sys2_jacobi;

  SYSTEM *system3 = createSystem(2);
  system3->_data[0] = sys3_func1;
  system3->_data[1] = sys3_func2;
  system3->_jacobi = sys3_jacobi;

  SYSTEM *system4 = createSystem(2);
  system4->_data[0] = sys4_func1;
  system4->_data[1] = sys4_func2;
  system4->_jacobi = sys4_jacobi;

  /* input number of tasks to solve */
  scanf("%d", &N);

  /* task loop */
  for (task = 0; task < N; ++task)
  {
    /* determine actual equation system */
    scanf("%d", &n);

    switch (n)
    {
    case 1:
      activeSystem = system1;
      break;

    case 2:
      activeSystem = system2;
      break;

    case 3:
      activeSystem = system3;
      break;

    case 4:
      activeSystem = system4;
      break;
    }

    /* input maxit and epsilon */
    scanf("%d", &maxit);
    scanf("%lf", &epsilon);

    /* input starting vector x0 */
    VECTOR *xk = createVector(activeSystem->_items);

    for (i = 0; i < xk->_items; ++i)
      scanf("%lf", &xk->_data[i]);

    /* start k loop */
    for (k = 1; k < maxit; ++k)
    {
      /* 1. calculate jacobi matrix and f vector --------------------------- */
      MATRIX *jacobixk = activeSystem->_jacobi(xk);
      VECTOR *fxk = calcSystem(activeSystem, xk);

      /* 2. solve jacobi * deltax = -f linear equation (PLU) --------------- */



      destroyVector(fxk);
      destroyMatrix(jacobixk);

    } /* end k loop */

    /* cleanup */
    destroyVector(xk);

  } /* end task loop */

  /* cleanup */
  destroySystem(system4);
  destroySystem(system3);
  destroySystem(system2);
  destroySystem(system1);

  return 0;
}
