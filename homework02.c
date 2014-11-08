#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct Matrix
{
  int _rows;
  int _cols;
  double **_data;
} MATRIX;

typedef struct Vector
{
  int _items;
  double *_data;
} VECTOR;

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

void swapMatrixRows(MATRIX *m, int a, int b)
{
  double *temp = m->_data[a];
  m->_data[a] = m->_data[b];
  m->_data[b] = temp;
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

  for (i = 0; i < result->_items; ++i)
  {
    result->_data[i] = 0;
  }

  for (i = 0; i < A->_rows; ++i)
  {
    /* result->_data[i] = 0; */

    for (j = 0; j < A->_cols; ++j)
      result->_data[i] += A->_data[i][j] * B->_data[j];
  }

  return result;
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

int main()
{
  /* INPUT ----------------------------------------------------------------- */
  int f, n, N, m, i, j, k, l, maxit;
  double c, epsilon, lambdaOld, lambdaNew, temp;

  /* matrixok szama */
  scanf("%d", &N);

  /* Matrixok */
  for (l = 0; l < N; ++l)
  {
    /* Matrix merete */
    scanf("%d", &n);

    /* A matrix */
    MATRIX *A = createMatrix(n, n);

    for (i = 0; i < n; ++i)
      for (j = 0; j < n; ++j)
        scanf("%lf", &A->_data[i][j]);

    /* matrixra vonatkozo feladatok szama */
    scanf("%d", &m);

    /* feladatok elkezdese ------------------------------------------------- */
    for (f = 0; f < m; ++f)
    {
      /* eltolasi parameter */
      scanf("%lf", &c);

      /* maximalis iteracio */
      scanf("%d", &maxit);

      /* leallasi feltetel epszilonja */
      scanf("%lf", &epsilon);

      /* y vektor */
      VECTOR *y = createVector(n);

      for (i = 0; i < n; ++i)
        scanf("%lf", &y->_data[i]);

      /* P vektor */
      VECTOR *P = createVector(n);

      for (i = 0; i < n; ++i)
        P->_data[i] = i;

      /* szingularis-e? */
      int singular = 0;

      /* Masolat keszitese A-rol */
      MATRIX *Aplu = copyMatrix(A);

      /* A-cE */
      for (i = 0; i < n; ++i)
        Aplu->_data[i][i] -= c;

      /* 1. A-cE PLU FELBONTASA ---------------------------------------------- */
      for (k = 0; k < n - 1; ++k)
      {
        /* legnagyobb abszolutertek kivalasztasa */
        double max = Aplu->_data[k][k];
        int maxRow = k;

        for (i = k + 1; i < n; ++i)
        {
          if (fabs(Aplu->_data[i][k]) > fabs(max))
          {
            max = Aplu->_data[i][k];
            maxRow = i;
          }
        }

        /* sorcsere legnagyobb aszolutertekure */
        if (fabs(max) < 1e-15)
        {
          singular = 1;
          break;
        }
        else if (maxRow != k)
        {
          swapMatrixRows(Aplu, maxRow, k);
          swapVectorItems(P, maxRow, k);
        }

        /* kinullazas */
        for (i = k + 1; i < n; ++i)
        {
          Aplu->_data[i][k] /= Aplu->_data[k][k];

          for (j = k + 1; j < n; ++j)
          {
            Aplu->_data[i][j] =
              Aplu->_data[i][j] - (Aplu->_data[i][k] * Aplu->_data[k][j]);
          }
        }
      }

      /* A-cE PLU FELBONTAS VEGE --------------------------------------------- */

      /* szingularitas ellenorzes */
      if (singular)
      {
        printf("%.8lf \n", c);
        destroyVector(y);
        destroyVector(P);
        destroyMatrix(Aplu);
        continue;
      }

      if (fabs(Aplu->_data[n - 1][n - 1]) < 1e-15)
      {
        printf("%.8lf \n", c);
        destroyVector(y);
        destroyVector(P);
        destroyMatrix(Aplu);
        continue;
      }

      /* 2. y0 norma kiszamitasa ------------------------------------------- */
      double ynorma = sqrt(innerProduct(y, y));

      if (fabs(ynorma) < 1e-15)
      {
        puts("kezdovektor");
        destroyVector(y);
        destroyVector(P);
        destroyMatrix(Aplu);
        continue;
      }

      /* 3. y normalas es elso lambda (Rayleigh-hanyados) ------------------ */
      for (i = 0; i < n; ++i)
        y->_data[i] /= ynorma;

      VECTOR *Ay = multiplyMatrixVector(A, y);

      lambdaOld = innerProduct(Ay, y);

      /* 4. INVERZ ITERACIO ------------------------------------------------ */
      for (k = 0; k < maxit; ++k)
      {
        VECTOR *yp = createVector(n);

        /* y vektor permutalasa */
        for (i = 0; i < n; ++i)
        {
          yp->_data[i] = y->_data[(int)P->_data[i]];
        }

        /* Lz = yp megoldasa (z y-ban tarolva) */
        for (i = 0; i < n; ++i)
        {
          temp = 0;

          for (j = 0; j < i; ++j)
          {
            temp += Aplu->_data[i][j] * y->_data[j];
          }

          y->_data[i] = yp->_data[i] - temp;
        }

        /* Ux = z megoldasa (z y-ban tarolva, x yp-ben) */
        for (i = n - 1; i >= 0; i--)
        {
          temp = 0;

          for (j = i + 1; j < n; j++)
          {
            temp += Aplu->_data[i][j] * yp->_data[j];
          }

          yp->_data[i] = (y->_data[i] - temp) / Aplu->_data[i][i];
        }

        /* kapott x (yp-ben tarolva) vektort normaljuk */
        double xnorma = sqrt(innerProduct(yp, yp));

        for (i = 0; i < n; ++i)
          y->_data[i] = yp->_data[i] / xnorma;

        /* sajatertek kozelitese */
        destroyVector(Ay);
        Ay = multiplyMatrixVector(A, y);

        lambdaNew = innerProduct(Ay, y);

        /* leallasi feltetel */
        if (fabs(lambdaNew - lambdaOld) <= epsilon * (1 + fabs(lambdaNew)))
        {
          destroyVector(yp);
          break;
        }

        lambdaOld = lambdaNew;
        destroyVector(yp);

      } /* INVERZ ITERACIO VEGE */

      /* 5. Ha elertuk maxit-et, hibaval kilepes --------------------------- */
      if (k == maxit)
      {
        puts("maxit");
        destroyVector(Ay);
        destroyVector(y);
        destroyVector(P);
        destroyMatrix(Aplu);
        continue;
      }

      /* 6. Sikerteszt ----------------------------------------------------- */
      VECTOR *lambday = createVector(n);

      /* lambda * y */
      for (i = 0; i < n; ++i)
        lambday->_data[i] = lambdaNew * y->_data[i];

      /* Ay - lambda * y (Ay-ben tarolva) */
      for (i = 0; i < n; ++i)
        Ay->_data[i] -= lambday->_data[i];

      double amount = innerProduct(Ay, Ay);

      /* sikerteszt elvegzese */
      if (amount <= epsilon)
        printf("siker ");
      else
        printf("sikertelen ");

      printf("%.8lf ", lambdaNew);
      printVector(y);
      printf("%.8lf %d\n", amount, k);

      /* takaritas */
      destroyVector(lambday);
      destroyVector(Ay);
      destroyVector(y);
      destroyVector(P);
      destroyMatrix(Aplu);

    } /* feladat vege */

    /* takaritas */
    destroyMatrix(A);

  } /* matrix vege */

  return 0;
}
