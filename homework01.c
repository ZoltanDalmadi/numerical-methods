#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
  int _rows;
  int _cols;
  double **_data;
} MATRIX;

typedef struct {
  int _items;
  double *_data;
} VECTOR;

MATRIX *createMatrix(int rows, int cols) {
  int i;
  MATRIX *m = (MATRIX *)malloc(sizeof(MATRIX));
  m->_rows = rows;
  m->_cols = cols;
  m->_data = (double **)malloc(rows * sizeof(double *));

  for (i = 0; i < rows; ++i)
    m->_data[i] = (double *)malloc(cols * sizeof(double));

  return m;
}

void destroyMatrix(MATRIX *m) {
  int i;

  for (i = 0; i < m->_rows; ++i)
    free(m->_data[i]);

  free(m->_data);
  free(m);
  m = NULL;
}

VECTOR *createVector(int items) {
  VECTOR *v = (VECTOR *)malloc(sizeof(VECTOR));
  v->_items = items;
  v->_data = (double *)malloc(items * sizeof(double));
  return v;
}

void destroyVector(VECTOR *v) {
  free(v->_data);
  free(v);
  v = NULL;
}

void swapMatrixRows(MATRIX *m, int a, int b) {
  double *temp = m->_data[a];
  m->_data[a] = m->_data[b];
  m->_data[b] = temp;
}

void swapVectorItems(VECTOR *v, int a, int b) {
  double temp = v->_data[a];
  v->_data[a] = v->_data[b];
  v->_data[b] = temp;
}

void printVector(VECTOR *v) {
  int i;

  for (i = 0; i < v->_items; ++i)
    printf("%.8lf ", v->_data[i]);

  printf("\n");
}

int main() {
  int n, m, i, j, k, l;
  scanf("%d", &m);

  for (l = 0; l < m; ++l) {
    scanf("%d", &n);

    /* A matrix */
    MATRIX *A = createMatrix(n, n);

    for (i = 0; i < n; ++i)
      for (j = 0; j < n; ++j)
        scanf("%lf", &A->_data[i][j]);

    /* b vektor */
    VECTOR *b = createVector(n);

    for (i = 0; i < n; ++i)
      scanf("%lf", &b->_data[i]);

    /* P vektor */
    VECTOR *P = createVector(n);

    for (i = 0; i < n; ++i)
      P->_data[i] = i;

    /* determinans elojel */
    int d = 1;

    /* szingularis-e? */
    int singular = 0;

    /* PLU FELBONTAS */
    for (k = 0; k < n - 1; ++k) {

      /* legnagyobb abszolutertek kivalasztasa */
      double max = A->_data[k][k];
      int maxRow = k;

      for (i = k + 1; i < n; ++i) {
        if (fabs(A->_data[i][k]) > fabs(max)) {
          max = A->_data[i][k];
          maxRow = i;
        }
      }

      /* sorcsere legnagyobb aszolutertekure */

      if (fabs(max) < 1e-15) {
        printf("szingularis\n");
        singular = 1;
        break;
      } else if (maxRow != k) {
        swapMatrixRows(A, maxRow, k);
        swapVectorItems(P, maxRow, k);
        d *= -1;
      }

      /* kinullazas */
      for (i = k + 1; i < n; ++i) {
        A->_data[i][k] /= A->_data[k][k];

        for (j = k + 1; j < n; ++j) {
          A->_data[i][j] = A->_data[i][j] - (A->_data[i][k] * A->_data[k][j]);
        }
      }
    }

    /* szingularitas ellenorzes */
    if (fabs(A->_data[n - 1][n - 1]) < 1e-15) {
      printf("szingularis\n");
      singular = 1;
    }

    if (singular) {
      destroyVector(P);
      destroyVector(b);
      destroyMatrix(A);
      continue;
    }

    VECTOR *nb = createVector(n);

    /* b vektor permutalasa */
    for (i = 0; i < n; ++i)
      nb->_data[i] = b->_data[(int)P->_data[i]];

    /* Ly = b megoldasa */

    double temp;

    for (i = 0; i < n; ++i) {
      temp = 0;

      for (j = 0; j < i; ++j) {
        temp += A->_data[i][j] * b->_data[j];
      }

      b->_data[i] = nb->_data[i] - temp;
    }

    /* Ux = y megoldasa */

    for (i = n - 1; i >= 0; i--) {
      temp = 0;

      for (j = i + 1; j < n; j++) {
        temp += A->_data[i][j] * nb->_data[j];
      }

      nb->_data[i] = (b->_data[i] - temp) / A->_data[i][i];
    }

    /* determinans */
    double det = 1;

    for (i = 0; i < n; ++i) {
      det *= A->_data[i][i];
    }

    det *= d;

    printf("%.8lf ", det);
    printVector(nb);

    destroyVector(nb);
    destroyVector(P);
    destroyVector(b);
    destroyMatrix(A);
  }

  return 0;
}
