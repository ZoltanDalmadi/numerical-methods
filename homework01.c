#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
  size_t _rows;
  size_t _cols;
  double **_data;
} MATRIX;

typedef struct {
  size_t _items;
  double *_data;
} VECTOR;

MATRIX *createMatrix(size_t rows, size_t cols) {
  size_t i;
  MATRIX *m = (MATRIX *)malloc(sizeof(MATRIX));
  m->_rows = rows;
  m->_cols = cols;
  m->_data = (double **)malloc(rows * sizeof(double *));

  for (i = 0; i < rows; ++i)
    m->_data[i] = (double *)calloc(cols, sizeof(double));

  return m;
}

void destroyMatrix(MATRIX *m) {
  size_t i;

  for (i = 0; i < m->_rows; ++i)
    free(m->_data[i]);

  free(m->_data);
  free(m);
  m = NULL;
}

VECTOR *createVector(size_t items) {
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

void swapMatrixRows(MATRIX *m, size_t a, size_t b) {
  double *temp = m->_data[a];
  m->_data[a] = m->_data[b];
  m->_data[b] = temp;
}

void swapVectorItems(VECTOR *v, size_t a, size_t b) {
  double temp = v->_data[a];
  v->_data[a] = v->_data[b];
  v->_data[b] = temp;
}

void printMatrix(MATRIX *m) {
  size_t i, j;

  for (i = 0; i < m->_rows; ++i) {
    printf("|");

    for (j = 0; j < m->_cols; ++j) {
      printf("%lf", m->_data[i][j]);

      if (j == m->_cols - 1)
        printf("|\n");
      else
        printf(" ");
    }
  }

  printf("\n");
}

void printVector(VECTOR *v) {
  size_t i;

  for (i = 0; i < v->_items; ++i)
    printf("|%lf|\n", v->_data[i]);

  printf("\n");
}

int main() {
  size_t n, i, j, k;
  scanf("%zu", &n);

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
  /* int singular = 0; */

  /* PLU FELBONTAS */
  for (k = 0; k < n - 1; ++k) {

    /* legnagyobb abszolutertek kivalasztasa */
    double max = A->_data[k][k];
    size_t maxRow = k;

    for (i = k + 1; i < n; ++i) {
      if (fabs(A->_data[i][k]) > fabs(max)) {
        max = A->_data[i][k];
        maxRow = i;
      }
    }

    /* sorcsere legnagyobb aszolutertekure */

    if (fabs(max) < 1e-15) {
      puts("szingularis 1");
      /* singular = 1; */
      break;
    } else if (maxRow != k) {
      swapMatrixRows(A, maxRow, k);
      swapVectorItems(P, maxRow, k);
      d *= -1;
    }

    /* kinullazas */
    for (i = k + 1; i < n; ++i) {
      /* divide whole row with its first item */
      A->_data[i][k] /= A->_data[k][k];

      for (j = k + 1; j < n; ++j) {
        A->_data[i][j] = A->_data[i][j] - (A->_data[i][k] * A->_data[k][j]);
      }
    }

    /* szingularitas ellenorzes */
    if (fabs(A->_data[n - 1][n - 1]) < 1e-15) {
      puts("szingularis 2");
      break;
    }
  }

  printVector(b);
  printVector(P);

  VECTOR *nb = createVector(n);

  /* b vektor permutalasa */
  for (i = 0; i < n; ++i)
    nb->_data[i] = b->_data[(size_t)P->_data[i]];

  /* swapVectorItems(b, i, (size_t)P->_data[i]); */

  printVector(nb);

  /* Ly = b megoldasa */
  VECTOR *y = createVector(n);

  double temp;

  for (i = 0; i < n; ++i) {
    temp = 0;

    for (j = 0; j < i; ++j) {
      /* y->_data[i] = b->_data[i] - A->_data[i][j] * y->_data[j]; */
      temp += A->_data[i][j] * y->_data[j];
    }

    y->_data[i] = nb->_data[i] - temp;
  }

  printVector(y);

  /* Ux = y megoldasa */
  VECTOR *x = createVector(n);

  for (i = n - 1; i >= 0; i--) {
    temp = 0;

    for (j = i + 1; j < n; ++j) {
      temp += A->_data[i][j] * x->_data[j];
    }

    x->_data[i] = (y->_data[i] - temp) / A->_data[i][i];
  }

  printVector(x);

  destroyVector(x);
  destroyVector(y);
  destroyVector(P);
  destroyVector(b);
  destroyMatrix(A);

  return 0;
}
