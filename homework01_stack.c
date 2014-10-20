#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void swapVectorItems(double *v, int a, int b) {
  double temp = v[a];
  v[a] = v[b];
  v[b] = temp;
}

void printVector(double *v, int n) {
  int i;

  for (i = 0; i < n; ++i)
    printf("%.8lf ", v[i]);

  printf("\n");
}

int main() {
  int n, m, i, j, k, l;

  scanf("%d", &m);

  for (l = 0; l < m; ++l) {
    scanf("%d", &n);

    /* A matrix */
    double A[n][n];

    for (i = 0; i < n; ++i)
      for (j = 0; j < n; ++j)
        scanf("%lf", &A[i][j]);

    /* b vektor */
    double b[n];

    for (i = 0; i < n; ++i)
      scanf("%lf", &b[i]);

    /* P vektor */
    double P[n];

    for (i = 0; i < n; ++i)
      P[i] = i;

    /* determinans elojel */
    int d = 1;

    /* szingularis-e? */
    int singular = 0;

    /* PLU FELBONTAS */
    for (k = 0; k < n - 1; ++k) {

      /* legnagyobb abszolutertek kivalasztasa */
      double max = A[k][k];
      int maxRow = k;

      for (i = k + 1; i < n; ++i) {
        if (fabs(A[i][k]) > fabs(max)) {
          max = A[i][k];
          maxRow = i;
        }
      }

      if (fabs(max) < 1e-15) {
        printf("szingularis\n");
        singular = 1;
        break;
      }

      /* sorcsere legnagyobb aszolutertekure */
      if (maxRow != k) {
        for (i = 0; i < n; ++i) {
          double temp = A[k][i];
          A[k][i] = A[maxRow][i];
          A[maxRow][i] = temp;
        }

        swapVectorItems(P, maxRow, k);
        d *= -1;
      }

      /* kinullazas */
      for (i = k + 1; i < n; ++i) {
        A[i][k] /= A[k][k];

        for (j = k + 1; j < n; ++j) {
          A[i][j] = A[i][j] - A[i][k] * A[k][j];
        }
      }
    }

    if (singular) {
      continue;
    }

    /* szingularitas ellenorzes */
    if (fabs(A[n - 1][n - 1]) < 1e-15) {
      printf("szingularis\n");
      continue;
    }

    /* b vektor permutalasa */
    double nb[n];

    for (i = 0; i < n; ++i)
      nb[i] = b[(int)P[i]];

    /* Ly = b megoldasa */
    double temp;

    for (i = 0; i < n; ++i) {
      temp = 0;

      for (j = 0; j < i; ++j) {
        temp += A[i][j] * b[j];
      }

      b[i] = nb[i] - temp;
    }

    /* Ux = y megoldasa */
    for (i = n - 1; i >= 0; i--) {
      temp = 0;

      for (j = i + 1; j < n; j++) {
        temp += A[i][j] * nb[j];
      }

      nb[i] = (b[i] - temp) / A[i][i];
    }

    /* determinans */
    double det = 1;

    for (i = 0; i < n; ++i) {
      det *= A[i][i];
    }

    det *= d;

    printf("%.8lf ", det);
    printVector(nb, n);
  }

  return 0;
}
