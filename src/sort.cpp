//
// Taken from functions in R/src/main/sort.c with small modifications
// (it does not seem that these functions are included in the standalone R library).
//
void rsort_with_index(double *x, int *indx, unsigned n)
{
    double v;
    int i, j, h, iv;

    for (h = 1; h <= n / 9; h = 3 * h + 1);
    for (; h > 0; h /= 3)
        for (i = h; i < n; i++) {
            v = x[i]; iv = indx[i];
            j = i;
            while (j >= h && x[j - h] > v)
            { x[j] = x[j - h]; indx[j] = indx[j-h]; j -= h; }
            x[j] = v; indx[j] = iv;
        }
}

void revsort(double *a, int *ib, int n)
{
/* Sort a[] into descending order by "heapsort";
 * sort ib[] alongside;
 * if initially, ib[] = 1...n, it will contain the permutation finally
 */

    int l, j, ir, i;
    double ra;
    int ii;

    if (n <= 1) return;

    a--; ib--;

    l = (n >> 1) + 1;
    ir = n;

    for (;;) {
        if (l > 1) {
            l = l - 1;
            ra = a[l];
            ii = ib[l];
        }
        else {
            ra = a[ir];
            ii = ib[ir];
            a[ir] = a[1];
            ib[ir] = ib[1];
            if (--ir == 1) {
                a[1] = ra;
                ib[1] = ii;
                return;
            }
        }
        i = l;
        j = l << 1;
        while (j <= ir) {
            if (j < ir && a[j] > a[j + 1]) ++j;
            if (ra > a[j]) {
                a[i] = a[j];
                ib[i] = ib[j];
                j += (i = j);
            }
            else
                j = ir + 1;
        }
        a[i] = ra;
        ib[i] = ii;
    }
}
