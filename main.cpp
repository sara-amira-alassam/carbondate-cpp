#include <iostream>
#define MATHLIB_STANDALONE
#include <Rmath.h>

int main() {
    double ans = dnorm4(1.0, 2.0, 0.5, 0);
    printf("ans %f\n", ans);
    std::cout << "Hello, World!" << std::endl;
    return 0;
}

