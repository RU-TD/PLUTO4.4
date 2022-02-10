#include <stdio.h>

void lkrm_c(double *, double *, double *);
void lkrm_init_c();

int main(){
  double n[2] = {0e0};
  double Tgas = 1e2;
  double dt = 1e9;

  lkrm_init_c();
  n[0] = 1e8;

  printf("%g %g %g %g\n", n[0], n[1], Tgas, dt);
  lkrm_c(n, &Tgas, &dt);
  printf("%g %g %g %g\n", n[0], n[1], Tgas, dt);

  return 0;
}
