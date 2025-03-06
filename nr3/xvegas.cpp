#include "nr3.h"
#include "ran.h"
#include "rebin.h"
#include "vegas.h"

Doub torusfunc(const VecDoub &x, const Doub wgt) {
  Doub den = exp(5. * x[2]);
  if (SQR(x[2]) + SQR(sqrt(SQR(x[0]) + SQR(x[1])) - 3.) <= 1.)
    return den;
  else
    return 0.;
}

int main() {
  Doub tgral, sd, chi2a;
  VecDoub regn(6);
  regn[0] = 1.;   regn[3] = 4.;
  regn[1] = -3.;  regn[4] = 4.;
  regn[2] = -1.;  regn[5] = 1.;
  vegas(regn, torusfunc, 0, 10000, 10, 0, tgral, sd, chi2a);
  vegas(regn, torusfunc, 1, 900000, 1, 0, tgral, sd, chi2a);
  return 0;
}