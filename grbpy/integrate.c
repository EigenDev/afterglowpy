#include <math.h>
#include "integrate.h"

#define KMAX 10

double trap(double (*f)(double, void *), double xa, double xb, int N, void *args)
{
    double dx = (xb - xa)/N;
    double I = 0.5*(f(xa, args) + f(xb, args));
    int i;
    for(i=1; i<N; i++)
        I += f(xa + i*dx, args);
    return I*dx;
}

double simp(double (*f)(double, void *), double xa, double xb, int N, void *args)
{
    if(N%2 == 1)
        N -= 1;

    double dx = (xb - xa)/N;
    double I1, I2, I3;
    int i;
    I1 = f(xa, args) + f(xb, args);
    I2 = 0.0;
    for(i=1; i<N; i+=2)
        I2 += f(xa + i*dx, args);
    I3 = 0.0;
    for(i=2; i<N; i+=2)
        I3 += f(xa + i*dx, args);
    return (I1 + 4*I2 + 2*I3) * dx / 3.0;
}

double romb(double (*f)(double, void *), double xa, double xb, int N, double atol, double rtol, void *args)
{
    double R[KMAX];

    int m, k, k0, fpm, Nk;
    double hk, Rp, err;

    hk = xb - xa;
    Nk = 1;
    R[KMAX-1] = 0.5*(xb-xa)*(f(xa, args) + f(xb, args));

    for(k=1; k<KMAX; k++)
    {
        k0 = KMAX-k-1;
        hk *= 0.5;
        Nk *= 2;

        Rp = 0.0;
        for(m=1; m<Nk; m+=2)
            Rp += f(xa + m*hk, args);
        R[k0] = 0.5*R[k0+1] + hk*Rp;

        fpm = 1;
        for(m=1; m<=k; m++)
        {
            fpm *= 4;
            R[k0+m] = (fpm*R[k0+m-1] - R[k0+m]) / (fpm - 1);
        }
        err = (R[KMAX-1] - R[0]) / (fpm - 1);
        R[0] = R[KMAX-1];

        if(fabs(err) < atol + rtol*fabs(R[0]))
            break;

        if(N > 1 && Nk >= N)
            break;
    }

    return R[0];
}
