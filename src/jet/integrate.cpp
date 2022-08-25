#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include "include/integrate.hpp"

using namespace afterglowpy::mesh;
constexpr int kmax = 20;
constexpr int max_gauss_n = 10;

static const double g4x[4] = {-0.86113631159405257522,
                              -0.33998104358485626480,
                              0.33998104358485626480,
                              0.86113631159405257522};
static const double g4w[4] = {0.34785484513745385737,
                              0.65214515486254614262,
                              0.65214515486254614262,
                              0.34785484513745385737};
static const double k9x[5] = {-0.97656025073757311153,
                              -0.64028621749630998240,
                              0,
                              0.64028621749630998240,
                              0.97656025073757311153};
static const double k9w[9] = {0.06297737366547301477,
                              0.17005360533572272680,
                              0.26679834045228444803,
                              0.32694918960145162956,
                              0.34644298189013636168,
                              0.32694918960145162956,
                              0.26679834045228444803,
                              0.17005360533572272680,
                              0.06297737366547301477};

static const double g7x[7] = {-9.491079123427585245261896840478513e-01,
                              -7.415311855993944398638647732807884e-01,
                              -4.058451513773971669066064120769615e-01,
                              0.0,
                              4.058451513773971669066064120769615e-01,
                              7.415311855993944398638647732807884e-01,
                              9.491079123427585245261896840478513e-01};
static const double g7w[7] = {1.294849661688696932706114326790820e-01,
                              2.797053914892766679014677714237796e-01,
                              3.818300505051189449503697754889751e-01,
                              4.179591836734693877551020408163265e-01,
                              3.818300505051189449503697754889751e-01,
                              2.797053914892766679014677714237796e-01,
                              1.294849661688696932706114326790820e-01};
static const double k15x[8] = {-9.914553711208126392068546975263285e-01,
                               -8.648644233597690727897127886409262e-01,
                               -5.860872354676911302941448382587296e-01,
                               -2.077849550078984676006894037732449e-01,
                               2.077849550078984676006894037732449e-01,
                               5.860872354676911302941448382587296e-01,
                               8.648644233597690727897127886409262e-01,
                               9.914553711208126392068546975263285e-01};
static const double k15w[15] = {2.293532201052922496373200805896959e-02,
                                6.309209262997855329070066318920429e-02,
                                1.047900103222501838398763225415180e-01,
                                1.406532597155259187451895905102379e-01,
                                1.690047266392679028265834265985503e-01,
                                1.903505780647854099132564024210137e-01,
                                2.044329400752988924141619992346491e-01,
                                2.094821410847278280129991748917143e-01,
                                2.044329400752988924141619992346491e-01,
                                1.903505780647854099132564024210137e-01,
                                1.690047266392679028265834265985503e-01,
                                1.406532597155259187451895905102379e-01,
                                1.047900103222501838398763225415180e-01,
                                6.309209262997855329070066318920429e-02,
                                2.293532201052922496373200805896959e-02};

static const double g10x[10] = {-9.739065285171717200779640120844521e-01,
                                -8.650633666889845107320966884234930e-01,
                                -6.794095682990244062343273651148736e-01,
                                -4.333953941292471907992659431657842e-01,
                                -1.488743389816312108848260011297200e-01,
                                1.488743389816312108848260011297200e-01,
                                4.333953941292471907992659431657842e-01,
                                6.794095682990244062343273651148736e-01,
                                8.650633666889845107320966884234930e-01,
                                9.739065285171717200779640120844521e-01};
static const double g10w[10] = {6.667134430868813759356880989333179e-02,
                                1.494513491505805931457763396576973e-01,
                                2.190863625159820439955349342281632e-01,
                                2.692667193099963550912269215694694e-01,
                                2.955242247147528701738929946513383e-01,
                                2.955242247147528701738929946513383e-01,
                                2.692667193099963550912269215694694e-01,
                                2.190863625159820439955349342281632e-01,
                                1.494513491505805931457763396576973e-01,
                                6.667134430868813759356880989333179e-02};
static const double k21x[11] = {-9.956571630258080807355272806890028e-01,
                                -9.301574913557082260012071800595083e-01,
                                -7.808177265864168970637175783450424e-01,
                                -5.627571346686046833390000992726941e-01,
                                -2.943928627014601981311266031038656e-01,
                                0.0,
                                2.943928627014601981311266031038656e-01,
                                5.627571346686046833390000992726941e-01,
                                7.808177265864168970637175783450424e-01,
                                9.301574913557082260012071800595083e-01,
                                9.956571630258080807355272806890028e-01};
static const double k21w[21] = {1.169463886737187427806439606219205e-02,
                                3.255816230796472747881897245938976e-02,
                                5.475589657435199603138130024458018e-02,
                                7.503967481091995276704314091619001e-02,
                                9.312545458369760553506546508336634e-02,
                                1.093871588022976418992105903258050e-01,
                                1.234919762620658510779581098310742e-01,
                                1.347092173114733259280540017717068e-01,
                                1.427759385770600807970942731387171e-01,
                                1.477391049013384913748415159720680e-01,
                                1.494455540029169056649364683898212e-01,
                                1.477391049013384913748415159720680e-01,
                                1.427759385770600807970942731387171e-01,
                                1.347092173114733259280540017717068e-01,
                                1.234919762620658510779581098310742e-01,
                                1.093871588022976418992105903258050e-01,
                                9.312545458369760553506546508336634e-02,
                                7.503967481091995276704314091619001e-02,
                                5.475589657435199603138130024458018e-02,
                                3.255816230796472747881897245938976e-02,
                                1.169463886737187427806439606219205e-02};


namespace afterglowpy
{
    namespace integrate
    {
        double trap(
            std::function<double(double, void *)> func, 
            double xa, 
            double xb, 
            int N,
            void *args, 
            std::function<int(void *)> errf
        )
        {
            double dx = (xb - xa)/N;
            double I  = 0.5 * (func(xa, args) + func(xb, args));
            for(int i = 1; i < N; i++)
            {
                I += func(xa + i*dx, args);
                if(errf(args))
                    return 0.0;
            }
            return I * dx;
        }

        double simp(
            std::function<double(double, void *)> func, 
            double xa, 
            double xb, 
            int N,
            void *args, 
            std::function<int(void *)> errf
        )
        {
            if(N % 2 == 1)
                N -= 1;

            double dx = (xb - xa) / N;
            int i;
            const double I1 = func(xa, args) + func(xb, args);
            if(errf(args))
                return 0.0;

            double I2 = 0.0;
            for(i=1; i<N; i+=2)
            {
                I2 += func(xa + i*dx, args);
                if(errf(args))
                    return 0.0;
            }

            double I3 = 0.0;
            for(i=2; i<N; i+=2)
            {
                I3 += func(xa + i*dx, args);
                if(errf(args))
                    return 0.0;
            }
            return (I1 + 4*I2 + 2*I3) * dx / 3.0;
        }

        double romb(
            std::function<double(double, void *)> func,
            double xa, 
            double xb, 
            int N,
            double atol, 
            double rtol, 
            void *args, 
            int *Neval, 
            double *eps,
            int verbose, 
            std::function<int (void *)> errf,
            double *pfa, 
            double *pfb)
        {
            double R[kmax] ;

            int m, k, k0, Nk;
            long fpm;
            double hk, Rp, err;

            double maxFracChange = 0.1;

            double fa, fb;
            if(!pfa)
            {
                fa = func(xa, args);
                if(errf(args))
                    return 0.0;
            }
            else
                fa = *pfa;

            if(!pfb)
            {
                fb = func(xb, args);
                if(errf(args))
                    return 0.0;
            }
            else
                fb = *pfb;

            hk = xb - xa;
            Nk = 1;
            R[kmax - 1] = 0.5*(xb-xa)*(fa + fb);

            if(Neval)
                *Neval = 2;

            for(k=1; k<kmax;  k++)
            {
                k0 = kmax - k - 1;
                hk *= 0.5;
                Nk *= 2;

                Rp = 0.0;
                for(m=1; m<Nk; m+=2)
                {
                    Rp += func(xa + m*hk, args);
                    if(errf(args))
                        return 0.0;
                }
                R[k0] = 0.5*R[k0+1] + hk*Rp;
                if(Neval != NULL)
                    *Neval += Nk/2;

                double lastVal = R[kmax - 1];
                fpm = 1;
                for(m=1; m<=k; m++)
                {
                    fpm *= 4;
                    err = (R[k0+m-1] - R[k0+m]) / (fpm - 1);
                    R[k0+m] = (fpm*R[k0+m-1] - R[k0+m]) / (fpm - 1);
                }
                //err = (R[kmax - 1] - R[0]) / (fpm - 1);

                double fracChange = fabs((R[kmax - 1] - lastVal) / lastVal);

                if(eps != NULL)
                    *eps = err;
                if(verbose)
                    printf("level %d:  Neval=%d  I=%.6lg  fracDelta=%.3lg"
                        " err=%.6lg  tol=%.6lg\n",
                            k, Nk+1, R[kmax - 1], fracChange, err, atol+rtol*fabs(R[0]));

                //printf("      k%d: I=%.6le err=%.6le frac=%.3le\n", k, R[0], err, 
                //        fabs(err) / (atol + rtol*fabs(R[0])));

                if((fabs(err) < atol + rtol*fabs(R[kmax - 1]))
                        && fracChange < maxFracChange)
                    break;

                if(N > 1 && Nk >= N)
                    break;
            }

            return R[kmax - 1];
        }

        double simp_adapt(
            std::function<double(double, void*)> func,
            double xa,
            double xb,
            int Nmax,
            double atol,
            double rtol, 
            void *args, 
            int *Neval,
            double *eps, 
            Mesh5 *mout, 
            int verbose,
            std::function<int(void*)> errf,
            double *pfa, 
            double *pfb
        )
        {
            double I = m5_adapt(func, xa, xb, Nmax, simpInitInterval,
                                simpProcessInterval, simpSplitInterval,
                                atol, rtol, args, Neval, eps, mout, verbose, errf,
                                pfa, pfb);
            return I;
        }

        double hybrid_adapt(
            std::function<double(double, void*)> func,
            double xa, 
            double xb,
            int Nmax,
            double atol, 
            double rtol, 
            void *args, 
            int *Neval,
            double *eps, 
            int verbose,
            std::function<int(void *)> errf,
            double *pfa, 
            double *pfb
        )
        {
            double I, fa, fb;

            //Compute f at endpoints
            if(!pfa)
            {
                fa = func(xa, args);
                if(errf(args))
                    return 0.0;
            }
            else
                fa = *pfa;

            if(!pfb)
            {
                fb = func(xb, args);
                if(errf(args))
                    return 0.0;
            }
            else
                fb = *pfb;

            double ratio = std::abs(fa/fb);

            double NLrtol = 1.0; //9.0e-5;

            // If there's a large gradient, use an adaptive scheme.
            if(ratio > 1e6|| ratio < 1e-6)
            {
                //Adaptive Simpson's rule is more efficient but less robust,
                // requires a small error tolerance
                if(rtol < NLrtol)
                    I = simp_adapt(func, xa, xb, Nmax, atol, rtol, args, Neval, eps,
                                NULL, verbose, errf, &fa, &fb);
                //The non-linear scheme is more robust, but less efficient. Useful
                //for large tolerances, where Simpson's rule under-estimates the error
                else
                    I = trapNL_adapt(func, xa, xb, Nmax, atol, rtol, args, Neval, eps,
                                    NULL, verbose, errf, &fa, &fb);
            }
            else
                //If the gradient is not too big, Romberg will converge the fastest.
                I = romb(func, xa, xb, Nmax, atol, rtol, args, Neval, eps, verbose, errf,
                        &fa, &fb);

            return I;
        }

        double cadre_adapt(
            std::function<double(double, void*)> func,
            double xa, 
            double xb,
            int Nmax,
            double atol, 
            double rtol, 
            void *args, 
            int *Neval,
            double *eps, 
            int verbose,
            std::function<int(void *)> errf,
            double *pfa, 
            double *pfb
        )
        {
            double I = m9_adapt(func, xa, xb, Nmax, cadreInitInterval,
                                cadreProcessInterval, cadreSplitInterval,
                                atol, rtol, args, Neval, eps, nullptr, verbose, errf,
                                pfa, pfb);
            return I;
        }

        double gk49_adapt(
            std::function<double(double, void*)> func,
            double xa, 
            double xb,
            int Nmax,
            double atol, 
            double rtol, 
            void *args, 
            int *Neval,
            double *eps, 
            int verbose,
            std::function<int(void *)> errf
        )
        {
            double I = m_adapt(func, xa, xb, Nmax,
                                gk49ProcessInterval, gk49SplitInterval,
                                atol, rtol, args, Neval, eps, NULL, verbose, errf);
            return I;
        }

        double gk715_adapt(
            std::function<double(double, void*)> func,
            double xa, 
            double xb,
            int Nmax,
            double atol, 
            double rtol, 
            void *args, 
            int *Neval,
            double *eps, 
            int verbose,
            std::function<int(void *)> errf
        )
        {
            double I = m_adapt(func, xa, xb, Nmax,
                                gk715ProcessInterval, gk715SplitInterval,
                                atol, rtol, args, Neval, eps, NULL, verbose, errf);
            return I;
        }

        double gk1021_adapt(
            std::function<double(double, void*)> func,
            double xa, 
            double xb,
            int Nmax,
            double atol, 
            double rtol, 
            void *args, 
            int *Neval,
            double *eps, 
            int verbose,
            std::function<int(void *)> errf
        )
        {
            double I = m_adapt(func, xa, xb, Nmax,
                                gk1021ProcessInterval, gk1021SplitInterval,
                                atol, rtol, args, Neval, eps, NULL, verbose, errf);
            return I;
        }

        double trap_adapt(
            std::function<double(double, void*)> func,
            double xa, 
            double xb, 
            int Nmax,
            double atol, 
            double rtol, 
            void *args, 
            int *Neval,
            double *eps, 
            Mesh3 *mout, 
            int verbose,
            std::function<int(void *)> errf, 
            double *pfa, 
            double *pfb
        )
        {
            double I = m3_adapt(func, xa, xb, Nmax, nullptr,
                                nullptr, nullptr,
                                atol, rtol, args, Neval, eps, mout, verbose, errf,
                                pfa, pfb);
            return I;
        }

        double trapNL_adapt(
            std::function<double(double, void*)> func,
            double xa, 
            double xb,
            int Nmax,
            double atol, 
            double rtol, 
            void *args, 
            int *Neval,
            double *eps, 
            Mesh5 *mout, 
            int verbose,
            std::function<int(void *)> errf,
            double *pfa, 
            double *pfb
        )
        {
            double I = m5_adapt(func, xa, xb, Nmax, trapNLInitInterval,
                                trapNLProcessInterval, trapNLSplitInterval,
                                atol, rtol, args, Neval, eps, mout, verbose, errf,
                                pfa, pfb);
            return I;
        }

        int trapInitInterval(
            std::function<double(double, void *)> func,
            void *args,
            Interval3 &i,
            std::function<int(void *)> errf,
            double *pfa,
            double *pfb
        )
        {
            if(!pfa)
            {
                i.fa = func(i.a, args);
                if(errf(args))
                    return 1;
            }
            else
                i.fa = *pfa;

            if(!pfb)
            {
                i.fb = func(i.b, args);
                if(errf(args))
                    return 2;
            }
            else
                i.fb = *pfb;

            return 2;
        }

        int trapProcessInterval(
            std::function<double(double, void *)> func,
            void *args,
            Interval3 &i,
            std::function<int(void *)> errf
        )
        {
            double fa = i.fa;
            double fb = i.fb;
            double fm = func(0.5*(i.a+i.b), args);
            if(errf(args))
                return 1;
            i.fm = fm;

            double h = 0.5*(i.b - i.a);

            double I0 = 0.5*(fa + fb) * 2*h;
            double I1 = 0.5*(fa + 2*fm + fb) * h;

            double err = (I1 - I0) / 3.0;
            i.err = std::abs(err);
            i.I = I1 + err;

            return 1;
        }

        int trapSplitInterval(
            std::function<double(double, void *)> func,
            void *args,
            Interval3 &i0,
            Interval3 &i1,
            Interval3 &i2,
            std::function<int(void *)> errf
        )
        {
            double x = 0.5*(i0.a + i0.b);
            i1.a = i0.a;
            i1.b = x;
            i2.a = x;
            i2.b = i0.b;

            i1.fa = i0.fa;
            i1.fb = i0.fm;
            i2.fa = i0.fm;
            i2.fb = i0.fb;

            int n = 0;
            n += trapProcessInterval(func, args, i1, errf);
            n += trapProcessInterval(func, args, i2, errf);

            return n;
        }

        int simpInitInterval(
            std::function<double(double, void *)> func,
            void *args,
            Interval5 &i,
            std::function<int(void *)> errf,
            double *pfa,
            double *pfb
        )
        {
            if(!pfa)
            {
                i.fa = func(i.a, args);
                if(errf(args))
                    return 1;
            }
            else
                i.fa = *pfa;

            if(!pfb)
            {
                i.fb = func(i.b, args);
                if(errf(args))
                    return 2;
            }
            else
                i.fb = *pfb;

            i.fm = func(0.5*(i.a + i.b), args);

            return 3;
        }

        int simpProcessInterval(
            std::function<double(double, void *)> func,
            void *args,
            Interval5 &i,
            std::function<int(void *)> errf
        )
        {
            double fa = i.fa;
            double fb = i.fb;
            double fm = i.fm;
            double fl = func(0.75*i.a+0.25*i.b, args);
            if(errf(args))
                return 0.0;
            double fr = func(0.25*i.a+0.75*i.b, args);
            if(errf(args))
                return 0.0;
            i.fl = fl;
            i.fr = fr;

            double h = 0.25*(i.b - i.a);

            double I0 = 2*h * (fa + 4*fm + fb)/3.0;
            double I1 = h * (fa + 4*fl + 2*fm + 4*fr + fb)/3.0;

            double err = (I1 - I0) / 15.0;
            i.err = std::abs(err);
            i.I = I1 + err;

            return 2;
        }

        int simpSplitInterval(
            std::function<double(double, void *)> func,
            void *args,
            Interval5 &i0,
            Interval5 &i1,
            Interval5 &i2,
            std::function<int(void *)> errf
        )
        {
            double xm = 0.5*(i0.a + i0.b);
            i1.a = i0.a;
            i1.b = xm;
            i2.a = xm;
            i2.b = i0.b;

            i1.fa = i0.fa;
            i1.fm = i0.fl;
            i1.fb = i0.fm;
            i2.fa = i0.fm;
            i2.fm = i0.fr;
            i2.fb = i0.fb;

            int n = 0;
            n += simpProcessInterval(func, args, i1, errf);
            n += simpProcessInterval(func, args, i2, errf);

            return n;
        }

        int trapNLInitInterval(
            std::function<double(double, void *)> func,
            void *args,
            Interval5 &i,
            std::function<int(void *)> errf,
            double *pfa,
            double *pfb
        )
        {
            if(!pfa)
            {
                i.fa = func(i.a, args);
                if(errf(args))
                    return 1;
            }
            else
                i.fa = *pfa;

            if(!pfb)
            {
                i.fb = func(i.b, args);
                if(errf(args))
                    return 2;
            }
            else
                i.fb = *pfb;

            i.fm = func(0.5*(i.a + i.b), args);

            return 3;
        }

        int trapNLProcessInterval(
            std::function<double(double, void *)> func,
            void *args,
            Interval5 &i,
            std::function<int(void *)> errf
        )
        {
            double fa = i.fa;
            double fb = i.fb;
            double fm = i.fm;
            double fl = func(0.75*i.a+0.25*i.b, args);
            if(errf(args))
                return 0.0;
            double fr = func(0.25*i.a+0.75*i.b, args);
            if(errf(args))
                return 0.0;
            i.fl = fl;
            i.fr = fr;

            double h = 0.25*(i.b - i.a);

            // First compute three trapezoid rule refinements
            // Assume error ~ epsilon * h^n, if f is smooth
            // then n = 2, but if we are under-resolved the effective n
            // could be less than two.
            double I0 = 2*h * (fa + fb);
            double I1 = h * (fa + 2*fm + fb);
            double I2 = 0.5*h * (fa + 2*(fl+fm+fr) + fb);

            double R = (I1 - I0) / (I2 - I1);
            double err = - (I2-I1)*(I2-I1) / (I2 - 2*I1 + I0);

            double Is0 = 2*h * (fa + 4*fm + fb)/3.0;
            double Is1 = h * (fa + 4*fl + 2*fm + 4*fr + fb)/3.0;
            double errs = (Is1 - Is0) / 15.0;

            /*
            double exact = (atan((i.b-((double *)args)[1])/((double *)args)[0])
                            - atan((i.a-((double *)args)[1])/((double *)args)[0]));
            printf("%lf\n", twopn);
            printf("   %lf %lf\n", i.a, i.b);
            printf("   %lf %lf %lf %lf %lf\n", fa, fl, fm, fr, fb);
            printf("   %lf %lf %lf\n", I0, I1, I2);
            printf("     %lf %le   (%lf, %le)\n", I2+err, err,
                    exact, std::abs(I2+err-exact));
            printf("   %lf %lf\n", Is0, Is1);
            printf("     %lf %le   (%lf, %le)\n", Is1+errs, errs,
                    exact, std::abs(I2+err-exact));
            */
            
            i.err = std::abs(err);
            i.I = I2 + err;

            if(R > 3.95 && R < 4.05)
            {
                i.err = std::abs(errs);
                i.I = Is1 + errs;
            }
            if(R > 4.05 || R < 1.95 || R != R)
            {
                double errt = (I2-I1)/3.0;
                i.err = std::abs(errt);
                i.I = I2 + errt;
            }

            return 2;
        }

        int trapNLSplitInterval(
            std::function<double(double, void *)> func,
            void *args,
            Interval5 &i0,
            Interval5 &i1,
            Interval5 &i2,
            std::function<int(void *)> errf
        )
        {
            double xm = 0.5*(i0.a + i0.b);
            i1.a = i0.a;
            i1.b = xm;
            i2.a = xm;
            i2.b = i0.b;

            i1.fa = i0.fa;
            i1.fm = i0.fl;
            i1.fb = i0.fm;
            i2.fa = i0.fm;
            i2.fm = i0.fr;
            i2.fb = i0.fb;

            int n = 0;
            n += trapNLProcessInterval(func, args, i1, errf);
            n += trapNLProcessInterval(func, args, i2, errf);

            return n;
        }

        int cadreInitInterval(
            std::function<double(double, void *)> func,
            void *args,
            Interval9 &i,
            std::function<int(void *)> errf,
            double *pfa,
            double *pfb
        )
        {
            int count = 0;
            if(!pfa)
            {
                i.fa = func(i.a, args);
                count++;
                if(errf(args))
                    return count;
            }
            else {
                i.fa = *pfa;
            }

            if(!pfb)
            {
                i.fb = func(i.b, args);
                count++;
                if(errf(args))
                    return count;
            }
            else
                i.fb = *pfb;

            i.fm = func(0.5*(i.a + i.b), args);
            count++;
            if(errf(args))
                return count;
            /*
            i.fl = func(0.75*i.a + 0.25*i.b, args);
            count++;
            if(errf(args))
                return count;
            i.fr = func(0.25*i.a + 0.75*i.b, args);
            count++;
            */

            return count;
        }

        int cadreProcessInterval(
            std::function<double(double, void *)> func,
            void *args,
            Interval9 &i,
            std::function<int(void *)> errf
        )
        {
            //printf("(%.3le, %.3le) %.3le +/- %.3le    %d\n",
            //       i.a, i.b, i.I, i.err, i.refinement);

            double scaleTol = 0.2;

            int count = 0;

            if(i.refinement < 4)
            {
                //This interval is not being refined, so start with the trapezoid
                //approach with non-linear error
                double fa = i.fa;
                double fm = i.fm;
                double fb = i.fb;

                double fl, fr;

                if(i.refinement == 0)
                {
                    fl = func(0.75*i.a+0.25*i.b, args);
                    count++;
                    if(errf(args))
                        return count;
                    fr = func(0.25*i.a+0.75*i.b, args);
                    count++;
                    if(errf(args))
                        return count;
                    
                    i.fl = fl;
                    i.fr = fr;
                }
                else
                {
                    fl = i.fl;
                    fr = i.fr;
                    i.refinement = 0;
                }

                double h = 0.25*(i.b - i.a);

                // Compute three trapezoid rule refinements
                // Assume error ~ epsilon * h^n, if f is smooth
                // then n = 2, but if we are under-resolved the effective n
                // could be less than two.
                double T0 = 2*h * (fa + fb);
                double T1 = 0.5*T0 + 2*h*fm;
                double T2 = 0.5*T1 + h*(fl + fr);

                double R0 = (T1 - T0) / (T2 - T1);

                if(R0 > 4-scaleTol && R0 < 4+scaleTol)
                {
                    // We might be in the asymptotic regime!  
                    // Compute an extra refinement level to make sure.
                
                    double fll = func(.875*i.a+0.125*i.b, args);
                    count++;
                    if(errf(args))
                        return count;
                    double flr = func(.625*i.a+0.375*i.b, args);
                    count++;
                    if(errf(args))
                        return count;
                    double frl = func(.375*i.a+0.625*i.b, args);
                    count++;
                    if(errf(args))
                        return count;
                    double frr = func(.125*i.a+0.875*i.b, args);
                    count++;
                    if(errf(args))
                        return count;

                    i.fll = fll;
                    i.flr = flr;
                    i.frl = frl;
                    i.frr = frr;
                    i.refinement = 1;
                    
                    double T3 = 0.5*T2 + 0.5*h*(fll + flr + frl + frr);
                    double R1 = (T2 - T1) / (T3 - T2);

                    if(R1 > 4-scaleTol && R1 < 4+scaleTol)
                    {
                        // We're in the asymptotic regime!
                        //
                        // Time to accelerate with
                        // Romberg.
                        // The interval entries "fx" are now going to store the
                        // Romberg tableau, and "refinement" the current refinement
                        // level.
                        //
                        // T0
                        // T1 R11
                        // T2 R21 R22
                        // T3 R31 R32 R33

                        double R11 = (4*T1 - T0) / 3.0;
                        double R21 = (4*T2 - T1) / 3.0;
                        double R22 = (16*R21 - R11) / 15.0;
                        i.fa = T3;
                        i.fll = (4*T3 - T2) / 3.0;
                        i.fl = (16*i.fll - R21) / 15.0;
                        i.flr = (64*i.fl - R22) / 63.0;

                        i.I = i.flr;
                        i.err = std::abs((i.fl - R22) / 63.0);
                        i.refinement = 4;

                        return count;
                    }
                    else
                    {
                        // Not asymptotic yet, just use trapezoid rule for now.
                        double errt = (T3-T2)/3.0;
                        i.err = std::abs(errt);
                        i.I = T3 + errt;

                        return count;
                    }
                }
                else if(R0 >= 4+scaleTol || R0 < 2 || R0 != R0)
                {
                    // We're in a strange regime.  Proceed with caution and
                    // just use vanilla Trapezoid rule with over-estimated error.
                    
                    double err10 = std::abs(T1-T0);
                    double err21 = std::abs(T2-T1);
                    double err = err10 > err21 ? err10 : err21;
                    i.err = err;
                    i.I = T2;//  + errt;
                    
                    return count;
                }
                else
                {
                    // 1 <~ n <~ 2
                    // We are near the n~1 regime, implying we're under resolved.
                    // Use the non-linear estimator.

                    double err = - (T2-T1)*(T2-T1) / (T2 - 2*T1 + T0);

                    i.err = std::abs(err);
                    if(std::abs(T2 + err) < 1.0e-14 * std::abs(T2))
                        i.I = 1.0e-14 * T2;
                    else
                        i.I = T2 + err;

                    return count;
                }

                return count;
            }

            // We're in Romberg mode! From now on this interval will not be split,
            // but will successively refine its estimate using romberg integration.
            // The Romberg tableau is stored in the interval's fx entries.
            
            if(i.refinement < 0)
                return 0;

            int n;
            int nint = 1;
            for(n=0; n<i.refinement; n++)
                nint *= 2;
            double h = (i.b - i.a) / nint;

            double sum = 0;
            for(n=1; n<nint; n += 2)
            {
                sum += func(i.a + h*n, args);
                count++;
                if(errf(args))
                    return count;
            }
        
            // m = 1
            double T = 0.5*i.fa + h*sum;
            double Rnm = i.fa;
            double Rnmp = i.fll;
            i.fa = T;
            i.fll = (4.0*i.fa - Rnm) / 3.0;
            if(i.refinement == 1)
            {
                i.I = i.fll;
                i.err = std::abs((i.fa - Rnm) / 3.0);
                i.refinement += 1;
                return count;
            }

            // m = 2
            Rnm = Rnmp;
            Rnmp = i.fl;
            i.fl = (16 * i.fll - Rnm) / 15.0;
            if(i.refinement == 2)
            {
                i.I = i.fl;
                i.err = std::abs((i.fll - Rnm) / 15.0);
                i.refinement += 1;
                return count;
            }

            // m = 3
            Rnm = Rnmp;
            Rnmp = i.flr;
            i.flr = (64 * i.fl - Rnm) / 63.0;
            if(i.refinement == 3)
            {
                i.I = i.flr;
                i.err = std::abs((i.fl - Rnm) / 63.0);
                i.refinement += 1;
                return count;
            }

            // m = 4
            Rnm = Rnmp;
            Rnmp = i.fm;
            i.fm = (256 * i.flr - Rnm) / 255.0;
            if(i.refinement == 4)
            {
                i.I = i.fm;
                i.err = std::abs((i.flr - Rnm) / 255.0);
                i.refinement += 1;
                return count;
            }

            // m = 5
            Rnm = Rnmp;
            Rnmp = i.frl;
            i.frl = (1024 * i.fm - Rnm) / 1023.0;
            if(i.refinement == 5)
            {
                i.I = i.frl;
                i.err = std::abs((i.fm - Rnm) / 1023.0);
                i.refinement += 1;
                return count;
            }

            // m = 6
            Rnm = Rnmp;
            Rnmp = i.fr;
            i.fr = (4096 * i.frl - Rnm) / 4095.0;
            if(i.refinement == 6)
            {
                i.I = i.fr;
                i.err = std::abs((i.frl - Rnm) / 4095.0);
                i.refinement += 1;
                return count;
            }

            // m = 7
            Rnm = Rnmp;
            Rnmp = i.frr;
            i.frr = (16384 * i.fr - Rnm) / 16383.0;
            if(i.refinement == 7)
            {
                i.I = i.frr;
                i.err = std::abs((i.fr - Rnm) / 16383.0);
                i.refinement += 1;
                return count;
            }

            // m = 8
            Rnm = Rnmp;
            Rnmp = i.fb;
            i.fb = (65536 * i.frr - Rnm) / 65535.0;
            if(i.refinement == 8)
            {
                i.I = i.fb;
                i.err = std::abs((i.frr - Rnm) / 65535.0);
                i.refinement += 1;
                return count;
            }

            // m = 9
            Rnm = Rnmp;
            i.I = (262144 * i.fb - Rnm) / 262143.0;
            i.err = std::abs((i.fb - Rnm) / 262143.0);
            i.refinement += 1;

            return count;
        }

        int cadreSplitInterval(
            std::function<double(double, void *)> func,
            void *args,
            Interval9 &i0,
            Interval9 &i1,
            Interval9 &i2,
            std::function<int(void *)> errf
        )
        {
            double xm = 0.5*(i0.a + i0.b);
            i1.a = i0.a;
            i1.b = xm;
            i2.a = xm;
            i2.b = i0.b;

            i1.fa = i0.fa;
            i1.fl = i0.fll;
            i1.fm = i0.fl;
            i1.fr = i0.flr;
            i1.fb = i0.fm;
            i1.refinement = i0.refinement;
            i2.fa = i0.fm;
            i2.fl = i0.frl;
            i2.fm = i0.fr;
            i2.fr = i0.frr;
            i2.fb = i0.fb;
            i2.refinement = i0.refinement;

            int n = 0;
            n += cadreProcessInterval(func, args, i1, errf);
            n += cadreProcessInterval(func, args, i2, errf);

            return n;
        }

        int gk49ProcessInterval(
            std::function<double(double, void *)> func,
            void *args,
            Interval &i,
            std::function<int(void *)> errf
        )
        {
            /*
            * z(a) = -1, z(b) = +1
            * z = 2/(b-a) * (x - a) - 1
            * x = (b-a)/2 * (z + 1) + a = (b-a)/2 * z + (a+b)/2
            * dx = (b-a)/2 * dz
            */

            double c = 0.5*(i.b - i.a);
            double z0 = 0.5*(i.a + i.b);

            int count = gk_compute(func, args, errf, c, z0, g4x, k9x, g4w, k9w, 4,
                                &(i.I), &(i.err));

            return count;
        }

        int gk49SplitInterval(
            std::function<double(double, void *)> func,
            void *args,
            Interval &i0,
            Interval &i1,
            Interval &i2,
            std::function<int(void *)> errf
        )
        {
            double x = 0.5*(i0.a + i0.b);
            i1.a = i0.a;
            i1.b = x;
            i2.a = x;
            i2.b = i0.b;

            int n = 0;
            n += gk49ProcessInterval(func, args, i1, errf);
            if(errf(args))
                return n;
            n += gk49ProcessInterval(func, args, i2, errf);

            return n;
        }

        int gk715ProcessInterval(
            std::function<double(double, void *)> func,
            void *args,
            Interval &i,
            std::function<int(void *)> errf
        )
        {
            /*
            * z(a) = -1, z(b) = +1
            * z = 2/(b-a) * (x - a) - 1
            * x = (b-a)/2 * (z + 1) + a = (b-a)/2 * z + (a+b)/2
            * dx = (b-a)/2 * dz
            */

            double c = 0.5*(i.b - i.a);
            double z0 = 0.5*(i.a + i.b);

            int count = gk_compute(func, args, errf, c, z0, g7x, k15x, g7w, k15w, 7,
                                &(i.I), &(i.err));

            return count;
        }

        int gk715SplitInterval(
            std::function<double(double, void *)> func,
            void *args,
            Interval &i0,
            Interval &i1,
            Interval &i2,
            std::function<int(void *)> errf
        )
        {
            double x = 0.5*(i0.a + i0.b);
            i1.a = i0.a;
            i1.b = x;
            i2.a = x;
            i2.b = i0.b;

            int n = 0;
            n += gk715ProcessInterval(func, args, i1, errf);
            if(errf(args))
                return n;
            n += gk715ProcessInterval(func, args, i2, errf);

            return n;
        }

        int gk1021ProcessInterval(
            std::function<double(double, void *)> func,
            void *args,
            Interval &i,
            std::function<int(void *)> errf
        )
        {

            double c = 0.5*(i.b - i.a);
            double z0 = 0.5*(i.a + i.b);

            int count = gk_compute(func, args, errf, c, z0, g10x, k21x, g10w, k21w, 10,
                                &(i.I), &(i.err));

            return count;
        }

        int gk1021SplitInterval(
            std::function<double(double, void *)> func,
            void *args,
            Interval &i0,
            Interval &i1,
            Interval &i2,
            std::function<int(void *)> errf
        )
        {
            double x = 0.5*(i0.a + i0.b);
            i1.a = i0.a;
            i1.b = x;
            i2.a = x;
            i2.b = i0.b;

            int n = 0;
            n += gk1021ProcessInterval(func, args, i1, errf);
            if(errf(args))
                return n;
            n += gk1021ProcessInterval(func, args, i2, errf);

            return n;
        }

        double m_adapt(
            std::function<double(double, void *)> func,
            double xa, 
            double xb, 
            int Nmax,
            std::function<int(std::function<double(double, void*)> func, void *, mesh::Interval&, std::function<int(void *)> errf)> processInterval,
            std::function<int(std::function<double(double, void*)> func, void *, mesh::Interval&, mesh::Interval&, mesh::Interval&, std::function<int(void *)> errf)> splitInterval,
            double atol, 
            double rtol, 
            void *args, 
            int *Neval,
            double *eps, 
            mesh::Mesh *mout, 
            int verbose,
            std::function<int(void *)> errf
        )
        {
            Mesh m;
            Interval i = {xa, xb, 0, 0};
            
            int n = processInterval(func, args, i, errf);
            if(errf(args))
            {
                return 0.0;
            }
            m.insert(i);

            double I = i.I;
            double err = i.err;
            int num_intervals = 1;
            int last_check = 1;

            while(n < Nmax && err >= atol + rtol*std::abs(I))
            {
                m.extract(i);
                Interval i1;
                Interval i2;
                n += splitInterval(func, args, i, i1, i2, errf);
                if(errf(args))
                {
                    return 0.0;
                }
                m.insert(i1);
                m.insert(i2);
                num_intervals++;

                if(num_intervals == 2*last_check)
                {
                    err        = m.totalError();
                    I          = m.totalIntergral();
                    last_check = num_intervals;
                }
                else
                {
                    err += i1.err + i2.err - i.err;
                    I += i1.I + i2.I - i.I;
                }
                
                if(verbose)
                    printf("Num Intervals: %d - I=%.12lg  err=%.3lg  tol=%.3lg"
                        "  meshOk=%d\n",
                        num_intervals, I, err, atol + rtol*std::abs(I), m.meshCheck());
            }

            I = m.totalIntergral();
            if(Neval)
                *Neval = n;

            if(eps)
            {
                err = m.totalError(); 
                *eps = err;
            }

            if(mout)
                *mout = m;

            return I;
        }


        double m3_adapt(
            std::function<double(double, void *)> func,
            double xa, 
            double xb, 
            int Nmax,
            std::function<int(std::function<double(double, void*)> func, void *, Interval3 &, std::function<int(void *)> errf, double *pfa, double *pfb)> initInterval,
            std::function<int(std::function<double(double, void*)> func, void *, Interval3 &, std::function<int(void *)> errf)> processInterval,
            std::function<int(std::function<double(double, void*)> func, void *, Interval3 &, Interval3 &, Interval3 &, std::function<int(void *)> errf)> splitInterval,
            double atol, 
            double rtol, 
            void *args, 
            int *Neval,
            double *eps, 
            Mesh3 *mout, 
            int verbose,
            std::function<int(void *)> errf,
            double *pfa, 
            double *pfb
        )
        {
            Mesh3 m;
            Interval3 i = {.a=xa, .b=xb, .I=0, .err=0, .fa=0, .fb=0, .fm=0};
            
            int n = initInterval(func, args, i, errf, pfa, pfb);
            if(errf(args))
            {
                return 0.0;
            }
            
            n += processInterval(func, args, i, errf);
            if(errf(args))
            {
                return 0.0;
            }

            m.insert(i);
            double I          = i.I;
            double err        = i.err;
            int num_intervals = 1;
            int last_check    = 1;

            while(n < Nmax && err >= atol + rtol*std::abs(I))
            {
                m.extract(i);
                Interval3 i1;
                Interval3 i2;
                n += splitInterval(func, args, i, i1, i2, errf);
                if(errf(args))
                {
                    return 0.0;
                }
                m.insert(i1);
                m.insert(i2);
                num_intervals++;

                if(num_intervals == 2*last_check)
                {
                    err        = m.totalError();
                    I          = m.totalIntergral();
                    last_check = num_intervals;
                }
                else
                {
                    err += i1.err + i2.err - i.err;
                    I   += i1.I + i2.I - i.I;
                }
                
                if(verbose)
                    printf("Num Intervals: %d - I=%.12lg  err=%.3lg  tol=%.3lg"
                        "  meshOk=%d\n",
                        num_intervals, I, err, atol + rtol*std::abs(I), m.meshCheck());
            }
            I = m.totalIntergral();

            if(Neval)
                *Neval = n;

            if(eps)
            {
                err = m.totalError();
                *eps = err;
            }

            if(mout)
                *mout = m;

            return I;
        }

        double m5_adapt(
            std::function<double(double, void *)> func,
            double xa, 
            double xb, 
            int Nmax,
            std::function<int(std::function<double(double, void*)> func, void *, Interval5 &, std::function<int(void *)> errf, double *pfa, double *pfb)> initInterval,
            std::function<int(std::function<double(double, void*)> func, void *, Interval5 &, std::function<int(void *)> errf)> processInterval,
            std::function<int(std::function<double(double, void*)> func, void *, Interval5 &, Interval5 &, Interval5 &, std::function<int(void *)> errf)> splitInterval,
            double atol, 
            double rtol, 
            void *args, 
            int *Neval,
            double *eps, 
            Mesh5 *mout, 
            int verbose,
            std::function<int(void *)> errf,
            double *pfa, 
            double *pfb
        )
        {
            Mesh5 m;
            Interval5 i = {.a=xa, .b=xb, .I=0, .err=0,
                        .fa=0, .fb=0, .fl=0, .fm=0, .fr=0};

            int n = initInterval(func, args, i, errf, pfa, pfb);
            if(errf(args))
            {
                return 0.0;
            }

            n += processInterval(func, args, i, errf);
            if(errf(args))
            {
                return 0.0;
            }
            m.insert(i);

            double I          = i.I;
            double err        = i.err;
            int num_intervals = 1;
            int last_check    = 1;

            while(n < Nmax && err > atol + rtol*std::abs(I))
            {
                m.extract(i);
                Interval5 i1;
                Interval5 i2;
                n += splitInterval(func, args, i, i1, i2, errf);
                if(errf(args))
                {
                    return 0.0;
                }
                m.insert(i1);
                m.insert(i2);
                num_intervals++;

                if(num_intervals == 2*last_check)
                {
                    err = m.totalError(); 
                    I   = m.totalIntergral(); 
                    last_check = num_intervals;
                }
                else
                {
                    err += i1.err + i2.err - i.err;
                    I += i1.I + i2.I - i.I;
                }
                
                if(verbose)
                    printf("Num Intervals: %d - I=%.12lg  err=%.3lg  tol=%.3lg"
                        "  meshOk=%d\n",
                        num_intervals, I, err, atol + rtol*std::abs(I), m.meshCheck());
            }

            I = m.totalIntergral();

            if(Neval)
                *Neval = n;

            if(eps)
            {
                err = m.totalError();
                *eps = err;
            }

            if(mout)
                *mout = m;

            return I;
        }

        double m9_adapt(
            std::function<double(double, void *)> func,
            double xa, 
            double xb, 
            int Nmax,
            std::function<int(std::function<double(double, void*)> func, void *, Interval9 &, std::function<int(void *)> errf, double *pfa, double *pfb)> initInterval,
            std::function<int(std::function<double(double, void*)> func, void *, Interval9 &, std::function<int(void *)> errf)> processInterval,
            std::function<int(std::function<double(double, void*)> func, void *, Interval9 &, Interval9 &, Interval9 &, std::function<int(void *)> errf)> splitInterval,
            double atol, 
            double rtol, 
            void *args, 
            int *Neval,
            double *eps, 
            Mesh9 *mout, 
            int verbose,
            std::function<int(void *)> errf,
            double *pfa, 
            double *pfb
        )
        {
            Mesh9 m;
            Interval9 i = {.a=xa, .b=xb, .I=0, .err=0,
                        .fa=0, .fll=0, .fl=0, .flr=0, .fm=0,
                        .frl=0, .fr=0, .frr=0, .fb=0, .refinement=0};
            int n = initInterval(func, args, i, errf, pfa, pfb);
            if(errf(args))
            {
                return 0.0;
            }

            n += processInterval(func, args, i, errf);
            if(errf(args))
            {
                return 0.0;
            }
            m.insert(i);

            double I           = i.I;
            double err         = i.err;
            int num_intervals  = 1;
            int num_iterations = 0;
            int last_check     = 0;

            if(verbose)
            {   
                printf("   Num Intervals: %d - I=%.12lg  err=%.3lg  tol=%.3lg"
                    "  meshOk=%d\n",
                    num_intervals, I, err, atol + rtol*std::abs(I), m.meshCheck());
                m.intervalWrite(i, stdout);
            }

            while(n < Nmax && err > atol + rtol*std::abs(I))
            {
                m.extract(i);

                if(verbose && num_iterations > 0)
                    m.intervalWrite(i, stdout);

                double newI, newErr;

                if(i.refinement >= 4)
                {
                    double oldI = i.I;
                    double olderr = i.err;

                    n += processInterval(func, args, i, errf);
                    if(errf(args))
                    {
                        return 0.0;
                    }
                    m.insert(i);

                    err += i.err - olderr;
                    I += i.I - oldI;

                    newI = i.I;
                    newErr = i.err;
                }
                else
                {
                    Interval9 i1;
                    Interval9 i2;
                    n += splitInterval(func, args, i, i1, i2, errf);
                    if(errf(args))
                    {
                        return 0.0;
                    }
                    m.insert(i1);
                    m.insert(i2);
                    num_intervals++;
                    
                    err += i1.err + i2.err - i.err;
                    I += i1.I + i2.I - i.I;
                    
                    newI = i1.I + i2.I;
                    newErr = i1.err + i2.err;
                }

                if(num_iterations == 2*last_check)
                {
                    err = m.totalError(); //(&m);
                    I   = m.totalIntergral(); //mesh9TotalIntegral(&m);
                    last_check = num_intervals;
                }
                
                if(verbose)
                {
                    printf("                   ---> %.12le +/- %.3le\n", newI, newErr);
                    printf("   Num Intervals: %d - I=%.12lg  err=%.3lg  tol=%.3lg"
                        "  meshOk=%d\n",
                        num_intervals, I, err, atol + rtol*std::abs(I), m.meshCheck());
                }
                num_iterations++;
            }

            // std::cout << m.N << ", " << Nmax << ", " << err << "\n";
            // std::cin.get();
            I = m.totalIntergral();

            if(Neval)
                *Neval = n;

            if(eps)
            {
                err  = m.totalError();
                *eps = err;
            }

            if(mout)
                *mout = m;
            
            std::cout << I << ", " << Nmax << "\n";
            return I;
        }

        int gk_compute(
            std::function<double(double, void *)> func,
            void *args,
            std::function<int(void*)> errf,
            double c,
            double z0, 
            const double xg[],
            const double xk[],
            const double wg[],
            const double wgk[],
            int ng,
            double *I,
            double *err
        )
        {
            /*
            * z(a) = -1, z(b) = +1
            * z = 2/(b-a) * (x - a) - 1
            * x = (b-a)/2 * (z + 1) + a = (b-a)/2 * z + (a+b)/2
            * dx = (b-a)/2 * dz
            */

            int count = 0;

            int nk = ng + 1;

            double fg[max_gauss_n];
            double fk[max_gauss_n];

            int n;

            for(n=0; n<ng; n++)
            {
                fg[n] = func(c * xg[n] + z0, args);
                count++;
                if(errf(args))
                    return count;
            }
            for(n=0; n<nk; n++)
            {
                fk[n] = func(c * xk[n] + z0, args);
                count++;
                if(errf(args))
                    return count;
            }

            double Ig = 0;
            double Igk = 0;
            double Iegk = 0;
            for(n=0; n<ng; n++)
            {
                Ig += wg[n] * fg[n];
                Igk += wgk[2*n+1] * fg[n];
            }
            for(n=0; n<nk; n++)
                Igk += wgk[2*n] * fk[n];

            double Kmean = 0.5 * Igk;
            for(n=0; n<ng; n++)
                Iegk += wgk[2*n+1] * std::abs(fg[n] - Kmean);
            for(n=0; n<nk; n++)
                Iegk += wgk[2*n] * std::abs(fk[n] - Kmean);

            double err2 = Iegk * pow(200 * std::abs(Igk - Ig) / Iegk, 1.5);
            if(err2 > Iegk)
                err2 = Iegk;
            
            *I = c * Igk;
            *err = c * err2;

            //*err = std::abs(Ik15 - Ig7);

            return count;
        }

    } // namespace integrate
    
} // namespace afterglowpy
