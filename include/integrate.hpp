#ifndef INTEGRATE_HPP
#define INTEGRATE_HPP

#include "interval.hpp"
#include <functional>
#include <vector>
#include <iostream>

/*
 * Various routines for integrating 1D functions.
 * trap() and simp() are fixed stencil implementations of the Trapezoid Rule
 * and Simpson's Rule respectively.
 *
 * romb() is an adaptive Romberg integrator.
 *
 * trap_adapt() and simp_adapt() are globally adaptive integrators based on
 *    the Trapezoid and Simpson's rule, respectively. They successively bisect
 *    the integration domain into subintervals, prioritizing the subintervals
 *    with largest (estimated) error, until the total absolute error estimate
 *    is within tolerance.
 */

/*
 * Integration routines
 */
namespace afterglowpy
{
    namespace integrate
    {
        /**
         * @brief Calculate integral using trapezoidal rule
         *
         * @param func function or integrand
         * @param xa   starting integration boundary
         * @param xb   final integration boundary
         * @param N    number of trapezoids 
         * @param args extra number of function arguments
         * @param errf error function for tolerance on integration
         * @return double
         */
        double trap(
            std::function<double(double, void *)> func, 
            double xa, 
            double xb, 
            int N,
            void *args, 
            std::function<int(void *)> errf);

        /**
         * @brief Calculate integral using Simpson's rule
         *
         * @param func function or integrand
         * @param xa   starting integration boundary
         * @param xb   final integration boundary
         * @param N    number of trapezoids 
         * @param args extra number of function arguments
         * @param errf error function tolerance on integration
         * @return double
         */
        double simp(
            std::function<double(double, void *)> func, 
            double xa, 
            double xb, 
            int N,
            void *args, 
            std::function<int(void *)> errf);

        /**
         * @brief Romberg integration method
         * 
         * @param func  fumction or integrand
         * @param xa    initial integral boundary
         * @param xb    final integral boundary
         * @param N     number of intervals under integrand
         * @param atol  absolute tolerance
         * @param rtol  relative tolerance
         * @param args  extra args in func integrand
         * @param Neval number of samples
         * @param eps   
         * @param verbose 
         * @param errf 
         * @param pfa 
         * @param pfb 
         * @return double 
         */
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
            double *pfb);

        /**
         * @brief Adaptive trapezoidal rule
         * 
         * @param func  The function or integrand
         * @param xa    Initial inegral boundary
         * @param xb    Final integral boundary
         * @param Nmax  Max number of subdivisions  
         * @param atol  absolute tolerance
         * @param rtol  relative tolerance
         * @param args  extra args for func
         * @param Neval array of samples
         * @param eps   
         * @param mout 
         * @param verbose 
         * @param errf 
         * @param pfa 
         * @param pfb 
         * @return double 
         */
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
            struct mesh::Mesh3 *mout, 
            int verbose,
            std::function<int(void *)> errf, 
            double *pfa, 
            double *pfb);

        /**
         * @brief Adaptive Simpson's rule
         * 
         * @param func  The function or integrand
         * @param xa    Initial inegral boundary
         * @param xb    Final integral boundary
         * @param Nmax  Max number of subdivisions  
         * @param atol  absolute tolerance
         * @param rtol  relative tolerance
         * @param args  extra args for func
         * @param Neval array of samples
         * @param eps 
         * @param mout 
         * @param verbose 
         * @param errf 
         * @param pfa 
         * @param pfb 
         * @return double 
         */
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
            struct mesh::Mesh5 *mout, 
            int verbose,
            std::function<int(void*)> errf,
            double *pfa, 
            double *pfb);

        //=========================================================================
        // Hereafter, assume all param names serve same function as commented above
        //=========================================================================
        /**
         * @brief Adpative Trapezoidal integral with NL...
         * 
         * @param func 
         * @param xa 
         * @param xb 
         * @param Nmax 
         * @param atol 
         * @param rtol 
         * @param args 
         * @param Neval 
         * @param eps 
         * @param mout 
         * @param verbose 
         * @param errf 
         * @param pfa 
         * @param pfb 
         * @return double 
         */
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
            struct mesh::Mesh5 *mout, 
            int verbose,
            std::function<int(void *)> errf,
            double *pfa, 
            double *pfb);

        /**
         * @brief 
         * 
         * @param func
         * @param xa 
         * @param xb 
         * @param Nmax 
         * @param atol 
         * @param rtol 
         * @param args 
         * @param Neval 
         * @param eps 
         * @param verbose 
         * @param errf 
         * @param pfa 
         * @param pfb 
         * @return double 
         */
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
            double *pfb);

        /**
         * @brief Integral method named Cadre...
         * 
         * @param func 
         * @param xa 
         * @param xb 
         * @param Nmax 
         * @param atol 
         * @param rtol 
         * @param args 
         * @param Neval 
         * @param eps 
         * @param verbose 
         * @param errf 
         * @param pfa 
         * @param pfb 
         * @return double 
         */
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
            double *pfb);

        /**
         * @brief 
         * 
         * @param func 
         * @param xa 
         * @param xb 
         * @param Nmax 
         * @param atol 
         * @param rtol 
         * @param args 
         * @param Neval 
         * @param eps 
         * @param verbose 
         * @param errf 
         * @return double 
         */
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
            std::function<int(void *)> errf);

        /**
         * @brief 
         * 
         * @param func 
         * @param xa 
         * @param xb 
         * @param Nmax 
         * @param atol 
         * @param rtol 
         * @param args 
         * @param Neval 
         * @param eps 
         * @param verbose 
         * @param errf 
         * @return double 
         */
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
            std::function<int(void *)> errf);

        /**
         * @brief 
         * 
         * @param func
         * @param xa 
         * @param xb 
         * @param Nmax 
         * @param atol 
         * @param rtol 
         * @param args 
         * @param Neval 
         * @param eps 
         * @param verbose 
         * @param errf 
         * @return double 
         */
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
            std::function<int(void *)> errf);

        /*
        * Internal functions for trap_adapt and simp_adapt.
        */

        //=====================================================================================
        /**
         * @brief Internal function for trapezoidal and Simpsons adaptive methods
         * 
         * @param func 
         * @param xa 
         * @param xb 
         * @param Nmax 
         * @param processInterval 
         * @param splitInterval 
         * @param atol 
         * @param rtol 
         * @param args 
         * @param Neval 
         * @param eps 
         * @param mout 
         * @param verbose 
         * @param errf 
         * @return double 
         */
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
            std::function<int(void *)> errf);

        /**
         * @brief 
         * 
         * @param func 
         * @param xa 
         * @param xb 
         * @param Nmax 
         * @param initInterval 
         * @param processInterval 
         * @param splitInterval 
         * @param atol 
         * @param rtol 
         * @param args 
         * @param Neval 
         * @param eps 
         * @param mout 
         * @param verbose 
         * @param errf 
         * @param pfa 
         * @param pfb 
         * @return double 
         */
        double m3_adapt(
            std::function<double(double, void *)> func,
            double xa, 
            double xb, 
            int Nmax,
            std::function<int(std::function<double(double, void*)> func, void *, mesh::Interval3&, std::function<int(void *)> errf, double *pfa, double *pfb)> initInterval,
            std::function<int(std::function<double(double, void*)> func, void *, mesh::Interval3&, std::function<int(void *)> errf)> processInterval,
            std::function<int(std::function<double(double, void*)> func, void *, mesh::Interval3&, mesh::Interval3&, mesh::Interval3&, std::function<int(void *)> errf)> splitInterval,
            double atol, 
            double rtol, 
            void *args, 
            int *Neval,
            double *eps, 
            mesh::Mesh3 *mout, 
            int verbose,
            std::function<int(void *)> errf,
            double *pfa, 
            double *pfb);

        /**
         * @brief 
         * 
         * @param func 
         * @param xa 
         * @param xb 
         * @param Nmax 
         * @param initInterval 
         * @param processInterval 
         * @param splitInterval 
         * @param atol 
         * @param rtol 
         * @param args 
         * @param Neval 
         * @param eps 
         * @param mout 
         * @param verbose 
         * @param errf 
         * @param pfa 
         * @param pbf 
         * @return double 
         */
        double m5_adapt(
            std::function<double(double, void *)> func,
            double xa, 
            double xb, 
            int Nmax,
            std::function<int(std::function<double(double, void*)> func, void *, mesh::Interval5&, std::function<int(void *)> errf, double *pfa, double *pfb)> initInterval,
            std::function<int(std::function<double(double, void*)> func, void *, mesh::Interval5&, std::function<int(void *)> errf)> processInterval,
            std::function<int(std::function<double(double, void*)> func, void *, mesh::Interval5&, mesh::Interval5&, mesh::Interval5&, std::function<int(void *)> errf)> splitInterval,
            double atol, 
            double rtol, 
            void *args, 
            int *Neval,
            double *eps, 
            mesh::Mesh5 *mout, 
            int verbose,
            std::function<int(void *)> errf,
            double *pfa, 
            double *pfb);

        /**
         * @brief 
         * 
         * @param func 
         * @param xa 
         * @param xb 
         * @param Nmax 
         * @param initInterval 
         * @param processInterval 
         * @param splitInterval 
         * @param atol 
         * @param rtol 
         * @param args 
         * @param Neval 
         * @param eps 
         * @param mout 
         * @param verbose 
         * @param errf 
         * @param pfa 
         * @param pbf 
         * @return double 
         */
        double m9_adapt(
            std::function<double(double, void *)> func,
            double xa, 
            double xb, 
            int Nmax,
            std::function<int(std::function<double(double, void*)> func, void *, mesh::Interval9&, std::function<int(void *)> errf, double *pfa, double *pfb)> initInterval,
            std::function<int(std::function<double(double, void*)> func, void *, mesh::Interval9&, std::function<int(void *)> errf)> processInterval,
            std::function<int(std::function<double(double, void*)> func, void *, mesh::Interval9&, mesh::Interval9&, mesh::Interval9&, std::function<int(void *)> errf)> splitInterval,
            double atol, 
            double rtol, 
            void *args, 
            int *Neval,
            double *eps, 
            mesh::Mesh9 *mout, 
            int verbose,
            std::function<int(void *)> errf,
            double *pfa, 
            double *pfb);

        //=====================================================================================
        /**
         * @brief 
         * 
         * @param func 
         * @param args 
         * @param i 
         * @param errf 
         * @param pfa 
         * @param pfb 
         * @return int 
         */
        int trapInitInterval(
            std::function<double(double, void *)> func,
            void *args,
            mesh::Interval3 &i,
            std::function<int(void *)> errf,
            double *pfa,
            double *pfb);

        /**
         * @brief 
         * 
         * @param func 
         * @param args 
         * @param i 
         * @param errf 
         * @return int 
         */
        int trapProcessInterval(
            std::function<double(double, void *)> func,
            void *args,
            mesh::Interval3 &i,
            std::function<int(void *)> errf);

        /**
         * @brief 
         * 
         * @param func
         * @param args 
         * @param i0 
         * @param i1 
         * @param i2 
         * @param errf 
         * @return int 
         */
        int trapSplitInterval(
            std::function<double(double, void *)> func,
            void *args,
            mesh::Interval3 &i0,
            mesh::Interval3 &i1,
            mesh::Interval3 &i2,
            std::function<int(void *)> errf);

        /**
         * @brief 
         * 
         * @param func
         * @param args 
         * @param i 
         * @param errf 
         * @param pfa 
         * @param pfb 
         * @return int 
         */
        int simpInitInterval(
            std::function<double(double, void *)> func,
            void *args,
            mesh::Interval5 &i,
            std::function<int(void *)> errf,
            double *pfa,
            double *pfb);

        /**
         * @brief 
         * 
         * @param func
         * @param args 
         * @param i 
         * @param errf 
         * @return int 
         */
        int simpProcessInterval(
            std::function<double(double, void *)> func,
            void *args,
            mesh::Interval5 &i,
            std::function<int(void *)> errf);

        /**
         * @brief 
         * 
         * @param func
         * @param args 
         * @param i0 
         * @param i1 
         * @param i2 
         * @param errf 
         * @return int 
         */
        int simpSplitInterval(
            std::function<double(double, void *)> func,
            void *args,
            mesh::Interval5 &i0,
            mesh::Interval5 &i1,
            mesh::Interval5 &i2,
            std::function<int(void *)> errf);

        /**
         * @brief 
         * 
         * @param func
         * @param args 
         * @param i 
         * @param errf 
         * @param pfa 
         * @param pfb 
         * @return int 
         */
        int trapNLInitInterval(
            std::function<double(double, void *)> func,
            void *args,
            mesh::Interval5 &i,
            std::function<int(void *)> errf,
            double *pfa,
            double *pfb);

        /**
         * @brief 
         * 
         * @param func
         * @param args 
         * @param i 
         * @param errf 
         * @return int 
         */
        int trapNLProcessInterval(
            std::function<double(double, void *)> func,
            void *args,
            mesh::Interval5 &i,
            std::function<int(void *)> errf);

        /**
         * @brief 
         * 
         * @param func
         * @param args 
         * @param i0 
         * @param i1 
         * @param i2 
         * @param errf 
         * @return int 
         */
        int trapNLSplitInterval(
            std::function<double(double, void *)> func,
            void *args,
            mesh::Interval5 &i0,
            mesh::Interval5 &i1,
            mesh::Interval5 &i2,
            std::function<int(void *)> errf);
        /**
         * @brief 
         * 
         * @param func
         * @param args 
         * @param i 
         * @param errf 
         * @param pfa 
         * @param pfb 
         * @return int 
         */
        int cadreInitInterval(
            std::function<double(double, void *)> func,
            void *args,
            mesh::Interval9 &i,
            std::function<int(void *)> errf,
            double *pfa,
            double *pfb);

        /**
         * @brief 
         * 
         * @param func
         * @param args 
         * @param i 
         * @param errf 
         * @return int 
         */
        int cadreProcessInterval(
            std::function<double(double, void *)> func,
            void *args,
            mesh::Interval9 &i,
            std::function<int(void *)> errf);

        /**
         * @brief 
         * 
         * @param func
         * @param args 
         * @param i0 
         * @param i1 
         * @param i2 
         * @param errf 
         * @return int 
         */
        int cadreSplitInterval(
            std::function<double(double, void *)> func,
            void *args,
            mesh::Interval9 &i0,
            mesh::Interval9 &i1,
            mesh::Interval9 &i2,
            std::function<int(void *)> errf);

        /**
         * @brief 
         * 
         * @param func
         * @param args 
         * @param i 
         * @param errf 
         * @return int 
         */
        int gk49ProcessInterval(
            std::function<double(double, void *)> func,
            void *args,
            mesh::Interval &i,
            std::function<int(void *)> errf);

        /**
         * @brief 
         * 
         * @param func
         * @param args 
         * @param i0 
         * @param i1 
         * @param i2 
         * @param errf 
         * @return int 
         */
        int gk49SplitInterval(
            std::function<double(double, void *)> func,
            void *args,
            mesh::Interval &i0,
            mesh::Interval &i1,
            mesh::Interval &i2,
            std::function<int(void *)> errf);
        /**
         * @brief 
         * 
         * @param func
         * @param args 
         * @param i 
         * @param errf 
         * @return int 
         */
        int gk715ProcessInterval(
            std::function<double(double, void *)> func,
            void *args,
            mesh::Interval &i,
            std::function<int(void *)> errf);

        /**
         * @brief 
         * 
         * @param func
         * @param args 
         * @param i0 
         * @param i1 
         * @param i2 
         * @param errf 
         * @return int 
         */
        int gk715SplitInterval(
            std::function<double(double, void *)> func,
            void *args,
            mesh::Interval &i0,
            mesh::Interval &i1,
            mesh::Interval &i2,
            std::function<int(void *)> errf);
        /**
         * @brief 
         * 
         * @param func
         * @param args 
         * @param i 
         * @param errf 
         * @return int 
         */
        int gk1021ProcessInterval(
            std::function<double(double, void *)> func,
            void *args,
            mesh::Interval &i,
            std::function<int(void *)> errf);

        /**
         * @brief 
         * 
         * @param func
         * @param args 
         * @param i0 
         * @param i1 
         * @param i2 
         * @param errf 
         * @return int 
         */
        int gk1021SplitInterval(
            std::function<double(double, void *)> func,
            void *args,
            mesh::Interval &i0,
            mesh::Interval &i1,
            mesh::Interval &i2,
            std::function<int(void *)> errf);

        /**
         * @brief 
         * 
         * @param func
         * @param args 
         * @param errf 
         * @param c 
         * @param z0 
         * @param xg 
         * @param xk 
         * @param wg 
         * @param wgk 
         * @param ng 
         * @param I 
         * @param err 
         * @return int 
         */
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
            double *err);
    } // namespace integrate 
    
} // namespace afterglowpy

#endif
