#ifndef SHOCKEVOLUTION_HPP
#define SHOCKEVOLUTION_HPP

static const double T0_inj = 1.0e3;

namespace afterglowpy
{
    /**
     * @brief Shock velocity in lab frame (four-velocity) 
     * 
     * @param u 
     * @return double 
     */
    double shockVel(double u);

    /**
     * @brief The injected energy
     * 
     * @param te 
     * @param L0 
     * @param q 
     * @param ts 
     * @return double 
     */
    double E_inj(double te, double L0, double q, double ts);

    /**
     * @brief The injected luminosity
     * 
     * @param te 
     * @param L0 
     * @param q 
     * @param ts 
     * @return double 
     */
    double L_inj(double te, double L0, double q, double ts);

    /**
     * @brief Callback for shock deceleration effects to take place
     * 
     * @param t0 
     * @param R0 
     * @param u0 
     * @param argv 
     */
    void shockInitDecel(double t0, double * R0, double * u0, void *argv); 

    /**
     * @brief Find the location of the shock front
     * 
     * @param t0 
     * @param R0 
     * @param u0 
     * @param tRes 
     * @param argv 
     */
    void shockInitFind(double t0, double * R0, double * u0, double tRes, void *argv);

    /**
     * @brief 
     * 
     * @param t 
     * @param x 
     * @param argv 
     * @param xdot 
     */
    void Rudot2D(double t, double * x, void *argv, double * xdot);
    void RuThdot3D(double t, double * x, void *argv, double * xdot, int spread);

    /**
     * @brief Evolve shock wave using RK4 ode solver
     * 
     * @param t 
     * @param R 
     * @param u 
     * @param N 
     * @param R0 
     * @param u0 
     * @param args 
     */
    void shockEvolveRK4(
        std::vector<double> &t,
        std::vector<double> &R, 
        std::vector<double> &u, 
        int N, 
        double R0, 
        double u0, 
        void *args);

    /**
     * @brief Evolve the shock wave using RK4 solver, but w/ spreading effects
     * 
     * @param t 
     * @param R 
     * @param u 
     * @param th 
     * @param N 
     * @param R0 
     * @param u0 
     * @param th0 
     * @param args 
     * @param spread 
     */
    void shockEvolveSpreadRK4(
        std::vector<double> &t, 
        std::vector<double> &R, 
        std::vector<double> &u, 
        std::vector<double> &th, 
        int N, 
        double R0, 
        double u0, 
        double th0, 
        void *args,
        int spread);
    
} // namespace afterglowpy

#endif
