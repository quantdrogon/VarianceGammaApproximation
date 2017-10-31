#ifndef JuZhong_hpp
#define JuZhong_hpp

#include <stdio.h>
#include "importantHeaders.h"
class JuZhongOption{
private:
    double S, K, r, q, sig, T, europeanPremium, thetaE, deltaE, gammaE, vegaE ;
    double AmerPrice, AmerDelta, AmerGamma, AmerVega, AmerTheta;
    double alpha, beta, hTau, lambdaH, qInfty, SInfty, hi, ah, theta, b, c, lambdaHDerivation;
    double theta;
    double initialGuess;
    double tol;
    std::string PutCall;
    int maxIterNewton;
    boost::math::normal Z;
    // derived numbers for calculations
      // derived numbers for calculations
    double denomD1, discountedRate, discountedDividend, c2;
    double n, k;
    //N = 2b/sig^2 , here n==N
    //M = 2r/sig^2, k = M/(1-exp(-rT)
    double dNdSig, dKdSig, tmp, q1,  dQ1dSig, q2, dQ2dSig;
    double leftside;
    double rightside;
    
    //void findGSxWhaley(double Sx, double K, double b1, double b2, double c2, double denomD1, double T, double qI, double phi);
    
    //std::tuple<double, double> findGSxWhaley(double Sx, double phi, double qI);
    
    double findSxViaNewton(double initialGuess, double K, double r, double q, double sig, double T, double Sx, double qI, std::string PutCall);
    double Sx, top, bottom, d1, d2 b1, b2, lhs, rhs;
public:
    JuZhongOption(double S, double K, double r, double q, double sig, double T, double europeanPremium, double deltaE, double gammaE, double vegaE, double thetaE, std::string PutCall, int maxIterNewton, double tol);
    void JuZhongPricer();
    double getPremium();
    double getDelta();
    double getGamma();
    double getVega();
    double getTheta();
    //void findGSxJuZhong(double Sx,);
    
};

#endif /* JuZhong_hpp */