//
//  whaley.hpp
//  Whaley
//
//  Created by Abhishek Sanghani on 10/6/17.
//  Copyright © 2017 Abhishek Sanghani. All rights reserved.
//

#ifndef whaley_hpp
#define whaley_hpp

#include <stdio.h>
#include "importantHeaders.h"
#include "normal_distribution.hpp"

class whaleyOption{
private:
    double S, K, r, q, sig, T, europeanPremium, deltaE, gammaE, vegaE ;
    double AmerPrice, AmerDelta, AmerGamma, AmerVega;
    double Sx;
    double tol;
    std::string PutCall;
    int maxIterNewton;
    //boost::math::normal Z;
    // derived numbers for calculations
    double denomD1, discountedRate, discountedDividend, c2;
    double n, k;
    //N = 2b/sig^2 , here n==N
    //M = 2r/sig^2, k = M/(1-exp(-rT)
    double dNdSig, dKdSig, tmp, q1,  dQ1dSig, q2, dQ2dSig;
    
    //void findGSxWhaley(double Sx, double K, double b1, double b2, double c2, double denomD1, double T, double qI, double phi);
    
    std::tuple<double, double> findGSxWhaley(double Sx, double phi, double qI);
    
    double findSxViaNewton(double Sx, double qI, double phi);
    
public:
    whaleyOption(double S, double K, double r, double q, double sig, double T, double europeanPremium, double deltaE, double gammaE, double vegaE, std::string PutCall, int maxIterNewton, double tol);
    void whaleyPricer();
    double getPremium();
    double getDelta();
    double getGamma();
    double getVega();
    //void findGSxWhaley(double Sx,);
    
};

#endif /* whaley_hpp */
