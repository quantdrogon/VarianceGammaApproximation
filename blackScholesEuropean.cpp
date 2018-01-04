//
//  blackScholesEuropean.cpp
//  Whaley
//
//  Created by Abhishek Sanghani on 10/2/17.
//  Copyright © 2017 Abhishek Sanghani. All rights reserved.
//

#include "blackScholesEuropean.hpp"

void BlackScholesEuropeanOption::calculateValues()
{
    double top = log(s/K) + (r - q + pow(v,2)/2) * T;
    double bottom = v * sqrt(T);
    double d1 = top / bottom;
    double d2 = d1 - v*sqrt(T);
    double b1 = exp(-q*T);
    double b2 = exp(-r*T);
    
    gamma = b1*norm_cdf(d1)/(s*bottom);
    vega = s*b1*norm_pdf(d1)/100;
    
    if (optionType == "C")
    {
        double nd1 = norm_cdf(d1);
        double nd2 = norm_cdf(d2);
        premium = s*b1*nd1 - K*b2*nd2;
        delta = b1*nd1;
        theta = -(b1*s*norm_pdf(d1)*v/(2*bottom)) - r*K*b2*norm_cdf(d2)+q*s*b1*norm_cdf(d1);
    }
    else if(optionType == "P")
    {
        double nNd1 = norm_cdf(-d1);
        double nNd2 = norm_cdf(-d2);
        premium = K*b2*nNd2 - s*b1*nNd1;
        delta = -b1*nNd2;
        theta = -(b1*s*norm_pdf(d1)*v/(2*bottom)) + r*K*b2*norm_cdf(-d2) - q*s*b1*norm_cdf(-d1);
    }
    
    
    
}

BlackScholesEuropeanOption::BlackScholesEuropeanOption(double s, double K, double r, double q, double v, double T, std::string optionType): s(s), K(K), r(r), q(q), v(v), T(T), optionType(optionType)
{
    calculateValues();
}

double BlackScholesEuropeanOption::getPremium()
{
    return premium;
}

double BlackScholesEuropeanOption::getDelta()
{
 return delta;
}

double BlackScholesEuropeanOption::getGamma()
{
    return gamma;
}

double BlackScholesEuropeanOption::getVega()
{
    return vega;
}

double BlackScholesEuropeanOption::getTheta()
{
    return theta;
}
