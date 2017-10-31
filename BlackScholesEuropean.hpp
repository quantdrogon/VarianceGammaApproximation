//
//  blackScholesEuropean.hpp
//  Whaley
//
//  Created by Abhishek Sanghani on 10/2/17.
//  Copyright © 2017 Abhishek Sanghani. All rights reserved.
// Amir Oskoui edit 10/30/2017 theta input

#ifndef blackScholesEuropean_hpp
#define blackScholesEuropean_hpp

#include "importantHeaders.h"

class BlackScholesEuropeanOption{
    private:
    double premium, delta, gamma, vega, theta;
    boost::math::normal Z;
    
    
    public :
    
    //s, K, r, q, v, T, optionType
    BlackScholesEuropeanOption(double s, double K, double r, double q, double v, double T, std::string optionType);
    double s, K, r, q, v, T;
    std::string optionType;
    void calculateValues();
    double getPremium();
    double getDelta();
    double getGamma();
    double getVega();
    double getTheta();
};


#endif /* blackScholesEuropean_hpp */