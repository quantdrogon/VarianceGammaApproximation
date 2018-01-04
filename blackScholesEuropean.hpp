//
//  blackScholesEuropean.hpp
//  Whaley
//
//  Created by Abhishek Sanghani on 10/2/17.
//  Copyright © 2017 Abhishek Sanghani. All rights reserved.
//

#ifndef blackScholesEuropean_hpp
#define blackScholesEuropean_hpp

#include "importantHeaders.h"

class BlackScholesEuropeanOption{
    private:
    double premium, delta, gamma, vega, theta;
   
    
    
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
