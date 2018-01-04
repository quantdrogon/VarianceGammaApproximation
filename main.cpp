//
//Barone-Adesi and Whaley quadratic approximation to vanilla options
//  main.cpp
//  Whaley
//
//  Created by Abhishek Sanghani on 10/2/17.
//  Copyright © 2017 Abhishek Sanghani. All rights reserved.
//
//Edited by Amir Oskoui 10/30/2017


#include <iostream>
#include "BlackScholesEuropean.hpp"
#include "whaley.hpp"
using namespace std;

int main(int argc, const char * argv[]) {
    // Input the option parameters
    double S = 40;                    // Spot Price
    double r = 0.0488;                 // Risk free rate
    double q = 0.0;                    // Dividend yield
    
    /*
    double K = 110;                   // Strike Price
    double v = 0.2;                   // Volatility
    double T = 0.05;                   // Maturity
    */
    std::string PutCall = "P";              // 'C'all or 'P'ut
    
    //int sizeK;
    
    std::vector<double> K_vec = {35,35,35,40,40,40,45,45,45,35,35,35,40,40,40,45,45,45,35,35,35,40,40,40,45,45,45};
    std::vector<double> v_vec = {0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4};
    std::vector<double> T_vec = {0.0833,0.3333,0.5833,0.0833,0.3333,0.5833,0.0833,0.3333,0.5833,0.0833,0.3333,0.5833,0.0833,0.3333,0.5833,0.0833,0.3333,0.5833,0.0833,0.3333,0.5833,0.0833,0.3333,0.5833,0.0833,0.3333,0.5833};
    
    int maxIterNewton = 100;
    double tol = pow(10,-6);
    
    for (int i =0; i < K_vec.size(); i++)
    {
        BlackScholesEuropeanOption bs_opt(S, K_vec.at(i), r, q, v_vec.at(i), T_vec.at(i), PutCall);
        //whaleyOption w_opt(S, K_vec.at(i), r, q, v_vec.at(i), T_vec.at(i), bs_opt.getPremium(), bs_opt.getDelta(), bs_opt.getGamma(), bs_opt.getVega(), PutCall, maxIterNewton, tol);
        JuZhongOption jz_opt(S, K_vec.at(i), r, q, v_vec.at(i), T_vec.at(i), bs_opt.getPremium(), bs_opt.getDelta(), bs_opt.getDelta(), bs_opt.getGamma(), bs_opt.getVega(), bs_opt.getTheta(), PutCall)
        //std::cout<<w_opt.getPremium()<<endl;
        std::cout<<jz_opt.getPremium()<<endl;
    }
    
    return 0;
}
