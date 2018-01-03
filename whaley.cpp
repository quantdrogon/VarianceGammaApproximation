//
//  whaley.cpp
//  Whaley
//
//  Created by Abhishek Sanghani on 10/6/17.
//  Copyright © 2017 Abhishek Sanghani. All rights reserved.
//

#include "whaley.hpp"

whaleyOption::whaleyOption(double S, double K, double r, double q, double sig, double T, double europeanPremium, double deltaE, double gammaE, double vegaE, std::string PutCall, int maxIterNewton, double tol) : S(S), K(K), r(r), q(q), sig(sig), T(T), europeanPremium(europeanPremium), deltaE(deltaE), gammaE(gammaE), vegaE(vegaE), PutCall(PutCall), maxIterNewton(maxIterNewton), tol(tol)
{
    // calculate the derived numbers
    discountedRate = exp(-r*T);
    discountedDividend = exp(-q*T);
    n = 2*(r-q)/pow(sig,2); //N = 2b/sig^2 , here n==N;
    k = 2*r/pow(sig,2)/(1-discountedRate); //M = 2r/sig^2, k = M/(1-exp(-rT);
    
    dNdSig = -2*n/sig;
    dKdSig = -2*k/sig;
    tmp = pow(n-1,2)+4*k;
    q1 = (1-n-sqrt(tmp))/2;
    dQ1dSig = 0.5*( -dNdSig - 0.5 *(2*(n-1)*dNdSig + 4*dKdSig) / sqrt(tmp));
    
    q2 = (1-n+sqrt(tmp))/2;
    dQ2dSig = 0.5*( -dNdSig + 0.5 *(2*(n-1)*dNdSig + 4*dKdSig) / sqrt(tmp));
    
    c2 = r-q+pow(sig,2)/2;
    denomD1 = sig*sqrt(T);
    whaleyPricer();
    
}

void whaleyOption::whaleyPricer()
{
    bool europeanStyle = false;
    
    double sNew;
    
    // Quadratic approximation
    if(PutCall == "C")
    {
        double phi = 1;
        
        double Sx1 = K;
        
        double gSx1, d1gSx1;
        std::tie(gSx1, d1gSx1) = findGSxWhaley(Sx1, phi, q2);
        
        double Sx2;
        if (S <= K)
            {Sx2 = 5.0*K;}
        else
            {Sx2 = 5.0*S;}
        double gSx2, d1gSx2;
        std::tie(gSx2, d1gSx2) = findGSxWhaley(Sx2, phi, q2);
        
        if ( (std::isnan(gSx1) == true) or (std::isnan(gSx2) == true) or ((gSx1*gSx2)>=0))
        {
            europeanStyle = true;
        }
        
        if (europeanStyle == true)
        {
            AmerPrice = europeanPremium;
            AmerVega = vegaE;
            AmerGamma = gammaE;
            Sx = std::numeric_limits<double>::quiet_NaN();
        }
        else
        {
            
            if (fabs(d1gSx2) > fabs(d1gSx1))
            {
                sNew = Sx2;
            }
            else
            {
               sNew = Sx1;
            }
            
            double sPrevious = sNew;
            int counter = 0;
            bool flag = true;
            
            while (flag == true)
            {
                counter++;
                
                sNew = findSxViaNewton(sNew, q2, phi);
                
                if (fabs(sNew - sPrevious)<tol)
                {   //disp(['counter =', num2str(counter), ' Newton-Raphson Converged']);
                    flag = false;
                    break;
                }
                if ((counter > maxIterNewton) or (sNew < 0))
                {   //%disp(['counter =', num2str(counter), ' break & switch to Bisection Method']);
                    break;
                }
                sPrevious = sNew;
            }
            
            counter = 0;
            
            while(flag == true)
            {
                counter++;
                
                double Sx_m = (Sx1+Sx2)/2;
                double gSx_m;
                std::tie(gSx_m, std::ignore) = findGSxWhaley(Sx_m, q2, phi);
                double check1 = gSx_m*gSx1;
                
                if (check1 >0)
                {
                    Sx1 = Sx_m;
                    gSx1 = gSx_m;
                    
                }
                else if (check1 <0)
                {
                    Sx2 = Sx_m;
                    //%gSx2 = gSx_m;
                }
                
                if (fabs(Sx1-Sx2) < tol)
                {
                    flag = false;
                    sNew = Sx_m;
                    //%disp(['counter =', num2str(counter), ' Bisection Method Converged']);
                }
            }
            
            Sx = sNew;
            double a2 = (log(Sx/K) + c2*T)/denomD1;
            double A2 = Sx*(1-discountedDividend*norm_cdf(a2))/q2;
            
            if (S<Sx)
            {
                AmerPrice = europeanPremium + A2*pow(S/Sx,q2);
                
                double d1 = (log(Sx/K) + c2*T)/denomD1;
                double dD1dSig = sqrt(T) - d1/sig;
                double dA2dSig = Sx*(-exp(-q*T)*norm_pdf(d1)*dD1dSig*q2 - (1-discountedDividend*norm_cdf(d1))*dQ2dSig)/pow(q2,2);
                AmerVega = vegaE + (dA2dSig *pow(S/Sx,q2) + A2 *pow(S/Sx,q2) * log(S/Sx) * dQ2dSig);
                
                AmerDelta = deltaE + A2*pow(S/Sx,q2) * q2/S;
                AmerGamma = gammaE + A2*pow(S/Sx,q2) * q2*(q2-1)/pow(S,2);
            }
            
            else
            {
                AmerPrice = fmax(S - K,0);
                AmerDelta = 1;
                AmerVega = 0;
                AmerGamma = 0;
                
            }
            
        }
    }
    
    else if (PutCall == "P")
    {
        double phi = -1;
        
        double Sx1;
        
        if (S >= K)
        {
            Sx1 = K/10.0;
            
        }
        else if (S < K)
        {
            Sx1 = S/10.0;
        }
        
        double gSx1, d1gSx1;
        std::tie(gSx1, d1gSx1) = findGSxWhaley(Sx1, phi, q1);
        
        double Sx2 = K;
        double gSx2, d1gSx2;
        std::tie(gSx2, d1gSx2) = findGSxWhaley(Sx2, phi,q1);
        
       
        
        
        if ( (std::isnan(gSx1) == true) or (std::isnan(gSx2) == true) or ((gSx1*gSx2)>=0))
        {
            europeanStyle = true;
        }
        
        if (europeanStyle == true)
        {
            AmerPrice = europeanPremium;
            AmerVega = vegaE;
            AmerGamma = gammaE;
            Sx = std::numeric_limits<double>::quiet_NaN();
        }
        else
        {
            
            if (fabs(d1gSx2) > fabs(d1gSx1))
            {
                sNew = Sx2;
            }
            else
            {
                sNew = Sx1;
            }
            
            double sPrevious = sNew;
            
            int counter = 0;
            bool flag = true;
            
            while (flag == true)
            {
                counter++;
                
                sNew = findSxViaNewton(sNew, q1, phi);
                
                if (fabs(sNew - sPrevious)<tol)
                {   //disp(['counter =', num2str(counter), ' Newton-Raphson Converged']);
                    flag = false;
                    break;
                }
                if ((counter > maxIterNewton) or (sNew < 0))
                {   //%disp(['counter =', num2str(counter), ' break & switch to Bisection Method']);
                    break;
                }
                sPrevious = sNew;
            }
            
            counter = 0;
            
            while(flag == true)
            {
                counter++;
                
                double Sx_m = (Sx1+Sx2)/2;
                double gSx_m;
                std::tie(gSx_m, std::ignore) = findGSxWhaley(Sx_m, q1, phi);
                double check1 = gSx_m*gSx1;
                
                if (check1 >0)
                {
                    Sx1 = Sx_m;
                    gSx1 = gSx_m;
                    
                }
                else if (check1 <0)
                {
                    Sx2 = Sx_m;
                    //%gSx2 = gSx_m;
                }
                
                if (fabs(Sx1-Sx2) < tol)
                {
                    flag = false;
                    sNew = Sx_m;
                    //%disp(['counter =', num2str(counter), ' Bisection Method Converged']);
                }
            }
            
            Sx = sNew;
            
            double a1 = -(log(Sx/K) + c2*T)/denomD1;
            double A1 = -Sx*(1-discountedDividend*norm_cdf(a1))/q1;
            
            if (S>Sx)
            {
                AmerPrice = europeanPremium + A1*pow(S/Sx,q1);
                //%disp([AmerPrice europeanPremium A1 S Sx q1]);
                
                double d1 = (log(Sx/K) + c2*T)/denomD1;
                double dD1dSig = sqrt(T) - d1/sig;
                double dA1dSig = -Sx*(exp(-q*T)*norm_pdf(-d1)*dD1dSig*q1 - (1-discountedDividend*norm_cdf(d1))*dQ1dSig)/pow(q1,2);
                AmerVega = vegaE + (dA1dSig *pow(S/Sx,q1) + A1 *pow(S/Sx,q1) * log(S/Sx) * dQ1dSig);
                
                AmerDelta = deltaE + A1*pow(S/Sx,q1) * q1/S;
                AmerGamma = gammaE + A1*pow(S/Sx,q1) * q1*(q1-1)/pow(S,2);
            }
            
            else
            {
                AmerPrice = fmax(K - S,0);
                AmerDelta = -1;
                AmerVega = 0;
                AmerGamma = 0;
                
            }
            
        }
    }
        
}


double whaleyOption::findSxViaNewton(double Sx, double qI, double phi)
{
    
    double tmp2 = Sx / qI;
    double c1 = log(Sx/K);
    double d1 = (c1 + c2*T) / denomD1;
    double d2 = d1 - denomD1;
    
    double gSx;
    double d1gSx;
    
    if (phi == 1)
    {
        double nd1 = norm_cdf(d1);
        double nd2 = norm_cdf(d2);
        double premium = Sx*discountedDividend*nd1 - K*discountedRate*nd2;
        double delta = discountedDividend*nd1;
        
        double a1 = (c1 + c2*T)/denomD1;
        double da1 = 1.0/ (Sx * denomD1);
        
        double tmp1 = 1.0-discountedDividend*norm_cdf(a1);
        double tmp3 = norm_pdf(a1);
        
        gSx = premium + tmp1*tmp2 - Sx + K;
        d1gSx = delta + tmp1/qI - discountedDividend * tmp3 *da1 * tmp2 - 1;
        //d2gSx  = gamma - (b1/qI)*tmp3*(2.0*daI + Sx*(aI*(daI^2)-d2aI));
        
        //d1y = 2*gSx*d1gSx;
        //d2y = 2*(d1gSx^2)+2*gSx*d2gSx;
        
    }
    
    else if (phi == -1)
    {
        double nNd1 = norm_cdf(-d1);
        double nNd2 = norm_cdf(-d2);
        double premium = K*discountedRate*nNd2 - Sx*discountedDividend*nNd1;
        double delta = -discountedDividend*nNd2;
        
        double a1 = -(c1 + c2*T)/denomD1;
        double da1 = -1.0/ (Sx * denomD1);
        
        double tmp1 = 1.0-discountedDividend*norm_cdf(a1);
        double tmp3 = norm_pdf(a1);
        
        gSx    = premium - tmp1*tmp2 + Sx - K;
        d1gSx  = delta - tmp1/qI + discountedDividend * tmp3 *da1 * tmp2 + 1;
        //d2gSx  = gamma + (b1/qI)*tmp3*(2.0*daI + Sx*(d2aI - aI*(daI^2)));
        
        //d1y = 2*gSx*d1gSx;
        //d2y = 2*(d1gSx^2)+2*gSx*d2gSx;
        
    }
    
    //Sx = Sx - d1y/d2y;
    Sx = Sx - gSx/d1gSx;

    return Sx;
    
}

std::tuple<double, double> whaleyOption::findGSxWhaley(double Sx, double phi, double qI)
{
    double tmp2 = Sx / qI;
    double c1 = log(Sx/K);
    double d1 = (c1 + c2*T) / denomD1;
    double d2 = d1 - denomD1;
    
    double gSx;
    double d1gSx;
    if (phi == 1)
    {
        double nd1 = norm_cdf(d1);
        double nd2 = norm_cdf(d2);
        double premium = Sx*discountedDividend*nd1 - K*discountedRate*nd2;
        double delta = discountedDividend*nd1;
        
        double a1 = (c1 + c2*T)/denomD1;
        double da1 = 1.0/ (Sx * denomD1);
        
        double tmp1 = 1.0-discountedDividend*norm_cdf(a1);
        double tmp3 = norm_pdf(a1);
        
        gSx    = premium + tmp1*tmp2 - Sx + K;
        d1gSx  = delta + tmp1/qI - discountedDividend * tmp3 *da1 * tmp2 - 1;
    }
    else if (phi == -1)
    {
        double nNd1 = norm_cdf(-d1);
        double nNd2 = norm_cdf(-d2);
        double premium = K*discountedRate*nNd2 - Sx*discountedDividend*nNd1;
        double delta = -discountedDividend*nNd2;
        
        double a1 = -(c1 + c2*T)/denomD1;
        double  da1 = -1.0/ (Sx * denomD1);
        
        double tmp1 = 1.0-discountedDividend*norm_cdf(a1);
        double tmp3 = norm_pdf(a1);
        
        gSx    = premium - tmp1*tmp2 + Sx - K;
        d1gSx  = delta - tmp1/qI + discountedDividend * tmp3 *da1 * tmp2 + 1;
    }
    
    return std::make_tuple(gSx, d1gSx);
  
}

//double AmerPrice, AmerDelta, AmerGamma, AmerVega;

double whaleyOption::getPremium()
{
    return AmerPrice;
}

double whaleyOption::getDelta()
{
    return AmerDelta;
}

double whaleyOption::getGamma()
{
    return AmerGamma;
}

double whaleyOption::getVega()
{
    return AmerVega;
}


