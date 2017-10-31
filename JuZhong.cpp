#include "JuZhong.hpp"


double JuZhongOption::findSxViaNewton(double initialGuess, double K, double r, double q, double sig, double T, double europeanPremium, double lambdaH, std::string PutCall, double maxIterNewton, double tol)
{
    int counter=0;
    bool flag = true;
    while (flag == true) and (counter<maxIterNewton)
    {
        counter++;
        if(PutCall == "C")
        {
            double phi = 1;
    
        }
        else {
            double phi =-1;
        }
        top = log(initialGuess/K) + (r - q + pow(sig,2)/2) * T;
        bottom = sig * sqrt(T);
        d1 = top / bottom;
        d2 = d1 - sig*sqrt(T);
        b1 = exp(-q*T);
        b2 = exp(-r*T);
        lhs = phi*initialGuess - lambdaH*(phi*(initialGuess-K));
        rhs = phi*initialGuess*b1*cdf(Z, phi*d1) - europeanPremium*lambdaH
  //Premium is the EuropeanPrice V_E
        if (abs((rhs-lhs/K)) <tol)
        {
            double Sx=initialGuess;
            flag = false;
            break;
        }

       else
       {
        double SlopeBi = b1*pdf(Z,phi*d1) / (sig*sqrt(T)) +( 1-lambdaH) * delta;
        double initialGuess=(lambdaH*K*phi+initialGuess*slopeBi -rhs)/(slopeBi-phi*(1-lambdaH));
        }
    }
    return Sx;
}


//void JuZhongOption::JuZhongPricer( double S, double K, double r, double q, double sig, double T, std::string PutCall) : S(S), K(K), r(r), q(q), sig(sig), T(T), PutCall(PutCall)
void JuZhongOption::JuZhongPricer()
{
     /*   alpha = 2 * r/(sig^2);
        beta = 2 * (r-q)/(pow(sig,2);
        hTau = 1 - exp(-r*T);
        lambdaH = (-(beta-1) + phi * sqrt(pow(beta-1,2) + 4 * alpha/hTau))/2;
    
        qInfty =  (1 - beta + phi * sqrt(pow(beta-1,2) + 4*alpha))/2;
        sInfty = K/(1 - 1/qInfty);
        hi = (-phi*(r-q)*T - 2*sig*sqrt(T)) * K / (phi * (sInfty - K));
        initialGuess = sInfty + (K - sInfty) * exp(hi);

        Sx = findSxViaJuZhong(initialGuess, K, r, q, sig, T, phi, lambdaH);
        ah = (phi * (Sx - K) - europeanPremium;
        //ah = (phi * (Sx - K) - BlackScholesEuropeanOption(Sx, K, r, q, sig, T, europeanPremium, lambdaH, PutCall)
        
    //NOTE: thetaE = BlackScholesEuropeanOption(Sx, K, r, q, sig, T, phi=phi, greekCal=True)[-1]
    //check if we can do this for theta
        lambdaHDerivation = -phi * alpha / (hTau^2 * sqrt(pow(beta-1,2) + 4*alpha/hTau));
        b = (1 - hTau) * alpha * lambdaHDerivation/(2*(2 * lambdaH + beta - 1));
        c = - (1 - hTau) * alpha / (2 * lambdaH + beta - 1) * (-thetaE/(hTau * ah * r * exp(-r*T)) + 1/hTau + lambdaHDerivation/(2 * lambdaH + beta - 1));
    // c requires theta from above so we need to check
        //euroPrice = BlackScholesEuropeanOption(S, K, r, q, sig, T, phi=phi);
        */
    if ((phi*(Sx -S))> 0)
    {
        AmerPrice = europeanPremium + (hTau * ah * pow((S/Sx),lambdaH))/(1 - b * (pow(log(S/Sx),2) - c * log(S/Sx));
    }

    else {
        AmerPrice = phi * (S - K);
    
    }
}
 double JuZhongOption::JuZhongOption(double S, double K, double r, double q, double sig, double T, double europeanPremium, double deltaE, double gammaE, double vegaE, double thetaE, std::string PutCall, int maxIterNewton, double tol): : S(S), K(K), r(r), q(q), sig(sig), T(T), europeanPremium(euroPremium), deltaE(deltaE), gammaE(gammaE), vegaE(vegaE), thetaE(thetaE), lambdaH(lambdaH), PutCall(PutCall), maxIterNewton(maxIterNewton), tol(tol);
{
        alpha = 2 * r/(sig^2);
        beta = 2 * (r-q)/(pow(sig,2);
        hTau = 1 - exp(-r*T);
        lambdaH = (-(beta-1) + phi * sqrt(pow(beta-1,2) + 4 * alpha/hTau))/2;
    
        qInfty =  (1 - beta + phi * sqrt(pow(beta-1,2) + 4*alpha))/2;
        sInfty = K/(1 - 1/qInfty);
        hi = (-phi*(r-q)*T - 2*sig*sqrt(T)) * K / (phi * (sInfty - K));
        initialGuess = sInfty + (K - sInfty) * exp(hi);

        Sx = findSxViaJuZhong(initialGuess, K, r, q, sig, T, phi, lambdaH);
        ah = (phi * (Sx - K) - europeanPremium)/hTau;
        //ah = (phi * (Sx - K) - BlackScholesEuropeanOption(Sx, K, r, q, sig, T, europeanPremium, lambdaH, PutCall)
        
    //NOTE: thetaE = BlackScholesEuropeanOption(Sx, K, r, q, sig, T, phi=phi, greekCal=True)[-1]
    //check if we can do this for theta
        lambdaHDerivation = -phi * alpha / (hTau^2 * sqrt(pow(beta-1,2) + 4*alpha/hTau));
        b = (1 - hTau) * alpha * lambdaHDerivation/(2*(2 * lambdaH + beta - 1));
        c = - (1 - hTau) * alpha / (2 * lambdaH + beta - 1) * (-thetaE/(hTau * ah * r * exp(-r*T)) + 1/hTau + lambdaHDerivation/(2 * lambdaH + beta - 1));
        JuZhongPricer();
}   
double JuZhongOption::getPremium()
{
    return AmerPrice;
}
