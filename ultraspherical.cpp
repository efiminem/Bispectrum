#include "ultraspherical.h"

double ULTRAS::GetNullValue (int j) {
    if (j==2) {return 1/sqrt(2);} // precomputed
    else if (j==3) {return sqrt(2/PI);}
    else if (j==4) {return sqrt(3)/2;}
    else {return sqrt(j)/(sqrt(2*PI)*GetNullValue(j-1));} // if not precomputed -> recursion
    }

std::vector<double> ULTRAS::CalculateAnotherMultiplier (int j, const std::size_t maxL, double sintheta) {
        std::vector<double> multiplier;
        double costheta = sqrt(1.-sintheta*sintheta);
        double temp=0;
        for (int h=0; h <= PT(maxL,maxL); h++) multiplier.push_back(0);
        temp = GetNullValue(j);
        multiplier[PT(0,0)] = temp;
        if (maxL > 0) {
            multiplier[PT(1,0)] = costheta*sqrt(j+1)*temp; //TODO: Get it from Hardcode
            temp = -sqrt((j+1.)/j) * sintheta * temp;
            multiplier[PT(1,1)] = temp;
            for (std::size_t l=2; l <= maxL ; l++) {
                for (std::size_t m =0; m <l-1; m++) {
                double A = sqrt((2*l+j-3.)*(2*l + j -1.)/((l-m)*(l+m+j-2.)));
                double B = -sqrt((l+m+j-3.)*(l-m-1)/((2*l+j-3.)*(2*l+j-5.)));
                multiplier[PT(l,m)]=A*(costheta*multiplier[PT(l-1,m)]+B* multiplier[PT(l-2,m)]);
                }
            multiplier[PT(l,l-1)] = costheta*sqrt(2*(l-1)+j+1)*temp;
            temp=-sqrt((2*l + j - 1.)/(2*l + j - 2.))*sintheta*temp;
            multiplier[PT(l,l)]=temp;
            }
         }
         return multiplier;
}

std::vector<complex> ULTRAS::CalculateFirstMultiplier(const std::size_t maxL, double sintheta) {
        std::vector<complex> multiplier;
        complex normalization_factor = {1/sqrt(2*PI), 1/sqrt(2*PI)};
        for (std::size_t l=0; l<= maxL; l++) multiplier.push_back(normalization_factor);
        double c1 = 1.0,c2 = sqrt(1.-sintheta*sintheta);
        double s1 = 0.0,s2 = -sintheta;
        double tc = 2.0*c2;
        for (std::size_t m =1; m <= maxL; m++) {
            double s = tc*s1-s2;
            double c = tc*c1-c2;
            s2 = s1; s1 = s; c2 = c1; c1 = c;
            multiplier[m].r*=c;
            multiplier[m].i*=s;
        }
        return multiplier;
}

void ULTRAS::calculate(std::size_t maxL) {
       multiplier = CalculateFirstMultiplier(maxL, Theta[0]);
       for (int j=2; j<=numofArguments; j++) {
          multipliers.push_back(CalculateAnotherMultiplier(j, maxL, Theta[j-1]));
       }

}
template<typename T = double> ULTRAS::complex Y(int t) {
       Index.push_back(t);
       numofIndices = Index.size();
       usphFunction = multiplier[Index[0]];
       if (numofArguments == numofIndices) {
          for (int h=2; h<=numofArguments; h++) {
             usphFunction.r*=multipliers[h-2][PT(Index[h-1], Index[h-2])];
             usphFunction.i*=multipliers[h-2][PT(Index[h-1], Index[h-2])];
          }
       }
       return usphFunction;
}
template<typename T = double, typename... Args> ULTRAS::complex Y(int t, Args... args) {
        Index.push_back(t);
        return Y(args...);
}
template<typename T = double> void ULTRAS::ULTRAS_REST(T t) {
    Theta.push_back(t);
    //end of loop
}
template<typename T = double, typename... Args> void ULTRAS::ULTRAS_REST(T t, Args... args) {
    Theta.push_back(t);
    ULTRAS_REST(args...);
}
template<typename T = double> ULTRAS::ULTRAS(T t) {
    Theta.push_back(t);
}
template<typename T = double, typename... Args> ULTRAS::ULTRAS(T t, Args... args) {
    Theta.push_back(t);
    numofArguments = sizeof...(args); // one lost
    numofArguments++;
    ULTRAS_REST(args...);
}
