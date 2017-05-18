#include <cstddef>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280

int PT(int l, int m) {return m+(l*(l+1)/2);}

struct complex {
double r, i;
};

class ULTRAS {
private:
    std::size_t numofArguments;
    std::size_t numofIndices;
    std::vector<double> Theta;
    std::vector<int> Index;

    double GetNullValue (int j) {
    
    }
    std::vector<double> CalculateAnotherMultiplier (int j, const std::size_t maxL, double sintheta) {
        std::vector<double> multiplier;
        double costheta = sqrt(1.-sintheta*sintheta);
        double temp=0;
        for (int h=0; h < PT(maxL,maxL); h++) multiplier.push_back(0);
        if (j==2) temp= 1/sqrt(2);
        if (j==3) temp=sqrt(2/PI); //TODO Null Value
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
            temp=-sqrt((2*l + j + 1.)/(2*l + j - 2.))*sintheta*temp;
            multiplier[PT(l,l)]=temp;
            }
         }
         return multiplier;
    }
    std::vector<complex> CalculateFirstMultiplier(const std::size_t maxL, double sintheta) {
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

public:
    complex usphFunction;
    std::vector<complex> multiplier;
    std::vector<std::vector<double>> multipliers;
    void calculate(std::size_t maxL) {
       multiplier = CalculateFirstMultiplier(maxL, Theta[0]); //j=1 complex number
       for (int j=2; j<=numofArguments; j++) {
          multipliers.push_back(CalculateAnotherMultiplier(j, maxL, Theta[j-1])); //j!=1 real numbers
       }

    }
    template<typename T = double> complex Y(int t) {
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
    template<typename T = double, typename... Args> complex Y(int t, Args... args) {
        Index.push_back(t);
        return Y(args...);
     }
    template<typename T = double> void ULTRAS_REST(T t) {
    Theta.push_back(t);
    //end of loop
    }
    template<typename T = double, typename... Args> void ULTRAS_REST(T t, Args... args) {
    Theta.push_back(t);
    ULTRAS_REST(args...);
    }
    template<typename T = double> ULTRAS(T t) {
    Theta.push_back(t);
    }
    template<typename T = double, typename... Args> ULTRAS(T t, Args... args) {
    Theta.push_back(t);
    numofArguments = sizeof...(args); // one lost
    numofArguments++;
    ULTRAS_REST(args...);
    }

};


int main() {
ULTRAS vec(0.707, 0.707, 0.707); //sintheta
vec.calculate(1000);
std::cout << vec.Y(0, 1, 2).r << '\n';
return 0;
}

