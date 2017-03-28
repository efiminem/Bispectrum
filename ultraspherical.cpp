#include <cstddef>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>

int PT(int l, int m) {return m+(l*(l+1)/2);}
int YR(int l, int m) {return m+l+(l*l);}

struct complex {
double r, i;
};

void printException (std::string message) {
    std::cout << message;
}

class Hardcode {
private:
    std::vector<double> A_v = std::vector<double>(100,0);
    std::vector<double> B_v = std::vector<double>(100,0);
    void expand (std::size_t maxL) {
        if (PT(maxL, maxL) >= 100) {
            for (int i=0; i<PT(maxL, maxL)- 100; i++) {
                A_v.push_back(0);
                B_v.push_back(0);
            }
        }
    }
public:
    double SQRT2 = 1.4142135623730950488;
    double SQRT3 = 1.7320508075688772935;
    double SQRT3DIV2 = -1.2247448713915890491;
    std::vector<double> nullValue;
    void update (std::size_t maxL, int j) {
        nullValue.push_back(0.282094791773878143474); // = 
sqrt(1/(4*PI))
        nullValue.push_back(0.797884560802865355879892); // = 
sqrt(2/(PI))
        expand(maxL);
        for (std::size_t l=2; l <= maxL ; l++) {
            double ls = l*l, lm1s = (l-1)*(l-1);
            for (std::size_t m=0; m <l-1; m++) {
                double ms=m*m;
                A_v[PT(l,m)]=sqrt((2*l+j-3)*(2*l + j - 
1)/((l-m)*(l+m+j-2)));
                
B_v[PT(l,m)]=-sqrt((l+m+j-3)*(l-m-1)/((2*l+j-3)*(2*l+j-5)));
            }
        }
    }
    double A(int l, int m, int j) {return A_v[PT(l,m)];}
    double B(int l, int m, int j) {return B_v[PT(l,m)];}
};

class ULTRAS {
private:
    std::size_t numofArguments;
    Hardcode tab; //table of already calculated factors
    std::vector<double> P_v; //vector of ASF
    std::vector<double> Y_v_r; //real spherical part
    std::vector<double> Y_v_i; //imaginary spherical part
    std::vector<double> H_v_r; //real spherical part
    std::vector<double> H_v_i; //imaginary spherical part
    std::vector<double> Theta; // vector of angles
    void computeALF (std::size_t L, double x, double j) {
        const double sintheta = sqrt(1.-x*x);
        double temp = tab.nullValue[j-2];
        /*if (PT(L, L) >= 100) {
            for (int i=0; i<PT(L, L)- 100; i++) {
                P_v.push_back(0);
                Y_v_r.push_back(0);
            }
        }*/
        P_v[PT(0,0)] = temp;
        if (L > 0) {
            P_v[PT(1,0)] = x*sqrt(j+1)*temp; //TODO: get it from 
Hardcode
            temp = -sqrt((j+1)/j) * sintheta * temp;
            P_v[PT(1,1)] = temp;
            for (std::size_t l=2; l <= L ; l++) {
                for (std::size_t m =0; m <l-1; m++) {
                P_v[PT(l,m)]=tab.A(l, m, j)*(x*P_v[PT(l-1,m)]+tab.B(l,m, 
j)* P_v[PT(l-2,m)]);
                }
            P_v[PT(l,l-1)] = x*sqrt(2*(l-1)+j+1)*temp;
            temp=-sqrt((2*l + j + 1)/(2*l + j - 2))*sintheta*temp;
            P_v[PT(l,l)]=temp;
            }
        }
    }
    void computeY (const std::size_t L, double y) {
        for (std::size_t l=0; l <= L; l++) Y_v_r[YR(l,0)] = 
P_v[PT(l,0)];
        double c1 = 1.0,c2 = y;
        double s1 = 0.0,s2 = -sqrt(1.-y*y);
        double tc = 2.0*c2;
        for (std::size_t m =1; m <= L; m++) {
            double s = tc*s1-s2;
            double c = tc*c1-c2;
            s2 = s1; s1 = s; c2 = c1; c1 = c;
            for (std::size_t l = m; l <= L; l++) {
                Y_v_r[YR(l,m)]=P_v[PT(l,m)]*c;
                Y_v_r[YR(l,-m)]=-Y_v_r[YR(l,m)];
                Y_v_i[YR(l,m)]=-P_v[PT(l,m)]*s;
                Y_v_i[YR(l,-m)]=Y_v_i[YR(l,m)];
            }
        }
    }
    void computeH (const std::size_t L, double y) {
        for (std::size_t l=0; l <= L; l++) H_v_r[YR(l,0)] = 
P_v[PT(l,0)];
        double c1 = 1.0,c2 = y;
        double s1 = 0.0,s2 = -sqrt(1.-y*y);
        double tc = 2.0*c2;
        for (std::size_t m =1; m <= L; m++) {
            double s = tc*s1-s2;
            double c = tc*c1-c2;
            s2 = s1; s1 = s; c2 = c1; c1 = c;
            for (std::size_t l = m; l <= L; l++) {
                H_v_r[YR(l,m)]=P_v[PT(l,m)]*c;
                H_v_r[YR(l,-m)]=-H_v_r[YR(l,m)];
                H_v_i[YR(l,m)]=-P_v[PT(l,m)]*s;
                H_v_i[YR(l,-m)]=H_v_i[YR(l,m)];
            }
        }
    }
public:
    complex usphFunction;
    void calculate(std::size_t maxL) {
    P_v = std::vector<double>(PT(maxL, maxL)+1,0);
    Y_v_r = std::vector<double>(YR(maxL, maxL)+1,0);
    Y_v_i = std::vector<double>(YR(maxL, maxL)+1,0);
    H_v_r = std::vector<double>(YR(maxL, maxL)+1,0);
    H_v_i = std::vector<double>(YR(maxL, maxL)+1,0);
    tab.update(maxL, 2);
    computeALF(maxL, Theta[1], 2);
    computeY(maxL, Theta[0]);
    if (numofArguments > 1) {tab.update(maxL, 3); computeALF(maxL, 
Theta[2], 3);computeH(maxL, Theta[0]);}
    }
    template<typename T = double, typename... Args> complex Y(int m, int 
l, int k, Args... args) {
        usphFunction.r = Y_v_r[YR(l,m)];
        std::cout << Y_v_r[YR(k,l)] << "|";
        usphFunction.i = Y_v_i[YR(l,m)];
        usphFunction.r *= P_v[YR(k,l)];
        std::cout << H_v_r[YR(k,l)] << "|";
        usphFunction.i *= H_v_i[YR(k,l)];
        return usphFunction;
     }

    template<typename T = double> void ULTRAS_REST(T t) {
    Theta.push_back(t);
    //end of loop

    }
    template<typename T = double, typename... Args> void ULTRAS_REST(T 
t, Args... args) {
    Theta.push_back(t);
    ULTRAS_REST(args...);
    }
    template<typename T = double> ULTRAS(T t) {
    printException("Error: you need two or more variables");
    }
    template<typename T = double, typename... Args> ULTRAS(T t, Args... 
args) {
    Theta.push_back(t);
    numofArguments = sizeof...(args);
    ULTRAS_REST(args...);
    }

};


int main() {
ULTRAS vec(0.707, 0.707, 0.707);
//TODO: what's wrong with 10000 ?
vec.calculate(100);
std::cout << vec.Y(1, 1, 1).r;
return 0;
}

