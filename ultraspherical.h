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

    double GetNullValue (int);
    std::vector<double> CalculateAnotherMultiplier (int, const std::size_t, double);
    std::vector<complex> CalculateFirstMultiplier(const std::size_t, double);

public:
    complex usphFunction;
    std::vector<complex> multiplier; //j=1 complex number
    std::vector<std::vector<double>> multipliers; //j!=1 real numbers
    void calculate(std::size_t);

    template<typename T = double> complex Y(int);
    template<typename T = double, typename... Args> complex Y(int, Args...);
    template<typename T = double> void ULTRAS_REST(T);
    template<typename T = double, typename... Args> void ULTRAS_REST(T, Args...);
    template<typename T = double> ULTRAS(T);
    template<typename T = double, typename... Args> ULTRAS(T, Args...);
};
