#include <cstddef>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280

struct complex {
    double r, i;
};

class ULTRAS {
private:
    std::size_t numofArguments;
    std::size_t numofIndices;
    std::vector<double> Theta;
    std::vector<int> Index;
    int PT(int, int);

    double GetNullValue (int);
    std::vector<double> CalculateAnotherMultiplier (int, const std::size_t, double);
    std::vector<complex> CalculateFirstMultiplier(const std::size_t, double);

public:
    complex usphFunction;
    std::vector<complex> multiplier; //j=1 complex number
    std::vector<std::vector<double>> multipliers; //j!=1 real numbers
    void calculate(std::size_t);
    complex get(std::vector<int>);
    void set(std::vector<double>);
    ULTRAS(int);
};
