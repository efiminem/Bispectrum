#include "ultraspherical.h"

int main() {
ULTRAS vec(4); //dimensions
std::vector<double> angles = {0.5, 0.5, 0.5};
std::vector<int> indices = {1,2,3};
vec.set(angles);
vec.calculate(100);
std::cout << vec.get(indices).r << '\n';
//std::cout << vec.multiplier[1].r << "*" << vec.multipliers[0][PT(2,1)]; << "*" << vec.multipliers[1][PT(0,0)]; //TEST
return 0;
}
