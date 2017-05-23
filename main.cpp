#include "ultraspherical.h"

int main() {
ULTRAS vec(0.5, 0.5, 0.5, 0.5, 0.5); //sintheta
vec.calculate(100);
std::cout << vec.Y(1, 2, 3, 4, 5).r << '\n';
//std::cout << vec.multiplier[1].r << "*" << vec.multipliers[0][PT(2,1)]; << "*" << vec.multipliers[1][PT(0,0)]; //TEST
return 0;
}
