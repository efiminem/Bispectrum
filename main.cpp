#include "ultraspherical.h"

int main() {
ULTRAS vec(4); //dimensions
vec.set({0.5, 0.5, 0.5});
vec.pack(100);
std::cout << vec.get({1,2,3}).r << '\n';
return 0;
}
