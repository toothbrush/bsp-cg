
#include "Matrix2HilbertCoordinates.hpp"
#include <iostream>
#include <sstream>

int main() {
	unsigned long int i;
	std::cout << "i: " << std::endl;
	std::cin >> i;
	unsigned long int j;
	std::cout << "j: " << std::endl;
	std::cin >> j;
	
	unsigned long int h1, h2;
	Matrix2HilbertCoordinates::IntegerToHilbert( i, j, h1, h2 );
	std::cout << "Most significant bits: \t" << h1 << ", ";
	std::cout << "least significant bits: \t" << h2 << std::endl;
}

