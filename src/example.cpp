#include <iostream>

#include <Eigen/Core>

#include "LeopardiPartition.hpp"
#include "lebedev/Lebedev.hpp"

int main() {
	auto [leo_weight, leo_sph_points] = leopardi_partition(20);
	std::cout << "Leopardi weight = " << leo_weight << std::endl;
	std::cout << "Leopardi points =\n" << leo_sph_points << std::endl;

	auto [leb_weight, leb_car_points] = lebedev(20);
	std::cout << "Lebedev weights =\n" << leb_weight << std::endl;
	std::cout << "Lebedev points =\n" << leb_car_points << std::endl;

}
