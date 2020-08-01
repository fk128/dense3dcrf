/*
 * -------------------- Class that holds the parameters of the CRF-----------------------
 */
#include<iostream>
#include <vector>

class Parameters {

public :
	int maxIterations;
	float posW;
	std::vector<float> posXYZStds;
	std::vector<float> bilateralXYZStds;
	std::vector<float> bilateralModsStds;
	int numberOfModalities;
	float bilateralW;

	Parameters() {
		numberOfModalities = 1;
		maxIterations = 2;
		posXYZStds = std::vector<float>(3, 3.0);
		posW = 3;
		std::vector<float> bilateralModsStds;
		bilateralXYZStds = std::vector<float>(3, 10.0);
		bilateralW = 5;
	}
	~Parameters() {

	}

	void setNumberOfModalities(int numbOfModalities) {
		numberOfModalities = numbOfModalities;
		bilateralModsStds = std::vector<float>(numbOfModalities, 3);
	}

	void print() {
		std::cout << "========== Parameters for the CRF: ============" << std::endl;
		std::cout << "maxIterations: " << maxIterations << std::endl;
		std::cout << "numberOfModalities: " << numberOfModalities << std::endl;
		for (int mod_i = 0 ; mod_i < numberOfModalities ; mod_i++) {
			std::cout << "bilateralModsStds[" << mod_i << "]: " << bilateralModsStds[mod_i] << std::endl;
		}
		std::cout << "posW: "     << posW    << std::endl;
		for (int mod_i = 0 ; mod_i < posXYZStds.size() ; mod_i++) {
			std::cout << "posXYZStds[" << mod_i << "]: " << posXYZStds[mod_i] << std::endl;
		}

		std::cout << "bilateralW: "     << bilateralW    << std::endl;



		for (int mod_i = 0 ; mod_i < bilateralXYZStds.size() ; mod_i++) {
			std::cout << "bilateralXYZStds[" << mod_i << "]: " << bilateralXYZStds[mod_i] << std::endl;
		}


	}
};