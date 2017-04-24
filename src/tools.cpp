#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
	VectorXd rmse(4);
	rmse << 0, 0, 0, 0;

	int est_size = estimations.size();

	// Check for valid vectors
	if (est_size != ground_truth.size() || est_size == 0) {
		std::cout << "Invalid vectors passed to CalculateRMSE." << std::endl;
		std::cout << "  est size: " << est_size << std::endl;
		std::cout << "  truth size: " << ground_truth.size() << std::endl;
		return rmse;
	}
	else {
		// Calculate RMSE of valid vectors
		for (int i = 0; i < est_size; i++) {
			VectorXd residual = estimations[i] - ground_truth[i];
			residual = residual.array() * residual.array();
			rmse += residual;
		}
		rmse = rmse / est_size;
		rmse = rmse.array().sqrt();
		return rmse;
	}
}
