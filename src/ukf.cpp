#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;  // Bike acceleration.  0-6 m/s in 3 seconds -> 2 m/s^2 or 0.2g

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 2.5; // Bike turn rate.  Complete circle in 5 seconds.  Opposite circle in 1 second -> 2.5 rad/s^2

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
	// first measurement flag
	is_initialized_ = false;

	// size of state vector
	n_x_ = x_.size();

	// size of augmented state vector
	n_aug_ = n_x_ + 2;  // +1 for std_a_ and +1 for std_yawdd_

	// define spreading parameters
	lambda_ = 3 - n_aug_;

	// initialize weights
	weights_ = VectorXd(2 * n_aug_ + 1);
	weights_.fill(0.0);

	// Initialize sigma points
	Xsig_pred_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
	Xsig_pred_.fill(0.0);

	// Initialize NIS for radar and lidar
	NIS_radar_ = 0.0;
	NIS_laser_ = 0.0;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
	// Initialize with first measurement
	if (!is_initialized_) {
		time_us_ = meas_package.timestamp_;

		// Initialize with Radar measurement
		if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
			// Set initial covariance
			P_ << 
				0.6, 0, 0, 0, 0,
				0, 0.6, 0, 0, 0,
				0, 0, 100, 0, 0,
				0, 0, 0, 100, 0,
				0, 0, 0, 0, 100;

			// Extract radar values for readability
			double range = meas_package.raw_measurements_[0];
			double bearing = meas_package.raw_measurements_[1];

			// convert to cartesion coordinates
			x_ <<
				range * cos(bearing),
				range * sin(bearing),
				0,
				0,
				0;
		}
		// Initialize with Lidar measurement
		else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
			// Set initial covariance
			P_ <<
				0.3, 0, 0, 0, 0,
				0, 0.3, 0, 0, 0,
				0, 0, 100, 0, 0,
				0, 0, 0, 100, 0,
				0, 0, 0, 0, 100;

			// Extract laser data for readability
			double x = meas_package.raw_measurements_[0];
			double y = meas_package.raw_measurements_[1];
			x_ << x, y, 0.0, 0.0, 0.0;
		}

		// done initializing, no need to predict or update
		is_initialized_ = true;
	}

	// After initialization, standard UKF process loop
	else {
		// Calculate time step
		double dt = meas_package.timestamp_ - time_us_;
		time_us_ = meas_package.timestamp_;

		// Predict
		Prediction(dt);

		// Update Radar
		if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
			UpdateRadar(meas_package);
		}
		// Update Lidar
		else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
			UpdateLidar(meas_package);
		}
	}

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */


	/**
	Generate augmented sigma points
	*/

	// Create augmented mean state
	VectorXd x_aug = VectorXd(n_aug_);
	x_aug.fill(0.0);
	x_aug.head(n_x_) = x_;

	// Create augmented covariance matrix
	MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
	P_aug.fill(0.0);
	P_aug.topLeftCorner(n_x_, n_x_) = P_;    // original covariance in upper left
	P_aug(n_x_, n_x_) = std_a_*std_a_;    // variance of acceleration
	P_aug(n_x_+1, n_x_+1) = std_yawdd_*std_yawdd_;  // variance of yaw acceleration

	// Calculate common values
	MatrixXd L = P_aug.llt().matrixL();  // sqrt of P_aug
	L *= sqrt(lambda_ + n_aug_);  // sigma point calculation sqrt term

	// Calculate sigma points
	MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
	Xsig_aug.fill(0.0);
	Xsig_aug.col(0) = x_aug;
	for (int i = 0; i < n_aug_; i++) {
		Xsig_aug.col(i + 1)          = x_aug + L.col(i);
		Xsig_aug.col(i + 1 + n_aug_) = x_aug - L.col(i);
	}


	/**
	Predict sigma points
	*/

	// Iterate through each sigma point and predict
	for (int i = 0; i < Xsig_aug.cols(); i++) {
		// Extract sigma point values for readability
		double px = Xsig_aug(0, i);					// x position
		double py = Xsig_aug(1, i);					// y position
		double v = Xsig_aug(2, i);					// heading velocity
		double psi = Xsig_aug(3, i);				// yaw
		double psi_d = Xsig_aug(4, i);			// yaw rate
		double nu_a = Xsig_aug(5, i);				// aceleration variance
		double nu_psi_dd = Xsig_aug(6, i);	// Yaw acceleration variance

		// Extract state values for readability
		double px_p = x_(0);		// x position
		double py_p = x_(1);		// y position
		double v_p = x_(2);			// heading velocity
		double psi_p = x_(3);		// yaw
		double psi_d_p = x_(4);	// yaw rate

		// Add process
		if (fabs(psi_d) > 0.0001) {
			double angle = psi + psi_d*delta_t;
			px_p += v / psi_d*(sin(angle) - sin(psi));
			py_p += v / psi_d*(-cos(angle) + cos(psi));
			psi_p += psi_d * delta_t;
		}
		// Avoid division by zero
		else {
			px_p += v * cos(psi) * delta_t;
			py_p += v * sin(psi) * delta_t;
		}

		// Add noise
		double dt2 = delta_t * delta_t;
		px_p += 0.5 * dt2 * cos(psi) * nu_a;
		py_p += 0.5 * dt2 * sin(psi) * nu_a;
		v_p += delta_t * nu_a;
		psi_p += 0.5 * dt2 * nu_psi_dd;
		psi_d_p += delta_t * nu_psi_dd;

		// write predicted sigma points into right column
		Xsig_pred_(0, i) = px_p;
		Xsig_pred_(1, i) = py_p;
		Xsig_pred_(2, i) = v_p;
		Xsig_pred_(3, i) = psi_p;
		Xsig_pred_(4, i) = psi_d_p;
	}


	/**
	Predicted state mean and covariance
	*/
	// Set weights
	weights_(0) = lambda_ / (lambda_ + n_aug_);
	for (int i = 1; i<Xsig_pred_.cols(); i++) {
		weights_(i) = 1 / (2 * (lambda_ + n_aug_));
	}

	// Predict state mean   
	for (int i = 0; i<Xsig_pred_.cols(); i++) {
		x_ += weights_(i)*Xsig_pred_.col(i);
	}

	// Predict state covariance matrix
	for (int i = 0; i<Xsig_pred_.cols(); i++) {
		// Difference vector
		VectorXd A = Xsig_pred_.col(i) - x_;

		// Normalize angle
		while (A(3) > M_PI) A(3) -= 2.0*M_PI;
		while (A(3) < -M_PI) A(3) += 2.0*M_PI;

		// Calculate covariance
		P_ += weights_(i)*A*A.transpose();
	}
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

	// Predict Lidar measurement
	// Update State

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

	// Predict Radar measurement
	// Update state
}
