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
	lambda_ = 3 - n_x_;

	// initialize weights
	weights_ = VectorXd(2 * n_aug_ + 1);
	weights_.fill(0.0);

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
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.

	first measurement
	- initialize x_
	- initialize P_
	- store time

	- calculate dt
	- send measurement to correct update process

  */

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
		// Process Radar measurement
		if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
		}
		// Proces Lidar measurement
		else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
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
}
