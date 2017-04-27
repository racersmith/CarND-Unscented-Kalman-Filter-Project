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
  std_yawdd_ = 0.2; // Bike turn rate.  Complete circle in 5 seconds.  Opposite circle in 1 second -> 2.5 rad/s^2

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

	// number of sigma points
	n_sig_ = 2 * n_aug_ + 1;

	// define spreading parameters
	lambda_ = 3 - n_aug_;

	// Calculate weights
	weights_ = VectorXd(n_sig_);
	weights_.fill(1 / (2 * (lambda_ + n_aug_)));
	weights_(0) = lambda_ / (lambda_ + n_aug_);

	// Initialize sigma points
	Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
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
				1, 0, 0, 0, 0,
				0, 1, 0, 0, 0,
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
				1, 0, 0, 0, 0,
				0, 1, 0, 0, 0,
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
		double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;     //dt - in seconds
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
	MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig_);
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
	for (int i = 0; i < n_sig_; i++) {
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

		// Write predicted sigma points into right column
		Xsig_pred_(0, i) = px_p;
		Xsig_pred_(1, i) = py_p;
		Xsig_pred_(2, i) = v_p;
		Xsig_pred_(3, i) = psi_p;
		Xsig_pred_(4, i) = psi_d_p;
	}


	/**
	Predicted state mean and covariance
	*/

	// Predicted state mean   
	for (int i = 0; i<n_sig_; i++) {
		x_ += weights_(i)*Xsig_pred_.col(i);
	}

	// Predicted state covariance matrix
	for (int i = 0; i<n_sig_; i++) {
		// Difference vector
		VectorXd A = Xsig_pred_.col(i) - x_;

		A(3) = NormalizeAngle(A(3));

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

	
	/**
	Predict Lidar measurement
	*/

	// Transform sigma points to measurement space
	int n_z = 2;    // Dimension of measurement
	MatrixXd Z_sig = Xsig_pred_.topRows(n_z);  // Transform is just the x, y terms

	// Calculate mean of predicted measurement
	VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0.0);
	for (int i = 0; i < n_sig_; i++) {
		z_pred += weights_(i) * Z_sig.col(i);
	}

	// Calculate covariance of predicted measurement
	MatrixXd S = MatrixXd(n_z, n_z);
	S.fill(0.0);
	for (int i = 0; i < n_sig_; i++) {
		// Residual
		VectorXd z_diff = Z_sig.col(i) - z_pred;

		// Calculate covariance
		S = S + weights_(i) * z_diff * z_diff.transpose();
	}


	/**
	Update State
	*/

	// Calculate cross correlation matrix
	MatrixXd Tc = MatrixXd(n_x_, n_z);
	Tc.fill(0.0);
	for (int i = 0; i < n_sig_; i++) {
		// Residual of state
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		x_diff(3) = NormalizeAngle(x_diff(3));  // heading

		// Residual of predicted measurment and sigma point
		VectorXd z_diff = Z_sig.col(i) - z_pred;

		Tc += weights_(i) * x_diff * z_diff.transpose();
	}
	// Kalman gain K
	MatrixXd K = Tc * S.inverse();

	// Residual of measurement to prediction
	VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

	// Update state and covariance
	x_ += K * z_diff;
	P_ -= K * S * K.transpose();
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

	// Transform sigma points to measurement space
	int n_z = 3;    // Dimension of measurement
	MatrixXd Z_sig = MatrixXd(n_z, n_sig_);
	for (int i = 0; i < n_sig_; i++) {
		
		// extract values for better readibility
		double p_x = Xsig_pred_(0, i);
		double p_y = Xsig_pred_(1, i);
		double v = Xsig_pred_(2, i);
		double yaw = Xsig_pred_(3, i);

		// calculate common values
		double v_cos_yaw = cos(yaw)*v;
		double v_sin_yaw = sin(yaw)*v;

		// measurement model
		Z_sig(0, i) = sqrt(p_x*p_x + p_y*p_y);												// r, Range

		// TODO domain error
		Z_sig(1, i) = atan2(p_y, p_x);																// phi, Bearing

		// TODO div0
		Z_sig(2, i) = (p_x*v_cos_yaw + p_y*v_sin_yaw) / Z_sig(0, i);	// r_dot, Range Rate
	}
	
	// Calculate mean of predicted measurement
	VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0.0);
	for (int i = 0; i < n_sig_; i++) {
		z_pred += weights_(i) * Z_sig.col(i); 
	}

	// Calculate covariance of predicted radar measurement
	MatrixXd S = MatrixXd(n_z, n_z);
	S.fill(0.0);

	for (int i = 0; i < n_sig_; i++) {
		// Residual
		VectorXd z_diff = Z_sig.col(i) - z_pred;

		// Angle normalization
		z_diff(1) = NormalizeAngle(z_diff(1));

		// Calculate covariance
		S += weights_(i) * z_diff * z_diff.transpose();
	}


	/**
	Update State
	*/

	// Calculate cross correlation matrix, Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z);
	Tc.fill(0.0);
	for (int i = 0; i < n_sig_; i++) {
		// Residual of state
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		x_diff(3) = NormalizeAngle(x_diff(3));  // heading

		// Residual of predicted measurment and sigma point
		VectorXd z_diff = Z_sig.col(i) - z_pred;
		z_diff(1) = NormalizeAngle(z_diff(1));  // bearing

		Tc += weights_(i) * x_diff * z_diff.transpose();
	}
	// Kalman gain K
	MatrixXd K = Tc * S.inverse();

	// Residual of measurement to prediction
	VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

	z_diff(1) = NormalizeAngle(z_diff(1));  // bearing

	// Update state and covariance
	x_ += K * z_diff;
	P_ -= K * S * K.transpose();
}

double UKF::NormalizeAngle(double angle) {
	std::cout << angle << " -> ";
	// Constrain to less than pi
	while (angle >  M_PI) angle -= 2.0*M_PI;

	// Constrain to greater than -pi
	while (angle < -M_PI) angle += 2.0*M_PI;
	std::cout << angle << std::endl;
	return angle;
}