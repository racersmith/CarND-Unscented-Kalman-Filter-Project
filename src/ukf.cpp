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
  std_a_ = 0.16;  // Bike acceleration.  0-6 m/s in 3 seconds -> 2 m/s^2 or 0.2g

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.29; // Bike turn rate.  Complete circle in 5 seconds.  Opposite circle in 1 second -> 2.5 rad/s^2

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
	weights_.fill(0.5/(lambda_ + n_aug_));
	weights_(0) = lambda_/(lambda_ + n_aug_);

	// Initialize sigma points
	Xsig_pred_ = MatrixXd(n_x_, n_sig_);
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
				0.2, 0, 0, 0, 0,
				0, 0.2, 0, 0, 0,
				0, 0, 2, 0, 0,
				0, 0, 0, 2, 0,
				0, 0, 0, 0, 2;

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
				0.15, 0, 0, 0, 0,
				0, 0.15, 0, 0, 0,
				0, 0, 2, 0, 0,
				0, 0, 0, 2, 0,
				0, 0, 0, 0, 2;

			// Extract laser data for readability
			double x = meas_package.raw_measurements_[0];
			double y = meas_package.raw_measurements_[1];
			x_ << x, 
				    y, 
						0.0, 
						0.0, 
						0.0;
		}

		// done initializing, no need to predict or update
		is_initialized_ = true;
	}

	// After initialization, standard UKF process loop
	else {
		// Calculate time step
		double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;     //dt - in seconds
		time_us_ = meas_package.timestamp_;
		//std::cout << "Process Measurment with dt =" << dt << std::endl << std::endl;

		// Predict
		Prediction(dt);

		// Update Radar
		if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
			//std::cout << "Radar" << std::endl;
			UpdateRadar(meas_package);
		}
		// Update Lidar
		else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
			//std::cout << "Lidar" << std::endl;
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


	// 1. Generate Augmented Sigma Points
	// ----------------------------------

	// Create augmented mean state
	VectorXd x_aug = VectorXd(n_aug_);
	x_aug.fill(0.0);
	x_aug.head(n_x_) = x_;

	// Create augmented covariance matrix
	MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
	P_aug.fill(0.0);
	P_aug.topLeftCorner(n_x_, n_x_) = P_;		// original covariance in upper left
	P_aug(5, 5) = std_a_*std_a_;						// variance of acceleration
	P_aug(6, 6) = std_yawdd_*std_yawdd_;		// variance of yaw acceleration

	// Calculate common values
	MatrixXd L = P_aug.llt().matrixL();		// sqrt of P_aug
	L *= sqrt(lambda_ + n_aug_);					// sigma point calculation sqrt term

	// Calculate sigma points
	MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig_);
	Xsig_aug.fill(0.0);
	Xsig_aug.col(0) = x_aug;
	for (int i = 0; i < n_aug_; i++) {
		Xsig_aug.col(i + 1)          = x_aug + L.col(i);
		Xsig_aug.col(i + 1 + n_aug_) = x_aug - L.col(i);
	}


	// 2. Predict Sigma Points
	// -----------------------

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

		// Predicted state values
		double px_p = px;		// x position
		double py_p = py;		// y position
		double v_p = v;			// heading velocity
		double psi_p = psi;		// yaw
		double psi_d_p = psi_d;	// yaw rate

		// Calculate common terms
		double cos_psi = cos(psi);
		double sin_psi = sin(psi);

		// Add process to predicted sigma point
		if (fabs(psi_d) > 0.0001) {
			// Calculate common terms
			double angle = psi + psi_d*delta_t;
			double v_psi_d = v / psi_d;

			// Add process to predicted sigma point
			px_p += v_psi_d * ( sin(angle) - sin_psi);
			py_p += v_psi_d * (-cos(angle) + cos_psi);
			psi_p += psi_d * delta_t;
		}
		// Avoid division by zero
		else {
			double v_dt = v * delta_t;
			px_p += v_dt * cos_psi;
			py_p += v_dt * sin_psi;
		}

		// Add noise to predicted sigma point
		double half_dt2 = 0.5 * delta_t * delta_t;
		px_p += half_dt2 * nu_a * cos_psi;
		py_p += half_dt2 * nu_a * sin_psi;
		v_p += delta_t * nu_a;
		psi_p += half_dt2 * nu_psi_dd;
		psi_d_p += delta_t * nu_psi_dd;

		// Write predicted sigma points into right column
		Xsig_pred_(0, i) = px_p;
		Xsig_pred_(1, i) = py_p;
		Xsig_pred_(2, i) = v_p;
		Xsig_pred_(3, i) = psi_p;
		Xsig_pred_(4, i) = psi_d_p;
	}
	//std::cout << "Prediction Xsig_pred_ (psi)" << std::endl << Xsig_pred_.row(3) << std::endl << std::endl;

	
	// 3. Predicted state mean and covariance
	// --------------------------------------

	// Predicted state mean
	x_.fill(0.0);
	for (int i = 0; i<n_sig_; i++) {
		x_ += weights_(i) * Xsig_pred_.col(i);
	}
	//std::cout << "Prediction x_" << std::endl << x_ << std::endl << std::endl;

	// Predicted state covariance matrix
	P_.fill(0.0);
	for (int i = 0; i<n_sig_; i++) {
		// Difference vector
		VectorXd x_diff = Xsig_pred_.col(i) - x_;

		//
		// TODO the heading angle is diverging.
		//
		x_diff(3) = NormalizeAngle(x_diff(3));

		// Calculate covariance
		P_ += weights_(i) * x_diff * x_diff.transpose();
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

	int n_z = 2;


	// 1. Transform Sigma Points to Measurement Space
	// ----------------------------------------------

	MatrixXd Z_sig = MatrixXd(n_z, n_sig_);
	Z_sig.row(0) = Xsig_pred_.row(0);
	Z_sig.row(1) = Xsig_pred_.row(1);
	//std::cout << "Z_sig" << std::endl << Z_sig << std::endl << std::endl;


	// 2. Calculate mean of predicted measurement
	// ------------------------------------------
	VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0.0);
	for (int i = 0; i < n_sig_; i++) {
		z_pred += weights_(i)*Z_sig.col(i);
	}
	
	//std::cout << "Lidar z_pred" << std::endl << z_pred << std::endl << std::endl;

	// Calculate covariance of predicted Lidar measurement
	MatrixXd S = MatrixXd(n_z, n_z);
	S.fill(0.0);
	for (int i = 0; i < n_sig_; i++) {
		VectorXd z_diff = Z_sig.col(i) - z_pred;
		S += weights_(i) * z_diff * z_diff.transpose();
	}

	// Add sensor noise
	S(0, 0) += std_laspx_*std_laspx_;
	S(1, 1) += std_laspy_*std_laspy_;

	//std::cout << "Lidar S" << std::endl << S << std::endl << std::endl;


	// 3. Update State
	// ---------------

	// Calculate cross correlation matrix, Tc
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

	// Calculate common term
	MatrixXd Si = S.inverse();

	// Kalman gain K
	MatrixXd K = Tc * Si;

	// Residual of measurement to prediction
	VectorXd z_diff = meas_package.raw_measurements_ - z_pred;
	
	// Update state and covariance
	x_ += K * z_diff;
	P_ -= K * S * K.transpose();


	// 4. Calculate Consistency NIS
	// ----------------------------
	NIS_laser_ = z_diff.transpose() * Si * z_diff;
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
	int n_z = 3;

	/**
	* Predict Radar Measurement
	*/

	// 1. Transform Sigma Points to Measurement Space
	// ----------------------------------------------

	MatrixXd Z_sig = MatrixXd(n_z, n_sig_);
	Z_sig.fill(0.0);

	for (int i = 0; i < n_sig_; i++) {
		// extract values for readibility
		double px = Xsig_pred_(0, i);
		double py = Xsig_pred_(1, i);
		double v = Xsig_pred_(2, i);
		double psi = Xsig_pred_(3, i);

		// measurement model
		// r, Range
		double range = sqrt(px*px + py*py);
		Z_sig(0, i) = range;	
		
		// phi, Bearing
		// Avoid atan2 domain errors at +/- inf
		if (fabs(px) > 0.0001) {
			Z_sig(1, i) = atan2(py, px);	
		}
		else if(py > 0.0001) {
			Z_sig(1, i) = M_PI / 2.0;
		}
		else if(py < -0.0001) {
			Z_sig(1, i) = -M_PI / 2.0;
		}
		else {
			Z_sig(1, i) = 0.0;
		}

		// r_dot, Range Rate
		// Avoid Div0 when range is 0
		if (range > 0.0001) {
			Z_sig(2, i) = (px*cos(psi)*v + py*sin(psi)*v) / range;
		}
		else {
			Z_sig(2, i) = 0.0;
		}
	}

	/*std::cout << 
		"Radar Z_sig" << std::endl << 
		Z_sig.row(0) << std::endl << 
		Z_sig.row(1) << std::endl << 
		Z_sig.row(2) << std::endl << std::endl;*/


	// 2. Predict Measurement
	// ----------------------

	// Calculate mean of predicted measurement
	VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0.0);
	for (int i = 0; i < n_sig_; i++) {
		z_pred += weights_(i) * Z_sig.col(i);
	}

	//std::cout << "Radar z_pred" << std::endl << z_pred << std::endl << std::endl;

	// Calculate covariance of predicted radar measurement
	MatrixXd S = MatrixXd(n_z, n_z);
	S.fill(0.0);
	for (int i = 0; i < n_sig_; i++) {
		VectorXd z_diff = Z_sig.col(i) - z_pred;
		z_diff(1) = NormalizeAngle(z_diff(1));    // normalize radar bearing, phi
		S += weights_(i) * z_diff * z_diff.transpose();
	}

	//std::cout << "Radar S" << std::endl << S << std::endl << std::endl;

	// Add sensor noise
	S(0, 0) += std_radr_*std_radr_;
	S(1, 1) += std_radphi_*std_radphi_;
	S(2, 2) += std_radrd_*std_radrd_;


	// 3. Update State
	// ---------------

	// Calculate cross correlation matrix, Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z);
	Tc.fill(0.0);
	for (int i = 0; i < n_sig_; i++) {
		// Residual of state
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		x_diff(3) = NormalizeAngle(x_diff(3));  // state heading

		// Residual of predicted measurment and sigma point
		VectorXd z_diff = Z_sig.col(i) - z_pred;
		z_diff(1) = NormalizeAngle(z_diff(1));  // radar bearing

		Tc += weights_(i) * x_diff * z_diff.transpose();
	}

	// Calculate Common Term
	MatrixXd Si = S.inverse();

	// Kalman gain K
	MatrixXd K = Tc * Si;

	// Residual of measurement to prediction
	VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

	//std::cout << "Radar z" << std::endl << meas_package.raw_measurements_ << std::endl << std::endl;

	//std::cout << z_diff(1) << " -> ";
	z_diff(1) = NormalizeAngle(z_diff(1));  // bearing
	//std::cout << z_diff(1) << std::endl;

	// Update state and covariance
	x_ += K * z_diff;
	P_ -= K * S * K.transpose();


	// 4. Calculate Consistency NIS
	// ----------------------------
	NIS_radar_ = z_diff.transpose() * Si * z_diff;
}

double UKF::NormalizeAngle(double angle) {
	//std::cout << angle << " -> ";
	// Constrain to less than pi
	while (angle >  M_PI) angle -= 2.0*M_PI;

	// Constrain to greater than -pi
	while (angle < -M_PI) angle += 2.0*M_PI;
	//std::cout << angle << std::endl;
	return angle;
}