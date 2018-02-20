#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  /*
   * initializing matrices
   */
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  /*
   * measurement covariance matrix - laser
   */
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  /*
   * measurement covariance matrix - radar
   */
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /*
   * initializing measurement noise matrices
   */
  H_laser_ << 1,0,0,0,
      0,1,0,0;

  Hj_ << 0,0,0,0,
      0,0,0,0,
      0,0,0,0;

}

/*
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {

    /*
     * Kalman Parameters Initialization
     */
    MatrixXd F, P, Q, R, H;
    VectorXd x(4);

    F = MatrixXd(4,4);
    P = MatrixXd(4,4);
    Q = MatrixXd(4,4);

    ekf_.F_ = MatrixXd(4,4);
    ekf_.P_ = MatrixXd(4,4);
    ekf_.Q_ = MatrixXd(4,4);

    ekf_.x_ = VectorXd(4);

    F << 1,0,0,0,
        0,1,0,0,
        0,0,1,0,
        0,0,0,1;

    P << 1,0,0,0,
        0,1,0,0,
        0,0,1000,0,
        0,0,0,1000;

    Q << 0,0,0,0,
        0,0,0,0,
        0,0,0,0,
        0,0,0,0;

    cout << "EKF: " << endl;

    previous_timestamp_= measurement_pack.timestamp_;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

      float rho = measurement_pack.raw_measurements_[0];
      float theta = measurement_pack.raw_measurements_[1];

      /*
       * Convert radar from polar to cartesian coordinates and initialize state.
       */
      float px = rho*cos(theta);
      float py = rho*sin(theta);

      x << px, py, 0.0, 0.0;

      R = MatrixXd(3,3);
      H = MatrixXd(3,4);

      ekf_.R_ = MatrixXd(3,3);
      ekf_.H_ = MatrixXd(3,4);

      R = R_radar_;
      H = Hj_;

      /*
       * Initializing Kalman Filter
       */
      ekf_.Init(x,P,F,H,R,Q);

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /*
       * Initialize state.
       */

      float px = measurement_pack.raw_measurements_[0];
      float py = measurement_pack.raw_measurements_[1];

      x << px, py, 0.0, 0.0;

      R = MatrixXd(2,2);
      H = MatrixXd(2,4);

      ekf_.R_ = MatrixXd(2,2);
      ekf_.H_ = MatrixXd(2,4);

      R = R_laser_;
      H = H_laser_;

      /*
       * Initializing Kalman Filter
       */
      ekf_.Init(x,P,F,H,R,Q);

    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }
  /*
   * Calculating delta time between two measurements
   */
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  float dt_2 = dt*dt;
  float dt_3 = dt_2*dt;
  float dt_4 = dt_3*dt;

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /*
   * Update F matrix to include delta time
   */
  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;

  /*
   * Update the process noise covariance matrix.
   */

  float noise_ax = 9.0;
  float noise_ay = 9.0;

  ekf_.Q_ << dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
             0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
             dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
             0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /*
   * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates

    ekf_.R_ = R_radar_;

    Hj_ = tools.CalculateJacobian(ekf_.x_);

    ekf_.H_ = Hj_;

    VectorXd z(3);
    z = measurement_pack.raw_measurements_;
    ekf_.UpdateEKF(z);

  } else {
    // Laser updates
    ekf_.R_ = R_laser_;

    ekf_.H_ = H_laser_;

    VectorXd z(2);
    z = measurement_pack.raw_measurements_;
    ekf_.Update(z);

  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
