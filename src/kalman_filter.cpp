#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  DONE:
    * predict the state
  */
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  DONE:
    * update the state by using Kalman Filter equations
  */
  MatrixXd Ht = H_.transpose();

  VectorXd y = z - (H_ * x_);
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd K = P_ * Ht * S.inverse();

  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  DONE:
    * update the state by using Extended Kalman Filter equations
  */

  // Define local variables to work on
  float px = x_[0];
  float py = x_[1];
  float vx = x_[2];
  float vy = x_[3];
  
  float rho = sqrt(pow(px,2)+pow(py,2));
  float phi = atan2(py,px); 
  float rho_dot;

  //Checking division by zero.
  if (fabs(rho) < 0.0001) {
    rho_dot = 0;
  } else {
    rho_dot = (vx*px + vy*py)/rho;
  }

  // Vector containing predicted values.
  VectorXd z_pred (3);
  z_pred << rho, phi, rho_dot;

  // Difference between actual radar values and prediction.
  VectorXd y = z - z_pred;

  //Ensuring that the phi value always stays between -Pi and +Pi.
  while (y[1] > M_PI){
     y[1] -= (2 * M_PI);
  }
  while (y[1] < - M_PI){
     y[1] += (2 * M_PI);
  }

  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}