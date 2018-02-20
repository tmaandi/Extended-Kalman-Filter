#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using namespace std;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /*
   * Calculating the RMSE here.
   */
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  bool validData = false;

  if ((estimations.size()!=0) && (estimations.size() == ground_truth.size()))
  {
      validData = true;
  }
  else
  {
      cout<<"CalculateRMSE(): Error - Incorrect vector length(s) \n";
  }

  if (validData == true)
  {
      VectorXd residual(4);
      VectorXd residualSum(4);
      residualSum << 0,0,0,0;
      /*
       * accumulate squared residuals
       */
      for(int i=0; i < estimations.size(); ++i){
            residual = (estimations[i] - ground_truth[i]);
            residual = residual.array().square();
            residualSum += residual;
      }
      /*
       * calculating mean of the residuals
       */
      residualSum = residualSum/estimations.size();

      /*
       * calculate the squared root
       */
      rmse = residualSum.array().sqrt();
  }

  return rmse;

}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /*
   * Calculate a Jacobian here.
   */
  MatrixXd Hj(3,4);

  /*
   * retrieve the state variables
   */
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  float t1 = px*px + py*py;
  float t2 = sqrt(t1);
  float t3 = t1*t2;

  /*
   * check division by zero
   */
  if (fabs(t1) < 0.0001)
  {
    cout<<"CalculateJacobian(): Error - Division by 0 \n";
  }
  else
  {
    Hj << px/t2, py/t2, 0, 0,
         -py/t1, px/t1, 0, 0,
          py*(vx*py - vy*px)/t3, px*(vy*px - vx*py)/t3, px/t2, py/t2;

  }

  return Hj;

}
