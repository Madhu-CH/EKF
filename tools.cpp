#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
  VectorXd rmse = VectorXd(4);
  rmse << 0,0,0,0;

  unsigned int t = estimations.size();
  
  if(estimations.size() == 0 || ground_truth.size() == 0){
     //cout << "Estimations or ground truth vector had size 0" << endl;
     return rmse;
  }

  if (t!= ground_truth.size()){

    return rmse;
  }
  
  for(unsigned int n =0; n < estimations.size(); n++){
      VectorXd residuals = estimations[n] - ground_truth[n];
      residuals = residuals.array() * residuals.array();
      rmse += residuals;
   }

  
  
   rmse =  rmse / t;
   rmse = rmse.array().sqrt();

   return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */ 
    MatrixXd Hj = MatrixXd::Zero(3,4);
  
   float px = x_state(0); // x pos
   float py = x_state(1); // y pos
   float vx = x_state(2); // x velocity
   float vy = x_state(3); // y velocity

   // calculate repetitive operations
  float c1 = px*px+py*py;
  float c2 = sqrt(c1);
  float c3 = (c1*c2);
  
   if((fabs(c1)<  0.00001) || (px == 0 || py == 0)){
      //cout << "Can't divide by 0" << endl;
      return Hj;
   }
	Hj << (px/c2), (py/c2), 0, 0,
       -(py/c1), (px/c1), 0, 0,
        py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;
  
   return Hj;
}


MatrixXd Tools::CalculateCovariant(const double dt, const double noise_ax, const double noise_ay){

   MatrixXd Q_(4,4);

   float dt_2 = dt * dt;
   float dt_3 = dt_2 * dt;
   float dt_4 = dt_3 * dt;

   Q_ = MatrixXd(4, 4);
   Q_ << (dt_4/4)*noise_ax, 0, ((dt_3)/2)*noise_ax, 0,
      0, ((dt_4)/4)*noise_ay, 0, ((dt_3)/2)*noise_ay,
      ((dt_3)/2)*noise_ax, 0, (dt_2)*noise_ax, 0,
      0, ((dt_3)/2)*noise_ay, 0, (dt_2)*noise_ay;

   return Q_;

}
