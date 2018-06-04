#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
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
  std_a_ = 1.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = M_PI/8;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
    is_initialized_=false;

    time_us_ = 0.0;
    n_x_ = 5;
    n_aug_ = 7;
    lambda_ = 3-n_aug_;

    weights_ =  VectorXd(2*n_aug_+1);

   previous_timestamp_ = 0.0;
//    double previous_timestamp_ = 0;
    NIS_radar = 0;
    NIS_laser = 0;
    //MatrixXd Xsig_pred = MatrixXd (n_x_, 2 * n_aug_+1);

    Xsig_pred = MatrixXd(n_x_, 2*n_aug_+1);

    P_ <<    1, 0, 0, 0, 0,
            0, 1, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;

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
  */
    if (!is_initialized_)
    {
        cout << "UKF: " <<endl;
        x_ <<1,1,1,1,1;
        previous_timestamp_ = meas_package.timestamp_;

        if (meas_package.sensor_type_== MeasurementPackage::RADAR)
        {
            float rho = meas_package.raw_measurements_[0];
            float phi = meas_package.raw_measurements_[1];
            float rho_dot = meas_package.raw_measurements_[2];
            float vx = rho_dot * cos(phi);
            float vy = rho_dot * sin(phi);
            x_[0] = rho * cos(phi);
            x_[1] = rho * sin(phi);
            x_[2] = sqrt(vx*vx+vy*vy);
            x_[3] = 0;
            x_[4] = 0;


        }
        else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
        {
            double px = meas_package.raw_measurements_[0];
            double py = meas_package.raw_measurements_[1];
            x_[0] = px;
            x_[1] = py;
            x_[2] = 0;
            x_[3] = 0;
            x_[4] = 0;
        }
        is_initialized_ = true;
        return;
    }
    double delta_t = (meas_package.timestamp_ - previous_timestamp_)/1000000.0;
    previous_timestamp_ = meas_package.timestamp_;
//Call Update step after prediction
    if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
        // Call Prediction step
        Prediction(delta_t);
        cout << "P1" <<endl;
        UpdateLidar(meas_package);
        cout << "P2" <<endl;
    }
    if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
        // Call Prediction step
        cout << delta_t <<endl;
        Prediction(delta_t);
        cout << "P3" <<endl;

        UpdateRadar(meas_package);
        cout << "P4" <<endl;
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

    // Generate sigma points

    VectorXd x_aug = VectorXd(n_aug_);
    MatrixXd P_aug = MatrixXd(n_aug_,n_aug_);
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_+1);
    Xsig_aug.fill(0.0);
    x_aug.head(n_x_) = x_;
    x_aug(5) = 0;
    x_aug(6) = 0;
    P_aug.fill(0.0);

    P_aug.topLeftCorner(5,5) = P_;
    P_aug(5,5) = std_a_ *std_a_;
    P_aug(6,6) = std_yawdd_ * std_yawdd_;
    MatrixXd L = P_aug.llt().matrixL();
    Xsig_aug.col(0) = x_aug;
    for (int i = 0; i<n_aug_; i++)
    {
        Xsig_aug.col(i+1) = x_aug +sqrt(lambda_+ n_aug_)* L.col(i);
        Xsig_aug.col(i+1+n_x_) = x_aug - sqrt(lambda_ + n_aug_) *L.col(i);
    }


    // Predict sigma points

    for (int i = 0; i<2*n_aug_+1; i++) {
        double p_x = Xsig_aug(0, i);
        double p_y = Xsig_aug(1, i);
        double v = Xsig_aug(2, i);
        double yaw = Xsig_aug(3, i);
        double yawd = Xsig_aug(4, i);
        double nu_a = Xsig_aug(5, i);
        double nu_yawdd = Xsig_aug(6, i);

        double px_p, py_p;
        if (fabs(yawd) > 0.001) {
            px_p = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
            py_p = p_y + v / yawd * (cos(yaw + yawd * delta_t));
        } else {
            px_p = p_x + v * delta_t * cos(yaw);
            py_p = p_y + v * delta_t * sin(yaw);
        }
        double v_p = v;
        double yaw_p = yaw + yawd * delta_t;
        double yawd_p = yawd;
        px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
        py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
        v_p = v_p + nu_a * delta_t;
        yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
        yawd_p = yawd_p + nu_yawdd * delta_t;

        Xsig_pred(0, i) = px_p;
        Xsig_pred(1, i) = py_p;
        Xsig_pred(2, i) = v_p;
        Xsig_pred(3, i) = yaw_p;
        Xsig_pred(4, i) = yawd_p;

    }


    //predict mean and covariance
     //set weights
         weights_(0) = lambda_/(lambda_+n_aug_);
        for (int i=1; i<2*n_aug_+1; i++)
        {  //2n+1 weights
            double a = 0.5/(n_aug_+lambda_);
            weights_(i) = a;
        }

        x_.fill(0.0);


        for (int i=0; i<2*n_aug_+1; i++)
        {
            x_ = x_+weights_(i)* Xsig_pred.col(i);
        }


        P_.fill(0.0);

    cout <<"P5"<<endl;


        for (int i = 0; i < 2 * n_aug_ + 1; i++)
        {  //iterate over sigma points

            // state difference
            cout << Xsig_pred <<endl;
            cout << endl;
            VectorXd x_diff = Xsig_pred.col(i) - x_;
            cout <<"P5.1" << endl;
            //angle normalization
            cout << x_diff(3);
            while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
            while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

            cout <<"P6"<<endl;
            P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
            cout <<"P7"<<endl;
            //cout << P_ <<endl;
            //cout <<"     "<<endl;
            //cout <<endl;
            //cout <<endl;
        }
    cout << "Prediction end " << endl;



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
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_+1);
    VectorXd z_pred = VectorXd(n_z);
    MatrixXd S = MatrixXd(n_z, n_z);
    MatrixXd R = MatrixXd(n_z, n_z);
    MatrixXd Tc = MatrixXd(n_x_, n_z);

    VectorXd z = meas_package.raw_measurements_;

    for (int i = 0; i<2*n_aug_+1; i++)
    {
        double p_x = Xsig_pred(0,i);
        double p_y = Xsig_pred(1,i);
        Zsig.col(i) << p_x, p_y;

    }


    z_pred.fill(0.0);
    for (int i=0; i<2*n_aug_+1; i++)
    {
        z_pred += weights_(i) * Zsig.col(i);
    }





    S.fill(0.0);
    for (int i=0; i<2*n_aug_+1; i++)
    {
        VectorXd z_diff = Zsig.col(i) - z_pred;
        S += weights_(i) * z_diff * z_diff.transpose();
    }



    R<<  std_laspx_*std_laspx_, 0,
         0, std_laspy_*std_laspy_;

    S += R;

    Tc.fill(0);
    for (int i=0; i<2 * n_aug_ + 1; i++)
    {
        VectorXd x_diff = Xsig_pred.col(i) - x_;

        VectorXd z_diff = Zsig.col(i) - z_pred;
        Tc += weights_(i) * x_diff * z_diff.transpose();
    }


    MatrixXd K = Tc * S.inverse();

    //update state mean and covariance matrix
    VectorXd z_diff_2 = z - z_pred;
    x_ += K*(z_diff_2);
    P_ -= K * S * K.transpose();

    NIS_laser = z_diff_2.transpose() * S.inverse() * z_diff_2;
    cout << "uUpdateLaser " << endl;

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
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
    VectorXd z_pred = VectorXd(n_z);
    MatrixXd S = MatrixXd(n_z,n_z);
    MatrixXd R = MatrixXd(n_z, n_z);
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    VectorXd z = meas_package.raw_measurements_;



    for (int i = 0; i<2*n_aug_+1; i++)
    {
        double p_x = Xsig_pred(0,i);
        double p_y = Xsig_pred(1,i);
        double v = Xsig_pred(2,i);
        double yaw = Xsig_pred(3,i);
        double rho = sqrt((p_x*p_x)+(p_y*p_y));
        double phi = 0.0;
        double rho_dot = 0.0;

        if ((fabs(p_x) < 0.001) && (fabs(p_y) < 0.001))
        {
            if (fabs(p_x) < 0.001)
            {p_x = 0.001;}
            else {p_y = 0.001;}
        }
        phi = atan2(p_y, p_x);

        if (fabs(rho) < 0.001)
        {
            rho = 0.001;
        }
        rho_dot = (p_x * cos(yaw) * v + p_y * sin(yaw) * v) / rho;

        Zsig(0,i)=rho;
        Zsig(1,i)= phi;
        Zsig(2,i)= rho_dot;

    }

    z_pred.fill(0.0);

    for (int i=0; i<2*n_aug_; i++)
    {
        z_pred = z_pred + weights_(i) * Zsig.col(i);

    }



    S.fill(0.0);
    for (int i=0; i<2*n_aug_+1; i++)
    {
        VectorXd z_diff = Zsig.col(i) - z_pred;

        //angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

        S += weights_(i) * z_diff * z_diff.transpose();
    }

    R.fill(0);
    R <<  std_radr_*std_radr_, 0, 0,
            0, std_radphi_*std_radphi_, 0,
            0, 0, std_radrd_*std_radrd_;

    //add measurement noise covariance matrix
    S += R;


    Tc.fill(0.0);
    for (int i= 0; i<2*n_aug_+1;i++)
    {
        VectorXd x_diff = Xsig_pred.col(i) - x_;
        while (x_diff(3)>M_PI)
        {
            x_diff(3) = x_diff(3)- 2.*M_PI;
        }
        while (x_diff(3)< -M_PI)
        {
            x_diff(3) = x_diff(3)+2.*M_PI;
        }

        VectorXd z_diff = Zsig.col(i) - z_pred;

        while (z_diff(1)>M_PI)
        {
            z_diff(1) = z_diff(1) - 2.*M_PI;
        }

        while (z_diff(1)< -M_PI)
        {
            z_diff(1) = z_diff(1) + 2.*M_PI;
        }

        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }

    MatrixXd K = Tc * S.inverse();
    VectorXd z_diff = z - z_pred;
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    x_ += K * z_diff;
    P_ -= K * S * K.transpose();

    NIS_radar = z_diff.transpose() * S.inverse() * z_diff;

    cout << "UpdateRadar " << endl;
}
