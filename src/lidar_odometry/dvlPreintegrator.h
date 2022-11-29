#pragma once

#include <gtsam/base/Vector.h>
#include <gtsam/geometry/Pose3.h>
#include <gtsam/geometry/Rot3.h>
#include <gtsam/geometry/Point3.h>

#include <array>
#include <bitset>
#include <cmath>
#include <numeric>
#include <utility>
#include <vector>
#include <Eigen/Dense>

inline Eigen::Matrix3d hat(const Eigen::Vector3d &w) {
    return (Eigen::Matrix3d() << 0, -w.z(), w.y(),
            w.z(), 0, -w.x(),
            -w.y(), w.x(), 0)
        .finished();
}

inline Eigen::Quaterniond expmap(const Eigen::Vector3d &w) {
    Eigen::AngleAxisd aa(w.norm(), w.stableNormalized());
    Eigen::Quaterniond q;
    q = aa;
    return q;
}

Eigen::Matrix3d right_jacobian(const Eigen::Vector3d &w) {
    static const double root2_eps = sqrt(std::numeric_limits<double>::epsilon());
    static const double root4_eps = sqrt(root2_eps);
    static const double qdrt720 = sqrt(sqrt(720.0));
    static const double qdrt5040 = sqrt(sqrt(5040.0));
    static const double sqrt24 = sqrt(24.0);
    static const double sqrt120 = sqrt(120.0);

    double angle = w.norm();
    double cangle = cos(angle);
    double sangle = sin(angle);
    double angle2 = angle * angle;

    double cos_term;
    // compute (1-cos(x))/x^2, its taylor expansion around 0 is 1/2-x^2/24+x^4/720+o(x^6)
    if (angle > root4_eps * qdrt720) {
        cos_term = (1 - cangle) / angle2;
    } else { // use taylor expansion to avoid singularity
        cos_term = 0.5;
        if (angle > root2_eps * sqrt24) { // we have to include x^2 term
            cos_term -= angle2 / 24.0;
        }
    }

    double sin_term;
    // compute (x-sin(x))/x^3, its taylor expansion around 0 is 1/6-x^2/120+x^4/5040+o(x^6)
    if (angle > root4_eps * qdrt5040) {
        sin_term = (angle - sangle) / (angle * angle2);
    } else {
        sin_term = 1.0 / 6.0;
        if (angle > root2_eps * sqrt120) { // we have to include x^2 term
            sin_term -= angle2 / 120.0;
        }
    }

    Eigen::Matrix3d hat_w = hat(w);
    return Eigen::Matrix3d::Identity() - cos_term * hat_w + sin_term * hat_w * hat_w;
}

struct DVLData {
  DVLData() = default;

  DVLData(double t, Eigen::Vector3d w, Eigen::Vector3d v)
    : t(t), w(w), v(v) {}

  double t;          // timestamp
  Eigen::Vector3d w; // gyro measurement
  Eigen::Vector3d v; // velocity measurement
};

class DVLPreintegrator
{
public:
    DVLPreintegrator(double noise_w, double noise_v, double noise_bg, double noise_bv) :
      cov_w(Eigen::Matrix3d::Identity()*noise_w*noise_w), cov_v(Eigen::Matrix3d::Identity()*noise_v*noise_v), 
      cov_bg(Eigen::Matrix3d::Identity()*noise_bg*noise_bg), cov_bv(Eigen::Matrix3d::Identity()*noise_bv*noise_bv)
    {
      reset();
    }

    struct Delta {
        double t;
        Eigen::Quaterniond q;
        Eigen::Vector3d p;
        Eigen::Matrix<double, 6, 6> cov; // only on q and p
    };

    struct Jacobian {
        Eigen::Matrix3d dq_dbg;
        Eigen::Matrix3d dp_dbg;
        Eigen::Matrix3d dp_dbv;
    };

    Eigen::Matrix3d cov_w; // continuous noise covariance
    Eigen::Matrix3d cov_v;
    Eigen::Matrix3d cov_bg; // continuous random walk noise covariance
    Eigen::Matrix3d cov_bv;

    Eigen::Vector3d bg{Eigen::Vector3d::Zero()}; // assume zero initial biases
    Eigen::Vector3d bv{Eigen::Vector3d::Zero()};

    Delta delta;
    Jacobian jacobian;

    void print(std::ostream &stream = std::cout) const
    {
        stream << "Preintegrated DVL Measurements: " << "\n";
        stream << "  deltaTij [ " << delta.t << " ]" << "\n";
        stream << "  deltaRij [ " << delta.q.vec().transpose() << " " << delta.q.w() << " ]'" << "\n";
        stream << "  deltaPij [ " << delta.p.transpose() << " ]'" << "\n";
        stream << "  gyrobias [ " << bg.transpose() << " ]'" << "\n";
        stream << "  velobias [ " << bv.transpose() << " ]'" << "\n";
        stream << "  dvlPreintMeasCov " << "\n";
        stream << " [ " << delta.cov << " ]" << "\n";
    }

    void reset(Eigen::Vector3d bg_new = Eigen::Vector3d::Zero(), Eigen::Vector3d bv_new = Eigen::Vector3d::Zero())
    {
      bg = bg_new;
      bv = bv_new;
      delta.t = 0;
      delta.q.setIdentity();
      delta.p.setZero();
      delta.cov.setZero();

      jacobian.dp_dbg.setZero();
      jacobian.dp_dbg.setZero();
      jacobian.dp_dbv.setZero();
    }

    void integrateMeasurement(double dt, const DVLData &data)
    {
        assert(("dt must > 0") && (dt >= 0));

        // correct measurements with biases
        Eigen::Vector3d w = data.w - bg;
        Eigen::Vector3d v = data.v - bv;

        Eigen::Quaterniond q_full(expmap(w * dt));
        Eigen::Quaterniond q_half(expmap(0.5 * w * dt));

        // 1. Calculate measurement covariance matrix iteratively
        Eigen::Matrix<double, 6, 6> A;
        A.setIdentity();
        A.block<3, 3>(0, 0) = q_full.conjugate().matrix();
        A.block<3, 3>(3, 0) = -dt * delta.q.matrix() * hat(v);

        Eigen::Matrix<double, 6, 6> B;
        B.setZero();
        B.block<3, 3>(0, 0) = dt * right_jacobian(w * dt);
        B.block<3, 3>(3, 3) = dt * delta.q.matrix();

        Eigen::Matrix<double, 6, 6> white_noise_cov;
        double inv_dt = 1.0 / std::max(dt, 1.0e-7);
        white_noise_cov.setZero();
        white_noise_cov.block<3, 3>(0, 0) = cov_w * inv_dt;
        white_noise_cov.block<3, 3>(3, 3) = cov_v * inv_dt;

        delta.cov = A * delta.cov * A.transpose() + B * white_noise_cov * B.transpose();

        // 2. Calculate jacobians
        jacobian.dp_dbg -= dt * delta.q.matrix() * hat(v) * jacobian.dq_dbg;
        jacobian.dp_dbv -= dt * delta.q.matrix();
        jacobian.dq_dbg = expmap(w * dt).conjugate().matrix() * jacobian.dq_dbg - dt * right_jacobian(w * dt);

        // 3. Calculate state deltas
        Eigen::Quaterniond q_mid = delta.q * q_half;
        Eigen::Vector3d v_mid = q_mid * v;

        delta.t = delta.t + dt;
        delta.p = delta.p + dt * v_mid;
        delta.q = (delta.q * q_full).normalized();
    }

    gtsam::Pose3 deltaPoseij()
    {
        // rotate deltaPij into world frame
        gtsam::Pose3 deltaPoseij = gtsam::Pose3(gtsam::Rot3(delta.q.w(), delta.q.x(), delta.q.y(), delta.q.z()), gtsam::Point3(delta.p.x(), delta.p.y(), delta.p.z()));
        return deltaPoseij;
    }

    gtsam::Vector6 Sigmasij()
    {
        gtsam::Vector6 Sigmasij = (gtsam::Vector(6) << 1e+18, 1e+18, 1e+18, sqrt(delta.cov.diagonal()[3]), sqrt(delta.cov.diagonal()[4]), sqrt(delta.cov.diagonal()[5])).finished();
        // gtsam::Vector6 Sigmasij = (gtsam::Vector(6) << 1e+18, 1e+18, 1e+18, 1e-10, 1e-10, 1e-10).finished();
        return Sigmasij;
    }
};
