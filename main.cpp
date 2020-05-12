#include "Eigen/Eigen"
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

class Solver
{
public:
    Solver(Eigen::MatrixXd _C, Eigen::MatrixXd _X, Eigen::MatrixXd _A, Eigen::MatrixXd _b) : C(_C), X(_X), A(_A), b(_b) {}
    bool run(const int _num_iter, const double t = 1.0)
    {
        Eigen::MatrixXd A_T = A.transpose();
        bool flag = true;
        for (int i = 0; i < _num_iter; ++i)
        {
            Eigen::VectorXd diag = A * X - b;
            for (int j = 0; j < diag.size(); ++j)
            {
                if (diag(j) > 0)
                {
                    std::cout << "Algorithm break due to the failure of the restriction check!" << std::endl;
                    flag = false;
                    best_point = Eigen::Infinity * Eigen::MatrixXd::Ones(X.rows(), X.cols());
                    break;
                }
            }
            if (!flag)
                break;
            // the first order deriv
            Eigen::VectorXd diag_upsidedown = diag;
            for (int j = 0; j < diag_upsidedown.size(); ++j)
                diag_upsidedown(j) = 1.0 / (diag_upsidedown(j) + 1e-6);
            Eigen::MatrixXd first_order_mat = t * C - A_T * diag_upsidedown;
            // the second order deriv and hessian mat
            Eigen::MatrixXd diagMat = Eigen::MatrixXd::Zero(A.rows(), A.rows());
            for (int j = 0; j < A.rows(); ++j)
                diagMat(j, j) = pow(diag_upsidedown(j), 2);
            Eigen::MatrixXd hessianMat = A_T * diagMat * A;
            Eigen::MatrixXd hessian_inv;
            bool inv_able;
            inv_able = invertible(hessianMat, hessian_inv);
            if (!inv_able)
            {
                std::cout << "Hessian Matrix can not be inversed!!!" << std::endl;
                break;
            }
            X = X - hessian_inv * first_order_mat;
            best_point = X;
        }
        return flag;
    }

    Eigen::MatrixXd get_last_iter_point()
    {
        return best_point;
    }

private:
    Eigen::MatrixXd C, X, A, b;
    Eigen::MatrixXd best_point;

    bool invertible(Eigen::MatrixXd mat, Eigen::MatrixXd &inv)
    {
        inv = mat.inverse();
        Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(mat.rows(), mat.cols());
        const double epsilon = 0.00001;
        return abs((identity - mat * inv).sum()) < epsilon;
    }
};

int main(int argc, char const *argv[])
{
    Eigen::MatrixXd C(2, 1), A(4, 2), b(4, 1), X(2, 1);
    C << 1, 1;
    A << 1, 2,
        2, 1,
        -1, 0,
        0, -1;
    b << 10, 10, 0, 0;
    double x0 = 3, y0 = 3, t = 0.1;
    bool flag = true, pre_reset = false;
    Eigen::MatrixXd pre_best(2, 1), valid_best_point(2, 1);
    pre_best << x0, y0;
    int iter_count = 0, max_loop = 8000;
    std::vector<std::pair<double,double>> points;
    points.push_back(std::make_pair(x0,y0));
    while (flag && max_loop-- > 0)
    {
        X << x0, y0;
        Solver solver(C, X, A, b);
        flag = solver.run(100, t = t);
        if (flag)
            valid_best_point = solver.get_last_iter_point();
        std::cout << "(x0,y0)=(" << x0 << "," << y0 << "),t=" << t << ", best point at:\n"
                  << solver.get_last_iter_point() << std::endl;
        if (pre_reset && !flag)
            break;
        if (iter_count > 1000 || (iter_count < 1000 && !flag))
        {
            // reset
            t = 0.1;
            x0 = valid_best_point(0, 0);
            y0 = valid_best_point(1, 0);
            flag = true;
            iter_count = 0;
            pre_reset = true;
            points.push_back(std::make_pair(x0,y0));
        }
        else
        {
            pre_reset = false;
        }
        ++iter_count;
        t += 0.1;
    }
    std::cout << "history start points change records:" << std::endl;
    for(std::pair<double,double> vals : points)
    {
        std::cout << vals.first << ", " << vals.second << std::endl;
    }
    return 0;
}
