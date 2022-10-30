#pragma once

#include <array>
#include <eigen3/Eigen/SVD>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>
#include <vector>

class Mpo {
   private:
    Eigen::Tensor<float, 3>* left;
    Eigen::Tensor<float, 3>* right;
    Eigen::Tensor<float, 4> mpoTensor;
    double enTemp;
    int numberOfmps;

    Eigen::Tensor<float, 1> applyMPO(Eigen::Tensor<float, 2> psi,
                                     Eigen::Tensor<float, 3> left,
                                     Eigen::Tensor<float, 4> firstMpo,
                                     Eigen::Tensor<float, 4> secondMpo,
                                     Eigen::Tensor<float, 3> right);
    Eigen::Tensor<float, 1> eigenLanczos(Eigen::Tensor<float, 1> psivec,
                      int index,
                      int maxit = 2,
                      int krydim = 4);
    static Eigen::Map<const Eigen::MatrixXf> tensorToMatrix(const Eigen::Tensor<float, 2> tensor);
    static Eigen::TensorMap<const Eigen::Tensor<float, 2>> matrixToTensor(const Eigen::MatrixXf mat);
    static Eigen::TensorMap<const Eigen::Tensor<float, 1>> vectorToTensor(const Eigen::VectorXf vec);

   public:
    Mpo(Eigen::Tensor<float, 3> mpsTensors[],
        int numberOfmps,
        Eigen::Tensor<float, 3> leftBoundry,
        Eigen::Tensor<float, 4> mpoTensor,
        Eigen::Tensor<float, 3> rightBoundry,
        int maximumDimension,
        int sweepsNumber = 10,
        int displayOn = 2,
        bool updateOn = true,
        int maxIteration = 2,
        int maxKrylovDimension = 4);
    virtual ~Mpo();
};