#include "Mpo.hpp"

Eigen::Tensor<float, 1> Mpo::applyMPO(Eigen::Tensor<float, 2> psi,
                                      Eigen::Tensor<float, 3> left,
                                      Eigen::Tensor<float, 4> firstMpo,
                                      Eigen::Tensor<float, 4> secondMpo,
                                      Eigen::Tensor<float, 3> right) {}

void Mpo::eigenLanczos(Eigen::Tensor<float, 1> psivec,
                       int maxit = 2,
                       int krydim = 4) {}

Mpo::Mpo(Eigen::Tensor<float, 3> mpsTensors[],
         int numberOfmps,
         Eigen::Tensor<float, 3> leftBoundry,
         Eigen::Tensor<float, 4> mpoTensor,
         Eigen::Tensor<float, 3> rightBoundry,
         int maximumDimension,
         int sweepsNumber = 10,
         int displayOn = 2,
         bool updateOn = true,
         int maxIteration = 2,
         int maxKrylovDimension = 4) {
    this->numberOfmps = numberOfmps;
    int chid = mpoTensor.dimensions()[2];

    this->left = new Eigen::Tensor<float, 3>[this->numberOfmps];
    this->left = new Eigen::Tensor<float, 3>[this->numberOfmps];

    this->left[0] = leftBoundry;
    this->right[this->numberOfmps - 1] = rightBoundry;

    for (int i = 0; i < this->numberOfmps - 1; i++) {
        
    }
}

Mpo::~Mpo() {
    delete[] this->left;
    delete[] this->right;
}