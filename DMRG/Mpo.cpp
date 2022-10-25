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
        int chil = mpsTensors[i].dimensions()[0];
        int chir = mpsTensors[i].dimensions()[2];

        std::array<int, 2> dims = {{chil * chid, chir}};
        Eigen::MatrixXf mat = Mpo::tensorToMatrix(mpsTensors[i].reshape(dims));

        Eigen::JacobiSVD<Eigen::MatrixXf, Eigen::ComputeFullU | Eigen::ComputeFullV> svd(mat);

        Eigen::Tensor<float, 1> stemp = Mpo::vectorToTensor(svd.singularValues());
        Eigen::Tensor<float, 2> utemp = Mpo::matrixToTensor(svd.matrixU());
        Eigen::Tensor<float, 2> vtemp = Mpo::matrixToTensor(svd.matrixV());

        std::array<int, 3> dims = {{chil, chid, chir}};
        mpsTensors[i] = utemp.reshape(dims);
        mpsTensors[i + 1] = ncon();
        this->left[i + 1] = ncon();
    }

    int chil = mpsTensors[this->numberOfmps - 1].dimensions()[0];
    int chir = mpsTensors[this->numberOfmps - 1].dimensions()[2];
}

Eigen::Map<const Eigen::MatrixXf> Mpo::tensorToMatrix(const Eigen::Tensor<float, 2> tensor) {
    return Eigen::Map<const Eigen::MatrixXf>(tensor.data());
}

Mpo::~Mpo() {
    delete[] this->left;
    delete[] this->right;
}