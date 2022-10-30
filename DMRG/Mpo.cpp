#include "Mpo.hpp"

Eigen::Tensor<float, 1> Mpo::applyMPO(Eigen::Tensor<float, 2> psi,
                                      Eigen::Tensor<float, 3> left,
                                      Eigen::Tensor<float, 4> firstMpo,
                                      Eigen::Tensor<float, 4> secondMpo,
                                      Eigen::Tensor<float, 3> right) {}

Eigen::Tensor<float, 1> Mpo::eigenLanczos(Eigen::Tensor<float, 1> psivec,
                                          int index,
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
    this->mpoTensor = mpoTensor;
    int chid = this->mpoTensor.dimensions()[2];

    this->left = new Eigen::Tensor<float, 3>[this->numberOfmps];
    this->left = new Eigen::Tensor<float, 3>[this->numberOfmps];

    this->left[0] = leftBoundry;
    this->right[this->numberOfmps - 1] = rightBoundry;

    for (int i = 0; i < this->numberOfmps - 1; i++) {
        int chil = mpsTensors[i].dimensions()[0];
        int chir = mpsTensors[i].dimensions()[2];

        std::array<int, 2> dims = {{chil * chid, chir}};
        Eigen::MatrixXf mat = Mpo::tensorToMatrix(mpsTensors[i].reshape(dims));

        Eigen::JacobiSVD<Eigen::MatrixXf, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(mat);

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

    std::array<int, 2> dims = {{chil * chid, chir}};
    Eigen::MatrixXf mat = Mpo::tensorToMatrix(mpsTensors[this->numberOfmps - 1].reshape(dims));

    Eigen::JacobiSVD<Eigen::MatrixXf, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(mat);

    Eigen::VectorXf stemp = svd.singularValues();
    Eigen::Tensor<float, 2> utemp = Mpo::matrixToTensor(svd.matrixU());
    Eigen::MatrixXf vtemp = svd.matrixV();

    std::array<int, 3> dims = {{chil, chid, chir}};
    mpsTensors[this->numberOfmps - 1] = utemp.reshape(dims);

    Eigen::MatrixXf* sWeights = new Eigen::MatrixXf[this->numberOfmps + 1];

    sWeights[this->numberOfmps] = stemp.asDiagonal() * vtemp / stemp.norm();

    std::vector<double> ekeep;
    Eigen::Tensor<float, 3>* b = new Eigen::Tensor<float, 3>[this->numberOfmps];

    for (int i = 1; i < sweepsNumber + 2; i++) {
        if (i == sweepsNumber + 1) {
            updateOn = false;
            displayOn = 0;
        }

        for (int j = this->numberOfmps - 2; j > -1; j--) {
            int chil = mpsTensors[i].dimensions()[0];
            int chir = mpsTensors[i].dimensions()[2];

            Eigen::Tensor<float, 1> psiGround = ncon();
            if (updateOn) {
                psiGround = this->eigenLanczos(psiGround, j, maxIteration, maxKrylovDimension);
                ekeep.push_back(this->enTemp);
            }

            std::array<int, 2> dims = {{chil * chid, chid * chir}};
            Eigen::MatrixXf mat = Mpo::tensorToMatrix(psiGround.reshape(dims));

            Eigen::JacobiSVD<Eigen::MatrixXf, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(mat);

            Eigen::VectorXf stemp = svd.singularValues();
            Eigen::Tensor<float, 2> utemp = Mpo::matrixToTensor(svd.matrixU());
            Eigen::Tensor<float, 2> vtemp = Mpo::matrixToTensor(svd.matrixV());

            int chitemp = std::min((int)stemp.size(), maximumDimension);
            
            std::array<Eigen::Index, 2> offsets = {0, 0};
            std::array<Eigen::Index, 2> extents = {utemp.dimensions()[0], chitemp};
            std::array<int, 3> dims = {{chil, chid, chitemp}};
            mpsTensors[j] = utemp.slice(offsets, extents).reshape(dims);

            sWeights[j + 1] = (stemp.head(chitemp) / stemp.head(chitemp).norm()).asDiagonal();
            
            std::array<Eigen::Index, 2> offsets = {0, 0};
            std::array<Eigen::Index, 2> extents = {chitemp, vtemp.dimensions()[0]};
            std::array<int, 3> dims = {{chitemp, chid, chir}};
            b[j + 1] = vtemp.slice(offsets, extents).reshape(dims);
        }
    }
}

Eigen::Map<const Eigen::MatrixXf> Mpo::tensorToMatrix(const Eigen::Tensor<float, 2> tensor) {
    return Eigen::Map<const Eigen::MatrixXf>(tensor.data());
}

Mpo::~Mpo() {
    delete[] this->left;
    delete[] this->right;
}