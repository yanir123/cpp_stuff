#include "Mpo.hpp"

Eigen::Tensor<float, 1> Mpo::applyMPO(Eigen::Tensor<float, 2> psi,
                                      Eigen::Tensor<float, 3>* left,
                                      Eigen::Tensor<float, 4> firstMpo,
                                      Eigen::Tensor<float, 4> secondMpo,
                                      Eigen::Tensor<float, 3>* right,
                                      int index) {
    return psi.reshape((Eigen::array<int, 4>){left[index].dimensions()[1], firstMpo.dimensions()[3], secondMpo.dimensions()[3], right[index].dimensions()[1]})
                .contract(left[index], (Eigen::array<Eigen::TensorIndexTupleOp<int>, 2>){(1, 3, 5, 7), (2, -1, 2)})
                .contract(firstMpo, (Eigen::array<Eigen::TensorIndexTupleOp<int>, 2>){(1), (2, 4, -2, 3)})
                .contract(secondMpo, (Eigen::array<Eigen::TensorIndexTupleOp<int>, 2>){(1), (4, 6, -3, 5)})
                .contract(right[index], (Eigen::array<Eigen::TensorIndexTupleOp<int>, 2>){(1), (6, -4, 7)})
                .reshape((Eigen::array<int, 1>){left[index].dimensions()[2] * firstMpo.dimensions()[3] * secondMpo.dimensions()[3] * right->dimensions()[2]});
}

Eigen::Tensor<float, 1> Mpo::eigenLanczos(Eigen::Tensor<float, 1> psivec,
                                          int index,
                                          int maxit = 2,
                                          int krydim = 4) {
    Eigen::VectorXf psiVector = Mpo::tensorToVector(psivec);

    if (psiVector.norm() == 0) {
        psiVector = Eigen::VectorXf::Random(psiVector.size());
    }

    Eigen::MatrixXf psi = Eigen::MatrixXf::Zero(psiVector.size(), krydim + 1);
    Eigen::MatrixXf a = Eigen::MatrixXf::Zero(krydim, krydim);
    float dval = 0;

    for (int i = 0; i < maxit; i++) {
        psi(Eigen::all, 0) = psiVector / std::max(psiVector.norm(), (float)1e-16);

        for (int j = 1; j < krydim + 1; j++) {
            psi(Eigen::all, j) = Mpo::tensorToVector(this->applyMPO(Mpo::matrixToTensor(psi(Eigen::all, j - 1)), this->left, this->mpoTensor, this->mpoTensor, this->right, index));

            for (int k = 0; k < j; k++) {
                a(j - 1, k) = psi(Eigen::all, j).dot(psi(Eigen::all, k));
                a(j, j - 1) = a.conjugate()(j - 1, k);
            }

            for (int k = 0; k < j; k++) {
                psi(Eigen::all, j) = psi(Eigen::all, j) - psi(Eigen::all, k).dot(psi(Eigen::all, j)) * psi(Eigen::all, k);
                psi(Eigen::all, j) = psi(Eigen::all, j) / std::max(psi(Eigen::all, j).norm(), (float)1e-16);
            }
        }

        Eigen::EigenSolver<Eigen::MatrixXf> solver(a);

        Eigen::VectorXf dtemp = solver.eigenvalues();
        Eigen::MatrixXf utemp = solver.eigenvectors();

        psiVector = psi(Eigen::all, Eigen::VectorXi::LinSpaced(krydim, 0, krydim - 1)) * utemp(Eigen::all, 0);
        dval = dtemp(0);
    }

    psiVector = psiVector / psiVector.norm();

    this->enTemp = dval;

    return Mpo::vectorToTensor(psiVector);
}

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
    this->mpsTensors = mpsTensors;
    this->numberOfmps = numberOfmps;
    this->mpoTensor = mpoTensor;
    int chid = this->mpoTensor.dimensions()[2];

    this->left = new Eigen::Tensor<float, 3>[this->numberOfmps];
    this->right = new Eigen::Tensor<float, 3>[this->numberOfmps];

    this->left[0] = leftBoundry;
    this->right[this->numberOfmps - 1] = rightBoundry;

    for (int i = 0; i < this->numberOfmps - 1; i++) {
        int chil = this->mpsTensors[i].dimensions()[0];
        int chir = this->mpsTensors[i].dimensions()[2];

        std::array<int, 2> dims = {{chil * chid, chir}};
        Eigen::MatrixXf mat = Mpo::tensorToMatrix(this->mpsTensors[i].reshape(dims));

        Eigen::JacobiSVD<Eigen::MatrixXf, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(mat);

        Eigen::VectorXf stemp = svd.singularValues();
        Eigen::Tensor<float, 2> utemp = Mpo::matrixToTensor(svd.matrixU());
        Eigen::MatrixXf vtemp = svd.matrixV();

        std::array<int, 3> dims = {{chil, chid, chir}};
        this->mpsTensors[i] = utemp.reshape(dims);
        this->mpsTensors[i + 1] = Mpo::matrixToTensor(stemp.asDiagonal() * vtemp).contract(this->mpsTensors[i + 1], (Eigen::array<Eigen::TensorIndexTupleOp<int>, 2>){(-1, 1), (1, -2, -3)}) / stemp.norm();
        this->left[i + 1] = this->left[i].contract(this->mpoTensor, (Eigen::array<Eigen::TensorIndexTupleOp<int>, 2>){(2, 1, 4), (2, -1, 3, 5)});
        this->left[i + 1] = this->left[i + 1].contract(this->mpsTensors[i], (Eigen::array<Eigen::TensorIndexTupleOp<int>, 2>){(1), (4, 5, -3)});
        this->left[i + 1] = this->left[i + 1].contract(this->mpsTensors[i].conjugate(), (Eigen::array<Eigen::TensorIndexTupleOp<int>, 2>){(1), (1, 3, -2)});
    }

    int chil = this->mpsTensors[this->numberOfmps - 1].dimensions()[0];
    int chir = this->mpsTensors[this->numberOfmps - 1].dimensions()[2];

    std::array<int, 2> dims = {{chil * chid, chir}};
    Eigen::MatrixXf mat = Mpo::tensorToMatrix(this->mpsTensors[this->numberOfmps - 1].reshape(dims));

    Eigen::JacobiSVD<Eigen::MatrixXf, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(mat);

    Eigen::VectorXf stemp = svd.singularValues();
    Eigen::Tensor<float, 2> utemp = Mpo::matrixToTensor(svd.matrixU());
    Eigen::MatrixXf vtemp = svd.matrixV();

    std::array<int, 3> dims = {{chil, chid, chir}};
    this->mpsTensors[this->numberOfmps - 1] = utemp.reshape(dims);

    this->sWeights = new Eigen::MatrixXf[this->numberOfmps + 1];

    this->sWeights[this->numberOfmps] = stemp.asDiagonal() * vtemp / stemp.norm();

    this->b = new Eigen::Tensor<float, 3>[this->numberOfmps];

    for (int i = 1; i < sweepsNumber + 2; i++) {
        if (i == sweepsNumber + 1) {
            updateOn = false;
            displayOn = 0;
        }

        for (int j = this->numberOfmps - 2; j > -1; j--) {
            int chil = this->mpsTensors[j].dimensions()[0];
            int chir = this->mpsTensors[j + 1].dimensions()[2];

            Eigen::Tensor<float, 1> psiGround = this->mpsTensors[j]
                                                    .contract(this->mpsTensors[j + 1], (Eigen::array<Eigen::TensorIndexTupleOp<int>, 2>){(-1, -2, 1), (1, -3, 2)})
                                                    .contract(this->sWeights[j + 2], (Eigen::array<Eigen::TensorIndexTupleOp<int>, 2>){(1), (2, -4)})
                                                    .reshape((Eigen::array<int, 1>){chil * chid * chid * chir});

            if (updateOn) {
                psiGround = this->eigenLanczos(psiGround, j, maxIteration, maxKrylovDimension);
                this->ekeep.push_back(this->enTemp);
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
            this->mpsTensors[j] = utemp.slice(offsets, extents).reshape(dims);

            this->sWeights[j + 1] = (stemp.head(chitemp) / stemp.head(chitemp).norm()).asDiagonal();

            std::array<Eigen::Index, 2> offsets = {0, 0};
            std::array<Eigen::Index, 2> extents = {chitemp, vtemp.dimensions()[1]};
            std::array<int, 3> dims = {{chitemp, chid, chir}};
            this->b[j + 1] = vtemp.slice(offsets, extents).reshape(dims);

            this->right[j] = this->mpoTensor
                                 .contract(this->right[j + 1], (Eigen::array<Eigen::TensorIndexTupleOp<int>, 2>){(-1, 2, 3, 5), (2, 1, 4)})
                                 .contract(this->b[j + 1], (Eigen::array<Eigen::TensorIndexTupleOp<int>, 2>){(1), (-3, 5, 4)})
                                 .contract(this->b[j + 1].conjugate(), (Eigen::array<Eigen::TensorIndexTupleOp<int>, 2>){(1), (-2, 3, 1)});

            if (displayOn) {
                std::cout << "Sweep: " << i << " of " << sweepsNumber << ", Loc: " << j << ", Energy: " << this->ekeep[-1] << std::endl;
            }
        }

        chil = this->mpsTensors[0].dimensions()[0];
        chir = this->mpsTensors[0].dimensions()[2];

        Eigen::Tensor<float, 2> Atemp = this->mpsTensors[0]
                                            .contract(this->sWeights[1], (Eigen::array<Eigen::TensorIndexTupleOp<int>, 2>){(-1, -2, 1), (1, -3)})
                                            .reshape((Eigen::array<int, 2>){chil, chid * chir});

        Eigen::MatrixXf mat = Mpo::tensorToMatrix(Atemp);

        Eigen::JacobiSVD<Eigen::MatrixXf, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(mat);

        Eigen::VectorXf stemp = svd.singularValues();
        Eigen::MatrixXf utemp = svd.matrixU();
        Eigen::Tensor<float, 2> vtemp = Mpo::matrixToTensor(svd.matrixV());

        std::array<int, 3> dims = {{chil, chid, chir}};
        this->b[0] = vtemp.reshape(dims);
        this->sWeights[0] = utemp * ((Eigen::MatrixXf)stemp.asDiagonal() / stemp.norm());

        for (int j = 0; j < this->numberOfmps - 1; j++) {
            chil = this->b[j].dimensions()[0];
            chir = this->b[j + 1].dimensions()[2];

            Eigen::Tensor<float, 1> psiGround = Mpo::matrixToTensor(this->sWeights[j])
                                                    .contract(this->b[j], (Eigen::array<Eigen::TensorIndexTupleOp<int>, 2>){(-1, 1), (1, -2, 2)})
                                                    .contract(this->b[j + 1], (Eigen::array<Eigen::TensorIndexTupleOp<int>, 2>){(1), (2, -3, -4)})
                                                    .reshape((Eigen::array<int, 1>){chil * chid * chid * chir});
            if (updateOn) {
                psiGround = eigenLanczos(psiGround, j, maxIteration, maxKrylovDimension);
                this->ekeep.push_back(this->enTemp);
            }

            std::array<int, 2> dims = {{chil * chid, chid * chir}};
            Eigen::MatrixXf mat = Mpo::tensorToMatrix(psiGround.reshape(dims));

            Eigen::JacobiSVD<Eigen::MatrixXf, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(mat);

            Eigen::VectorXf stemp = svd.singularValues();
            Eigen::Tensor<float, 2> utemp = Mpo::matrixToTensor(svd.matrixU());
            Eigen::Tensor<float, 2> vtemp = Mpo::matrixToTensor(svd.matrixV());

            int chitemp = std::min((int)stemp.size(), maximumDimension);

            std::array<int, 3> dims = {{chil, chid, chitemp}};
            std::array<int, 2> offsets = {0, 0};
            std::array<int, 2> extents = {utemp.dimensions()[0], chitemp};
            this->mpsTensors[j] = utemp.slice(offsets, extents).reshape(dims);

            this->sWeights[j + 1] = (stemp.head(chitemp) / stemp.head(chitemp).norm()).asDiagonal();

            std::array<int, 3> dims = {{chitemp, chid, chir}};
            std::array<int, 2> offsets = {0, 0};
            std::array<int, 2> extents = {chitemp, vtemp.dimensions()[1]};
            this->b[j + 1] = vtemp.slice(offsets, extents).reshape(dims);

            this->left[j + 1] = this->left[j]
                                    .contract(this->mpoTensor, (Eigen::array<Eigen::TensorIndexTupleOp<int>, 2>){(2, 1, 4), (2, -1, 3, 5)})
                                    .contract(this->mpsTensors[j], (Eigen::array<Eigen::TensorIndexTupleOp<int>, 2>){(1), (4, 5, -3)})
                                    .contract(this->mpsTensors[j].conjugate(), (Eigen::array<Eigen::TensorIndexTupleOp<int>, 2>){(1), (1, 3, -2)});

            if (displayOn) {
                std::cout << "Sweep: " << i << " of " << sweepsNumber << ", Loc: " << j << ", Energy: " << this->ekeep[-1] << std::endl;
            }
        }

        int chil = this->b[this->numberOfmps - 1].dimensions()[0];
        int chir = this->b[this->numberOfmps - 1].dimensions()[2];

        Eigen::Tensor<float, 2> Atemp = this->b[this->numberOfmps - 1]
                                            .contract(this->sWeights[this->numberOfmps - 1], (Eigen::array<Eigen::TensorIndexTupleOp<int>, 2>){(1, -2, -3), (-1, 1)})
                                            .reshape((Eigen::array<int, 2>){chil * chid, chir});

        Eigen::MatrixXf mat = Mpo::tensorToMatrix(Atemp);

        Eigen::JacobiSVD<Eigen::MatrixXf, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(mat);

        Eigen::VectorXf stemp = svd.singularValues();
        Eigen::Tensor<float, 2> urtemp = Mpo::matrixToTensor(svd.matrixU());
        Eigen::MatrixXf vrtemp = svd.matrixV();

        std::array<int, 3> dims = {{chil, chid, chir}};

        this->mpsTensors[this->numberOfmps - 1] = urtemp.reshape(dims);
        this->sWeights[this->numberOfmps] = (stemp / stemp.norm()) * vrtemp;

        if (displayOn) {
            std::cout << "Sweep: " << i << " of " << sweepsNumber << ", Energy: " << this->ekeep[-1] << ", Bond dim: " << maximumDimension << std::endl;
        }
    }
}

Eigen::Map<const Eigen::MatrixXf> Mpo::tensorToMatrix(const Eigen::Tensor<float, 2> tensor) {
    return Eigen::Map<const Eigen::MatrixXf>(tensor.data());
}

Eigen::Map<const Eigen::VectorXf> tensorToVector(const Eigen::Tensor<float, 1> tensor) {
    return Eigen::Map<const Eigen::VectorXf>(tensor.data());
}

Eigen::TensorMap<const Eigen::Tensor<float, 2>> matrixToTensor(const Eigen::MatrixXf mat) {
    return Eigen::TensorMap<const Eigen::Tensor<float, 2>>(mat.data());
}

Eigen::TensorMap<const Eigen::Tensor<float, 1>> vectorToTensor(const Eigen::VectorXf vec) {
    return Eigen::TensorMap<const Eigen::Tensor<float, 1>>(vec.data());
}

Mpo::~Mpo() {
    delete[] this->left;
    delete[] this->right;
    delete[] this->sWeights;
    delete[] this->b;
}