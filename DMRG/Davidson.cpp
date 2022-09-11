#include "Davidson.hpp"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>

#include "Lapack.hpp"

Davidson::Davidson(const int vectorLength, const int maxVectorsNumber, const int keepVectorsNumber, const double residualTolerence, const double diagCutoff, const bool debugPrint, const char problemType = 'E') {
    assert((problemType == 'E') || (problemType == 'L'));

    this->debugPrint = debugPrint;
    this->vectorLength = vectorLength;
    this->maxVectorsNumber = maxVectorsNumber;
    this->keepVectorsNumber = keepVectorsNumber;
    this->problemType = problemType;
    this->diagCutoff = diagCutoff;
    this->residualTolerence = residualTolerence;

    this->state = 'I';
    this->multipicationsNumber = 0;

    this->vectorsNumber = 0;
    this->vectors = new double*[maxVectorsNumber];
    this->hVectors = new double*[maxVectorsNumber];
    this->numAllocated = 0;

    this->mXMatrix = new double[maxVectorsNumber * maxVectorsNumber];
    this->mXMatrixEigens = new double[maxVectorsNumber];
    this->mXMatrixVectors = new double[maxVectorsNumber];
    this->mXMatrixWorkLength = 3 * maxVectorsNumber - 1;
    this->mXMatrixWork = new double[this->mXMatrixWorkLength];
    this->mXMatrixRhs = problemType == 'L' ? new double[maxVectorsNumber] : nullptr;

    this->diag = new double[this->vectorLength];
    this->tVector = new double[this->vectorLength];
    this->uVector = new double[this->vectorLength];
    this->workVector = new double[this->vectorLength];
    this->rhs = problemType == 'L' ? new double[this->vectorLength] : nullptr;

    this->ReorthoLowdin = nullptr;
    this->ReorthoEigensVectors = nullptr;
    this->ReorthoOverlap = nullptr;
    this->ReorthoOverlapEigens = nullptr;
}

Davidson::~Davidson() {
    for (int i = 0; i < this->numAllocated; i++) {
        delete[] this->vectors[i];
        delete[] this->hVectors[i];
    }

    delete[] this->vectors;
    delete[] this->hVectors;

    delete[] this->mXMatrix;
    delete[] this->mXMatrixEigens;
    delete[] this->mXMatrixVectors;
    delete[] this->mXMatrixWork;

    if (this->mXMatrixRhs != nullptr) {
        delete[] this->mXMatrixRhs;
    }

    delete[] this->diag;
    delete[] this->tVector;
    delete[] this->uVector;
    delete[] this->workVector;

    if (this->rhs != nullptr) {
        delete[] this->rhs;
    }

    if (this->ReorthoEigensVectors != nullptr) {
        delete[] this->ReorthoEigensVectors;
    }

    if (this->ReorthoLowdin != nullptr) {
        delete[] this->ReorthoLowdin;
    }

    if (this->ReorthoOverlap != nullptr) {
        delete[] this->ReorthoOverlap;
    }

    if (this->ReorthoOverlapEigens != nullptr) {
        delete[] this->ReorthoOverlapEigens;
    }
}

char Davidson::fetchInstruction(double** const pointers) {
    if (this->state == 'I') {
        pointers[0] = this->tVector;
        pointers[1] = this->uVector;

        if (problemType == 'L') {
            pointers[2] = this->rhs;
        }
        this->state = 'U';
        return 'A';
    }

    if (state == 'U') {
        this->SafetyCheckGuess();
        this->AddNewVector();
        pointers[0] = this->vectors[this->vectorsNumber];
        pointers[1] = this->hVectors[this->vectorsNumber];
        this->multipicationsNumber++;
        this->state = 'N';
        return 'B';
    }

    if (this->state == 'N') {
        const double rnorm = this->DiagonalizeSmallMatrixAndCalcResidual();

        if (rnorm > this->residualTolerence) {
            this->CalculateNewVector();
            if (this->vectorsNumber == this->maxVectorsNumber) {
                this->Deflation();
                pointers[0] = this->vectors[this->vectorsNumber];
                pointers[1] = this->hVectors[this->vectorsNumber];
                this->multipicationsNumber++;
                this->vectorsNumber++;
                this->state = 'F';
                return 'B';
            }

            this->AddNewVector();
            pointers[0] = this->vectors[this->vectorsNumber];
            pointers[1] = this->hVectors[this->vectorsNumber];
            this->multipicationsNumber++;
            this->state = 'N';
            return 'B';
        }

        this->state = 'C';
        pointers[0] = this->uVector;
        pointers[1] = this->workVector;

        if (this->problemType == 'E') {
            this->workVector[0] = this->mXMatrixEigens[0];
        }

        if (this->problemType == 'L') {
            this->workVector[0] = rnorm;
        }

        return 'C';
    }

    if (this->state == 'F') {
        if (this->vectorsNumber == this->keepVectorsNumber) {
            this->MxMafterDeflation();
            this->AddNewVector();

            pointers[0] = this->vectors[this->vectorsNumber];
            pointers[1] = this->hVectors[this->vectorsNumber];
            this->multipicationsNumber++;
            this->state = 'N';
            return 'B';
        }

        pointers[0] = this->vectors[this->vectorsNumber];
        pointers[1] = this->hVectors[this->vectorsNumber];
        this->multipicationsNumber++;
        this->vectorsNumber++;
        this->state = 'F';
        return 'B';
    }

    return 'D';
}

int Davidson::getNumMultiplications() const {
    return this->multipicationsNumber;
}

double Davidson::FrobeniusNorm(double* currentVector) {
    char frobenius = 'F';
    int increment = 1;
    return dlange_(&frobenius, &this->vectorLength, &increment, currentVector, &this->vectorLength, nullptr);
}

void Davidson::SafetyCheckGuess() {
    const double norm = this->FrobeniusNorm(this->tVector);

    if (norm == 0.0) {
        for (int i = 0; i < this->vectorLength; i++) {
            this->tVector[i] = ((double)rand()) / RAND_MAX;
        }

        if (this->debugPrint) {
            std::cout << "WARNING AT DAVIDSON : Initial guess was a zero-vector. Now it is overwritten with random numbers." << std::endl;
        }
    }
}

void Davidson::AddNewVector() {
    int increment = 1;

    for (int i = 0; i < this->vectorsNumber; i++) {
        double minusOverlap = -ddot_(&this->vectorLength, this->tVector, &increment, this->vectors[i], &increment);
        daxpy_(&this->vectorLength, &minusOverlap, this->vectors[i], &increment, tVector, &increment);
    }

    double alpha = 1.0 / FrobeniusNorm(this->tVector);
    dscal_(&this->vectorLength, &alpha, this->tVector, &increment);

    if (this->vectorsNumber < this->numAllocated) {
        double* temp = this->vectors[this->vectorsNumber];
        this->vectors[this->vectorsNumber] = this->tVector;
        this->tVector = temp;
    } else {
        this->vectors[this->numAllocated] = this->tVector;
        this->hVectors[this->numAllocated] = new double[this->vectorLength];
        this->numAllocated++;
    }
}

double Davidson::DiagonalizeSmallMatrixAndCalcResidual() {
    int increment = 1;

    if (problemType == 'E') {
        for (int i = 0; i < this->vectorsNumber; i++) {
            this->mXMatrix[i + this->maxVectorsNumber * this->vectorsNumber] = ddot_(&this->vectorLength, this->vectors[this->vectorsNumber], &increment, this->hVectors[i], &increment);
            this->mXMatrix[this->vectorsNumber + this->maxVectorsNumber * i] = this->mXMatrix[i + this->maxVectorsNumber * this->vectorsNumber];
        }

        this->mXMatrix[this->vectorsNumber + this->maxVectorsNumber * this->vectorsNumber] = ddot_(&this->vectorLength, this->vectors[this->vectorsNumber], &increment, this->hVectors[this->vectorsNumber], &increment);
    } else {
        for (int i = 0; i < this->vectorsNumber; i++) {
            this->mXMatrix[i + this->maxVectorsNumber * this->vectorsNumber] = ddot_(&this->vectorLength, this->hVectors[this->vectorsNumber], &increment, this->hVectors[i], &increment);
            this->mXMatrix[this->vectorsNumber + this->maxVectorsNumber * i] = this->mXMatrix[i + this->maxVectorsNumber * this->vectorsNumber];
        }

        this->mXMatrix[this->vectorsNumber + this->maxVectorsNumber * this->vectorsNumber] = ddot_(&this->vectorLength, this->hVectors[this->vectorsNumber], &increment, this->hVectors[this->vectorsNumber], &increment);
        this->mXMatrixRhs[this->vectorsNumber] = ddot_(&vectorLength, this->hVectors[this->vectorsNumber], &increment, this->rhs, &increment);
    }

    this->vectorsNumber++;

    char jobz = 'V';
    char uplo = 'U';
    int info;

    for (int i = 0; i < this->vectorsNumber; i++) {
        for (int j = 0; j < this->vectorsNumber; j++) {
            mXMatrixVectors[i + this->maxVectorsNumber * j] = this->mXMatrix[i + this->maxVectorsNumber * j];
        }
    }

    dsyev_(&jobz, &uplo, &this->vectorsNumber, this->mXMatrixVectors, &this->maxVectorsNumber, this->mXMatrixEigens, this->mXMatrixWork, &this->mXMatrixWorkLength, &info);

    if (this->problemType == 'L') {
        double one = 1.0;
        double set = 0.0;
        char trans = 'T';
        char notra = 'N';
        dgemm_(&trans, &notra, &this->vectorsNumber, &increment, &this->vectorsNumber, &one, this->mXMatrixVectors, &this->maxVectorsNumber, this->mXMatrixRhs, &this->maxVectorsNumber, &set, this->mXMatrixWork, &this->maxVectorsNumber);

        for (int i = 0; i < this->vectorsNumber; i++) {
            double currentEigenValue = this->mXMatrixEigens[i];
            if (fabs(currentEigenValue) < this->diagCutoff) {
                currentEigenValue = this->diagCutoff * (currentEigenValue < 0 ? -1 : 1);
                if (this->debugPrint) {
                    std::cout << "WARNING AT DAVIDSON : The eigenvalue " << this->mXMatrixEigens[i] << " to solve Ax = b has been overwritten with " << currentEigenValue << "." << std::endl;
                }
            }

            this->mXMatrixWork[i] = this->mXMatrixWork[i] / currentEigenValue;
        }

        dgemm_(&notra, &notra, &this->vectorsNumber, &increment, &this->vectorsNumber, &one, this->mXMatrixVectors, &this->maxVectorsNumber, this->mXMatrixWork, &this->vectorsNumber, &set, this->mXMatrixWork + this->maxVectorsNumber, &this->maxVectorsNumber);

        for (int i = 0; i < this->vectorsNumber; i++) {
            this->mXMatrixVectors[i] = this->mXMatrixWork[this->maxVectorsNumber + i];
        }
    }

    for (int i = 0; i < this->vectorLength; i++) {
        this->tVector[i] = 0.0;
        this->uVector[i] = 0.0;
    }

    for (int i = 0; i < this->vectorsNumber; i++) {
        double alpha = this->mXMatrixVectors[i];
        daxpy_(&this->vectorLength, &alpha, this->hVectors[i], &increment, this->tVector, &increment);
        daxpy_(&this->vectorLength, &alpha, this->vectors[i], &increment, this->uVector, &increment);
    }

    if (this->problemType == 'E') {
        double alpha = -this->mXMatrixEigens[0];
        daxpy_(&this->vectorLength, &alpha, this->uVector, &increment, this->tVector, &increment);
    } else {
        double alpha = -1.0;
        daxpy_(&this->vectorLength, &alpha, this->rhs, &increment, this->tVector, &increment);
    }

    return this->FrobeniusNorm(this->tVector);
}

void Davidson::CalculateNewVector() {
    int increment = 1;
    const double shift = this->problemType == 'E' ? this->mXMatrixEigens[0] : 0.0;

    for (int i = 0; i < this->vectorLength; i++) {
        const double difference = this->diag[i] - shift;
        const double fabsDifference = fabs(difference);

        if (fabsDifference > this->diagCutoff) {
            this->workVector[i] = this->uVector[i] / difference;
        } else {
            this->workVector[i] = this->uVector[i] / this->diagCutoff;
        }

        if (this->debugPrint) {
            std::cout << "WARNING AT DAVIDSON : fabs( precon[" << i << "] ) = " << fabsDifference << std::endl;
        }
    }

    double alpha = ddot_(&this->vectorLength, this->workVector, &increment, this->tVector, &increment) / ddot_(&this->vectorLength, this->workVector, &increment, this->uVector, &increment);
    daxpy_(&this->vectorLength, &alpha, this->uVector, &increment, this->tVector, &increment);

    for (int i = 0; i < this->vectorLength; i++) {
        const double difference = this->diag[i] - shift;
        const double fabsDifference = fabs(difference);

        if (fabsDifference > this->diagCutoff) {
            this->tVector[i] = -this->tVector[i] / difference;
        } else {
            this->tVector[i] = -this->tVector[i] / this->diagCutoff;
        }
    }
}

void Davidson::Deflation() {
    int increment = 1;

    if (this->keepVectorsNumber <= 1) {
        double alpha = 1.0 / this->FrobeniusNorm(this->uVector);
        dscal_(&this->vectorLength, &alpha, this->uVector, &increment);
        dcopy_(&this->vectorLength, uVector, &increment, this->vectors[0], &increment);
    } else {
        if (this->problemType == 'L') {
            this->SolveLinearSystemDeflation(this->keepVectorsNumber);
        }

        this->ReorthoEigensVectors = this->ReorthoEigensVectors == nullptr ? new double[this->vectorLength * this->keepVectorsNumber] : this->ReorthoEigensVectors;
        this->ReorthoOverlap = this->ReorthoOverlap == nullptr ? new double[this->keepVectorsNumber * this->keepVectorsNumber] : this->ReorthoOverlap;
        this->ReorthoOverlapEigens = this->ReorthoOverlapEigens == nullptr ? new double[this->keepVectorsNumber] : this->ReorthoOverlapEigens;
        this->ReorthoLowdin = this->ReorthoLowdin == nullptr ? new double[this->keepVectorsNumber * this->keepVectorsNumber] : this->ReorthoLowdin;

        dcopy_(&this->vectorLength, this->uVector, &increment, this->ReorthoEigensVectors, &increment);
        for (int i = 0; i < this->keepVectorsNumber; i++) {
            for (int j = 0; j < this->vectorLength; j++) {
                this->ReorthoEigensVectors[j + this->vectorLength * i] = 0.0;
                for (int k = 0; k < this->vectorLength; k++) {
                    this->ReorthoEigensVectors[j + this->vectorLength * i] += this->vectors[k][j] * this->mXMatrixVectors[k + this->maxVectorsNumber * i];
                }
            }
        }

        char trans = 'T';
        char notr = 'N';
        double one = 1.0;
        double zero = 0.0;  //set
        dgemm_(&trans, &notr, &this->keepVectorsNumber, &this->keepVectorsNumber, &this->vectorLength, &one, this->ReorthoEigensVectors, &this->vectorLength, this->ReorthoEigensVectors, &this->vectorLength, &zero, this->ReorthoOverlap, &this->keepVectorsNumber);

        char jobz = 'V';
        char uplo = 'U';
        int info;
        dsyev_(&jobz, &uplo, &this->keepVectorsNumber, this->ReorthoOverlap, &this->keepVectorsNumber, &this->ReorthoOverlapEigens, this->mXMatrixWork, &this->mXMatrixWorkLength, &info);
        for (int i = 0; i < this->keepVectorsNumber; i++) {
            this->ReorthoOverlapEigens[i] = pow(this->ReorthoOverlapEigens[i], -0.25);
            dscal_(&this->keepVectorsNumber, this->ReorthoOverlapEigens + i, this->ReorthoOverlap + this->keepVectorsNumber * i, &increment);
        }

        dgemm_(&notr, &trans, &this->keepVectorsNumber, &this->keepVectorsNumber, &this->keepVectorsNumber, &one, this->ReorthoOverlap, &this->keepVectorsNumber, this->ReorthoOverlap, &this->keepVectorsNumber, &zero, this->ReorthoLowdin, &this->keepVectorsNumber);

        for (int i = 0; i < this->keepVectorsNumber; i++) {
            for (int j = 0; j < this->vectorLength; j++) {
                this->vectors[i][j] = 0.0;
            }

            for (int j = 0; j < this->keepVectorsNumber; j++) {
                daxpy_(&this->vectorLength, this->ReorthoLowdin + j + this->keepVectorsNumber * i, this->ReorthoEigensVectors + this->vectorLength * j, &increment, this->vectors[i], &increment);
            }
        }
    }

    this->vectorsNumber = 0;
}

void Davidson::MxMafterDeflation() {
    int increment = 1;

    if (this->problemType == 'E') {
        for (int i = 0; i < this->keepVectorsNumber; i++) {
            for (int j = i; j < this->keepVectorsNumber; j++) {
                this->mXMatrix[i + this->maxVectorsNumber * j] = ddot_(&this->vectorLength, this->vectors[i], &increment, this->hVectors[j], &increment);
                this->mXMatrix[j + this->maxVectorsNumber * i] = this->mXMatrix[i + this->maxVectorsNumber * j];
            }
        }
    } else {
        for (int i = 0; i < this->keepVectorsNumber; i++) {
            for (int j = i; j < this->keepVectorsNumber; j++) {
                this->mXMatrix[i + this->maxVectorsNumber * j] = ddot_(&this->vectorLength, this->hVectors[i], &increment, this->hVectors[j], &increment);
                this->mXMatrix[j + this->maxVectorsNumber * i] = this->mXMatrix[i + this->maxVectorsNumber * j];
            }
        }
        for (int i = 0; i < this->keepVectorsNumber; i++) {
            this->mXMatrixRhs[i] = ddot_(&this->vectorLength, this->hVectors[i], &increment, this->rhs, &increment);
        }
    }
}

void Davidson::SolveLinearSystemDeflation(const int SOLUTIONS_NUMBER) {
    assert(this->problemType == 'L');
    assert(this->vectorsNumber == this->maxVectorsNumber);
    assert(SOLUTIONS_NUMBER <= this->maxVectorsNumber);
    assert(SOLUTIONS_NUMBER >= 2);

    double* firstWork = new double[this->maxVectorsNumber * this->maxVectorsNumber];
    double* secondWork = new double[this->maxVectorsNumber * this->maxVectorsNumber];
    double* thirdWork = new double[this->maxVectorsNumber * SOLUTIONS_NUMBER];

    for (int solution = 0; solution < SOLUTIONS_NUMBER; solution++) {
        for (int i = 0; i < this->maxVectorsNumber * this->maxVectorsNumber; i++) {
            firstWork[i] = 0.0;
        }

        for (int diag = 0; diag < this->maxVectorsNumber; diag++) {
            firstWork[diag * (1 + this->maxVectorsNumber)] = 1.0;
        }

        for (int prev = 0; prev < solution; prev++) {
            for (int i = 0; i < this->maxVectorsNumber; i++) {
                for (int j = 0; j < this->maxVectorsNumber; j++) {
                    firstWork[i + this->maxVectorsNumber * j] -= secondWork[i + this->maxVectorsNumber * prev] * secondWork[j + this->maxVectorsNumber * prev];
                }
            }
        }

        {
            double one = 1.0;
            double set = 0.0;
            char notrans = 'N';
            int inc1 = 1;
            dgemm_(&notrans, &notrans, &this->maxVectorsNumber, &this->maxVectorsNumber, &this->maxVectorsNumber, &one, firstWork, &this->maxVectorsNumber, this->mXMatrix, &this->maxVectorsNumber, &set, thirdWork, &this->maxVectorsNumber);
            dgemm_(&notrans, &notrans, &this->maxVectorsNumber, &this->maxVectorsNumber, &this->maxVectorsNumber, &one, thirdWork, &this->maxVectorsNumber, firstWork, &this->maxVectorsNumber, &set, this->mXMatrixVectors, &this->maxVectorsNumber);
            dgemm_(&notrans, &notrans, &this->maxVectorsNumber, &inc1, &this->maxVectorsNumber, &one, firstWork, &this->maxVectorsNumber, this->mXMatrixRhs, &this->maxVectorsNumber, &set, thirdWork, &this->maxVectorsNumber);
        }

        {
            char jobz = 'V';
            char uplo = 'U';
            int info;
            dsyev_(&jobz, &uplo, &this->maxVectorsNumber, this->mXMatrixVectors, &this->maxVectorsNumber, this->mXMatrixEigens, this->mXMatrixWork, &this->mXMatrixWorkLength, &info);
        }

        {
            double one = 1.0;
            double set = 0.0;
            char trans = 'T';
            char notrans = 'N';
            int inc1 = 1;
            dgemm_(&trans, &notrans, &this->maxVectorsNumber, &inc1, &this->maxVectorsNumber, &one, this->mXMatrixVectors, &this->maxVectorsNumber, thirdWork, &this->maxVectorsNumber, &set, this->mXMatrixWork, &this->maxVectorsNumber);
            for (int diag = 0; diag < this->maxVectorsNumber; diag++) {
                if (diag < solution) {
                    this->mXMatrixWork[diag] = 0.0;  // PSEUDOINVERSE
                } else {
                    double current_eigenvalue = this->mXMatrixEigens[diag];
                    if (fabs(current_eigenvalue) < this->diagCutoff) {
                        current_eigenvalue = this->diagCutoff * ((current_eigenvalue < 0.0) ? -1 : 1);
                        if (this->debugPrint) {
                            std::cout << "WARNING AT DAVIDSON : The eigenvalue " << this->mXMatrixEigens[diag] << " to solve Ax = b has been overwritten with " << current_eigenvalue << "." << std::endl;
                        }
                    }
                    this->mXMatrixWork[diag] = this->mXMatrixWork[diag] / current_eigenvalue;
                }
            }
            dgemm_(&notrans, &notrans, &this->maxVectorsNumber, &inc1, &this->maxVectorsNumber, &one, this->mXMatrixVectors, &this->maxVectorsNumber, this->mXMatrixWork, &this->maxVectorsNumber, &set, secondWork + this->maxVectorsNumber * solution, &this->maxVectorsNumber);
        }

        {
            int inc1 = 1;
            double* ptr = secondWork + this->maxVectorsNumber * solution;
            const double twonorm = sqrt(ddot_(&this->maxVectorsNumber, ptr, &inc1, ptr, &inc1));
            double factor = 1.0 / twonorm;
            dscal_(&this->maxVectorsNumber, &factor, ptr, &inc1);
        }

        {
            int increment = 1;
            int size = this->maxVectorsNumber * SOLUTIONS_NUMBER;
            dcopy_(&size, secondWork, &increment, this->mXMatrixVectors, &increment);
        }
    }

    delete[] firstWork;
    delete[] secondWork;
    delete[] thirdWork;
}
