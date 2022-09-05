#include "Tilde.hpp"

Tilde::Tilde(Indices* handler) {
    this->handler = handler;

    this->occupiedElements = new int[this->handler->getIrrepsNumber()];
    for (int irrep = 0; irrep < this->handler->getIrrepsNumber(); irrep++) {
        this->occupiedElements[irrep] = this->handler->getOccupiedOrbitalsPerIrrep(irrep) + this->handler->getActiveOrbitalsPerIrrep(irrep);
    }

    this->tildeMatrix = new double***[this->handler->getIrrepsNumber()];
    for (int firstIrrep = 0; firstIrrep < this->handler->getIrrepsNumber(); firstIrrep++) {
        this->tildeMatrix[firstIrrep] = new double**[this->handler->getIrrepsNumber()];
        for (int secondIrrep = 0; secondIrrep < this->handler->getIrrepsNumber(); secondIrrep++) {
            const unsigned int firstSizeBlock = this->occupiedElements[firstIrrep] * this->occupiedElements[secondIrrep];
            const unsigned int secondSizeBlock = this->handler->getOrbitalsPerIrrep(firstIrrep) * this->handler->getOrbitalsPerIrrep(secondIrrep);
            this->tildeMatrix[firstIrrep][secondIrrep] = new double*[firstSizeBlock];
            for (unsigned int combinedIrrep = 0; combinedIrrep < firstSizeBlock; combinedIrrep++) {
                this->tildeMatrix[firstIrrep][secondIrrep][combinedIrrep] = new double[secondSizeBlock];
            }
        }
    }
}

Tilde::~Tilde() {
    for (int firstIrrep = 0; firstIrrep < this->handler->getIrrepsNumber(); firstIrrep++) {
        for (int secondIrrep = 0; secondIrrep < this->handler->getIrrepsNumber(); secondIrrep++) {
            const unsigned int sizeBlock = this->occupiedElements[firstIrrep] * this->occupiedElements[secondIrrep];
            for (unsigned int combinedIrrep = 0; combinedIrrep < sizeBlock; combinedIrrep++) {
                delete[] this->tildeMatrix[firstIrrep][secondIrrep][combinedIrrep];
            }
            delete[] this->tildeMatrix[firstIrrep][secondIrrep];
        }
        delete[] this->tildeMatrix[firstIrrep];
    }
    delete[] this->tildeMatrix;
    delete[] this->occupiedElements;
}

void Tilde::clear() {
    for (int firstIrrep = 0; firstIrrep < this->handler->getIrrepsNumber(); firstIrrep++) {
        for (int secondIrrep = 0; secondIrrep < this->handler->getIrrepsNumber(); secondIrrep++) {
            const unsigned int firstSizeBlock = this->occupiedElements[firstIrrep] * this->occupiedElements[secondIrrep];
            const unsigned int secondSizeBlock = this->handler->getOrbitalsPerIrrep(firstIrrep) * this->handler->getOrbitalsPerIrrep(secondIrrep);
            for (unsigned int thirdIrrep = 0; thirdIrrep < firstSizeBlock; thirdIrrep++) {
                for (unsigned int fourthIrrep = 0; fourthIrrep < secondSizeBlock; fourthIrrep++) {
                    this->tildeMatrix[firstIrrep][secondIrrep][thirdIrrep][fourthIrrep] = 0.0;
                }
            }
        }
    }
}

void Tilde::set(const int firstIrrep, const int secondIrrep, const int first, const int second, const int third, const int fourth, const double value) {
    this->tildeMatrix[firstIrrep][secondIrrep][first + this->occupiedElements[firstIrrep] * third][second + this->handler->getOrbitalsPerIrrep(firstIrrep) * fourth] = value;
}

double Tilde::get(const int firstIrrep, const int secondIrrep, const int first, const int second, const int third, const int fourth) const {
    return this->tildeMatrix[firstIrrep][secondIrrep][first + this->occupiedElements[firstIrrep] * third][second + this->handler->getOrbitalsPerIrrep(firstIrrep) * fourth];
}

double* Tilde::getBlock(const int firstIrrep, const int secondIrrep, const int first, const int third) {
    return this->tildeMatrix[firstIrrep][secondIrrep][first + this->occupiedElements[firstIrrep] * third];
}
