#include "Indices.hpp"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>

Indices::Indices(const int orbitals, const int group, int* const occupiedOrbitalsPerIrrep, int* const activeOrbitalsPerIrrep, int* const virtualOrbitalsPerIrrep) {
    this->orbitals = orbitals;
    symmInfo.setGroup(group);
    this->irrepsNumber = symmInfo.numberOfIrreps();

    this->orbitalsPerIrrep = new int[this->irrepsNumber];
    this->occupiedOrbitalsPerIrrep = new int[this->irrepsNumber];
    this->activeOrbitalsPerIrrep = new int[this->irrepsNumber];
    this->virtualOrbitalsPerIrrep = new int[this->irrepsNumber];
    this->cumulativeOrbitalsPerIrrep = new int[this->irrepsNumber + 1];
    this->cumulativeActiveOrbitalsPerIrrep = new int[this->irrepsNumber + 1];

    int totalNumberOfOrbs = 0;
    this->cumulativeOrbitalsPerIrrep[0] = 0;
    this->cumulativeActiveOrbitalsPerIrrep[0] = 0;

    for (int irrep = 0; irrep < this->irrepsNumber; irrep++) {
        assert(occupiedOrbitalsPerIrrep[irrep] >= 0);
        assert(activeOrbitalsPerIrrep[irrep] >= 0);
        assert(virtualOrbitalsPerIrrep[irrep] >= 0);

        this->orbitalsPerIrrep[irrep] = occupiedOrbitalsPerIrrep[irrep] + activeOrbitalsPerIrrep[irrep] + virtualOrbitalsPerIrrep[irrep];
        this->occupiedOrbitalsPerIrrep[irrep] = occupiedOrbitalsPerIrrep[irrep];
        this->activeOrbitalsPerIrrep[irrep] = activeOrbitalsPerIrrep[irrep];
        this->virtualOrbitalsPerIrrep[irrep] = virtualOrbitalsPerIrrep[irrep];

        totalNumberOfOrbs += this->orbitalsPerIrrep[irrep];

        this->cumulativeOrbitalsPerIrrep[irrep + 1] = this->cumulativeOrbitalsPerIrrep[irrep] + this->orbitalsPerIrrep[irrep];
        this->cumulativeActiveOrbitalsPerIrrep[irrep + 1] = this->cumulativeActiveOrbitalsPerIrrep[irrep] + this->activeOrbitalsPerIrrep[irrep];
    }

    assert(totalNumberOfOrbs == this->orbitals);

    this->irrepPerDMRG = new int[this->cumulativeActiveOrbitalsPerIrrep[this->irrepsNumber]];
    this->irrepPerOrbital = new int[this->orbitals];

    for (int irrep = 0; irrep < this->irrepsNumber; irrep++) {
        for (int count = 0; count < this->activeOrbitalsPerIrrep[irrep]; count++) {
            this->irrepPerDMRG[this->cumulativeActiveOrbitalsPerIrrep[irrep] + count] = irrep;
        }
        for (int count = 0; count < this->orbitalsPerIrrep[irrep]; count++) {
            this->irrepPerOrbital[this->cumulativeOrbitalsPerIrrep[irrep] + count] = irrep;
        }
    }
}

Indices::~Indices() {
    delete[] this->orbitalsPerIrrep;
    delete[] this->activeOrbitalsPerIrrep;
    delete[] this->occupiedOrbitalsPerIrrep;
    delete[] this->virtualOrbitalsPerIrrep;
    delete[] this->irrepPerDMRG;
    delete[] this->irrepPerOrbital;
    delete[] this->cumulativeActiveOrbitalsPerIrrep;
    delete[] this->cumulativeOrbitalsPerIrrep;
}

int Indices::getOrbitals() const {
    return this->orbitals;
}

int Indices::groupNumber() const {
    return symmInfo.getGroupNumber();
}

int Indices::getIrrepsNumber() const {
    return this->irrepsNumber;
}

int Indices::getOrbitalsPerIrrep(const int irrep) const {
    return this->orbitalsPerIrrep[irrep];
}

int Indices::getOccupiedOrbitalsPerIrrep(const int irrep) const {
    return this->occupiedOrbitalsPerIrrep[irrep];
}

int Indices::getActiveOrbitalsPerIrrep(const int irrep) const {
    return this->activeOrbitalsPerIrrep[irrep];
}

int Indices::getVirtualOrbitalsPerIrrep(const int irrep) const {
    return this->virtualOrbitalsPerIrrep[irrep];
}

int Indices::getCumulativeActiveOrbitalsPerIrrep(const int irrep) const {
    return this->cumulativeActiveOrbitalsPerIrrep[irrep];
}

int Indices::getOriginalOccupiedOrbitalsPerIrrep(const int irrep) const {
    return this->cumulativeOrbitalsPerIrrep[irrep];
}

int Indices::getOriginalActiveOrbitalsPerIrrep(const int irrep) const {
    return this->cumulativeOrbitalsPerIrrep[irrep] + this->occupiedOrbitalsPerIrrep[irrep];
}

int Indices::getOriginalVirtualOrbitalsPerIrrep(const int irrep) const {
    return this->cumulativeOrbitalsPerIrrep[irrep + 1] - this->virtualOrbitalsPerIrrep[irrep];
}

int* Indices::getIrrepPerDMRG() {
    return this->irrepPerDMRG;
}

int Indices::getIrrepPerOrbital(const int index) const {
    return this->irrepPerOrbital[index];
}

int Indices::getOccupiedOrbitalsPerIrrepSum() const {
    int sum = 0;
    for (int irrep = 0; irrep < this->irrepsNumber; irrep++) {
        sum += this->getOccupiedOrbitalsPerIrrep(irrep);
    }
    return sum;
}

int Indices::getActiveOrbitalsPerIrrepSum() const {
    int sum = 0;
    for (int irrep = 0; irrep < this->irrepsNumber; irrep++) {
        sum += this->getActiveOrbitalsPerIrrep(irrep);
    }
    return sum;
}

int Indices::getVirtualOrbitalsPerIrrepSum() const {
    int sum = 0;
    for (int irrep = 0; irrep < this->irrepsNumber; irrep++) {
        sum += this->getVirtualOrbitalsPerIrrep(irrep);
    }
    return sum;
}

void Indices::prnint() const {
    std::cout << "NORB  = [ ";
    for (int irrep = 0; irrep < this->irrepsNumber - 1; irrep++) {
        std::cout << this->orbitalsPerIrrep[irrep] << " , ";
    }
    std::cout << this->orbitalsPerIrrep[this->irrepsNumber - 1] << " ]" << std::endl;

    std::cout << "NOCC  = [ ";
    for (int irrep = 0; irrep < this->irrepsNumber - 1; irrep++) {
        std::cout << this->occupiedOrbitalsPerIrrep[irrep] << " , ";
    }
    std::cout << this->occupiedOrbitalsPerIrrep[this->irrepsNumber - 1] << " ]" << std::endl;

    std::cout << "NDMRG = [ ";
    for (int irrep = 0; irrep < this->irrepsNumber - 1; irrep++) {
        std::cout << this->activeOrbitalsPerIrrep[irrep] << " , ";
    }
    std::cout << this->activeOrbitalsPerIrrep[this->irrepsNumber - 1] << " ]" << std::endl;

    std::cout << "NVIRT = [ ";
    for (int irrep = 0; irrep < this->irrepsNumber - 1; irrep++) {
        std::cout << this->virtualOrbitalsPerIrrep[irrep] << " , ";
    }
    std::cout << this->virtualOrbitalsPerIrrep[this->irrepsNumber - 1] << " ]" << std::endl;
}