#include "Integrals.hpp"

#include <cassert>

Integrals::Integrals(Indices* indicesHandler) {
    this->irrepsNumber = indicesHandler->getIrrepsNumber();
    this->coreOrbitalsPerIrrep = new int[this->irrepsNumber];
    this->virtualOrbitalsPerIrrep = new int[this->irrepsNumber];
    this->totalOrbitalsPerIrrep = new int[this->irrepsNumber];

    for (int irrep = 0; irrep < this->irrepsNumber; irrep++) {
        this->coreOrbitalsPerIrrep[irrep] = indicesHandler->getOccupiedOrbitalsPerIrrep(irrep) + indicesHandler->getActiveOrbitalsPerIrrep(irrep);
        this->virtualOrbitalsPerIrrep[irrep] = indicesHandler->getVirtualOrbitalsPerIrrep(irrep);
        this->totalOrbitalsPerIrrep[irrep] = indicesHandler->getOrbitalsPerIrrep(irrep);
    }

    this->coulombSize = calcNumCoulomb(true);
    this->exchangeSize = calcNumExchange(true);

    this->coulombArray = new double[this->coulombSize];
    this->exchangeArray = new double[this->exchangeSize];
}

Integrals::~Integrals() {
    delete[] this->coulombArray;
    delete[] this->exchangeArray;

    calcNumExchange(false);
    calcNumCoulomb(false);

    delete[] this->coreOrbitalsPerIrrep;
    delete[] this->totalOrbitalsPerIrrep;
    delete[] this->virtualOrbitalsPerIrrep;
}

void Integrals::clear() {
    for (long long counter = 0; counter<this->coulombSize> this->exchangeSize ? this->coulombSize : this->exchangeSize; counter++) {
        if (counter < this->coulombSize) {
            this->coulombArray[counter] = 0.0;
        }

        if (counter < this->exchangeSize) {
            this->exchangeArray[counter] = 0.0;
        }
    }
}

void Integrals::setCoulomb(const int firstIrrepCore,
                           const int secondIrrepCore,
                           const int firstIrrepVirtualCore,
                           const int secondIrrepVirtualCore,
                           const int firstCore,
                           const int secondCore,
                           const int firstVirtualCore,
                           const int secondVirtualCore,
                           const double value) {
    this->coulombArray[this->getCoulombPointer(firstIrrepCore,
                                               secondIrrepCore,
                                               firstIrrepVirtualCore,
                                               secondIrrepVirtualCore,
                                               firstCore,
                                               secondCore,
                                               firstVirtualCore,
                                               secondVirtualCore)] = value;
}

void Integrals::addCoulomb(const int firstIrrepCore,
                           const int secondIrrepCore,
                           const int firstIrrepVirtualCore,
                           const int secondIrrepVirtualCore,
                           const int firstCore,
                           const int secondCore,
                           const int firstVirtualCore,
                           const int secondVirtualCore,
                           const double value) {
    this->coulombArray[this->getCoulombPointer(firstIrrepCore,
                                               secondIrrepCore,
                                               firstIrrepVirtualCore,
                                               secondIrrepVirtualCore,
                                               firstCore,
                                               secondCore,
                                               firstVirtualCore,
                                               secondVirtualCore)] += value;
}

double Integrals::getCoulomb(const int firstIrrepCore,
                             const int secondIrrepCore,
                             const int firstIrrepVirtualCore,
                             const int secondIrrepVirtualCore,
                             const int firstCore,
                             const int secondCore,
                             const int firstVirtualCore,
                             const int secondVirtualCore) const {
    return this->coulombArray[this->getCoulombPointer(firstIrrepCore,
                                                      secondIrrepCore,
                                                      firstIrrepVirtualCore,
                                                      secondIrrepVirtualCore,
                                                      firstCore,
                                                      secondCore,
                                                      firstVirtualCore,
                                                      secondVirtualCore)];
}

void Integrals::setExchange(const int firstIrrepCore,
                            const int secondIrrepCore,
                            const int firstIrrepVirtual,
                            const int secondIrrepVirtual,
                            const int firstCore,
                            const int secondCore,
                            const int firstVirtual,
                            const int secondVirtual,
                            const double value) {
    this->exchangeArray[this->getExchangePointer(firstIrrepCore,
                                                 secondIrrepCore,
                                                 firstIrrepVirtual,
                                                 secondIrrepVirtual,
                                                 firstCore,
                                                 secondCore,
                                                 firstVirtual,
                                                 secondVirtual)] = value;
}

void Integrals::addExchange(const int firstIrrepCore,
                            const int secondIrrepCore,
                            const int firstIrrepVirtual,
                            const int secondIrrepVirtual,
                            const int firstCore,
                            const int secondCore,
                            const int firstVirtual,
                            const int secondVirtual,
                            const double value) {
    this->exchangeArray[this->getExchangePointer(firstIrrepCore,
                                                 secondIrrepCore,
                                                 firstIrrepVirtual,
                                                 secondIrrepVirtual,
                                                 firstCore,
                                                 secondCore,
                                                 firstVirtual,
                                                 secondVirtual)] += value;
}

double Integrals::getExchange(const int firstIrrepCore,
                              const int secondIrrepCore,
                              const int firstIrrepVirtual,
                              const int secondIrrepVirtual,
                              const int firstCore,
                              const int secondCore,
                              const int firstVirtual,
                              const int secondVirtual) const {
    return this->exchangeArray[this->getExchangePointer(firstIrrepCore,
                                                        secondIrrepCore,
                                                        firstIrrepVirtual,
                                                        secondIrrepVirtual,
                                                        firstCore,
                                                        secondCore,
                                                        firstVirtual,
                                                        secondVirtual)];
}

double Integrals::FourIndexAPI(const int firstIrrep,
                               const int secondIrrep,
                               const int thirdIrrep,
                               const int fourthIrrep,
                               const int firstIndex,
                               const int secondIndex,
                               const int thirdIndex,
                               const int fourthIndex) const {
    assert(Irreps::directProduct(firstIrrep, secondIrrep) == Irreps::directProduct(thirdIrrep, fourthIrrep));

    const bool core1 = (firstIndex < this->coreOrbitalsPerIrrep[firstIrrep]) ? true : false;
    const bool core2 = (secondIndex < this->coreOrbitalsPerIrrep[secondIrrep]) ? true : false;
    const bool core3 = (thirdIndex < this->coreOrbitalsPerIrrep[thirdIrrep]) ? true : false;
    const bool core4 = (fourthIndex < this->coreOrbitalsPerIrrep[fourthIrrep]) ? true : false;

    const int numCore = ((core1) ? 1 : 0) + ((core2) ? 1 : 0) + ((core3) ? 1 : 0) + ((core4) ? 1 : 0);
    assert(numCore >= 2);

    if (numCore == 4) {
        return this->getCoulomb(firstIrrep, thirdIrrep, secondIrrep, fourthIrrep, firstIndex, thirdIndex, secondIndex, fourthIndex);
    }

    if (numCore == 3) {
        if ((!core1) || (!core3)) {
            return this->getCoulomb(secondIrrep, fourthIrrep, firstIrrep, thirdIrrep, secondIndex, fourthIndex, firstIndex, thirdIndex);
        }
        if ((!core2) || (!core4)) {
            return this->getCoulomb(firstIrrep, thirdIrrep, secondIrrep, fourthIrrep, firstIndex, thirdIndex, secondIndex, fourthIndex);
        }
    }

    if (numCore == 2) {
        if (!core1) {
            if (!core2) {
                return this->getExchange(thirdIrrep, fourthIrrep, firstIrrep, secondIrrep, thirdIndex, fourthIndex, firstIndex, secondIndex);
            }
            if (!core3) {
                return this->getCoulomb(secondIrrep, fourthIrrep, firstIrrep, thirdIrrep, secondIndex, fourthIndex, firstIndex, thirdIndex);
            }
            if (!core4) {
                return this->getExchange(thirdIrrep, secondIrrep, firstIrrep, fourthIrrep, thirdIndex, secondIndex, firstIndex, fourthIndex);
            }
        }
        if (!core2) {
            if (!core3) {
                return this->getExchange(fourthIrrep, firstIrrep, secondIrrep, thirdIrrep, fourthIndex, firstIndex, secondIndex, thirdIndex);
            }
            if (!core4) {
                return this->getCoulomb(firstIrrep, thirdIrrep, secondIrrep, fourthIrrep, firstIndex, thirdIndex, secondIndex, fourthIndex);
            }
        }
        return this->getExchange(firstIrrep, secondIrrep, thirdIrrep, fourthIrrep, firstIndex, secondIndex, thirdIndex, fourthIndex);
    }

    assert(0 == 1);
    return 0.0;
}

long long Integrals::calcNumCoulomb(const bool allocate) {
    long long objectSize = 0;

    this->coulombPointer = allocate ? new long long***[this->irrepsNumber] : this->coulombPointer;

    for (int coulombCore = 0; coulombCore < this->irrepsNumber; coulombCore++) {
        this->coulombPointer[coulombCore] = allocate ? new long long**[this->irrepsNumber] : this->coulombPointer[coulombCore];

        for (int firstIrrepCore = 0; firstIrrepCore < this->irrepsNumber; firstIrrepCore++) {
            const int secondIrrepCore = Irreps::directProduct(coulombCore, firstIrrepCore);

            if (this->coreOrbitalsPerIrrep[firstIrrepCore] > 0 && this->coreOrbitalsPerIrrep[secondIrrepCore] > 0 && firstIrrepCore <= secondIrrepCore) {
                this->coulombPointer[coulombCore][firstIrrepCore] = allocate ? new long long*[this->irrepsNumber] : this->coulombPointer[coulombCore][firstIrrepCore];

                for (int firstIrrepVirtualCore = 0; firstIrrepVirtualCore < this->irrepsNumber; firstIrrepVirtualCore++) {
                    const int secondIrrepVirtualCore = Irreps::directProduct(coulombCore, firstIrrepVirtualCore);

                    if (this->totalOrbitalsPerIrrep[firstIrrepVirtualCore] > 0 && this->totalOrbitalsPerIrrep[secondIrrepVirtualCore] && firstIrrepVirtualCore <= secondIrrepVirtualCore) {
                        if (coulombCore == 0) {
                            if (allocate) {
                                const long long coreTriangle = (this->coreOrbitalsPerIrrep[firstIrrepCore] * (this->coreOrbitalsPerIrrep[firstIrrepCore] + 1)) / 2;
                                const long long virtualCoreTriangle = (this->totalOrbitalsPerIrrep[firstIrrepVirtualCore] * (this->totalOrbitalsPerIrrep[firstIrrepVirtualCore] + 1)) / 2;

                                this->coulombPointer[coulombCore][firstIrrepCore][firstIrrepVirtualCore] = new long long[coreTriangle];

                                for (int combinedCore = 0; combinedCore < coreTriangle; combinedCore++) {
                                    this->coulombPointer[coulombCore][firstIrrepCore][firstIrrepVirtualCore][combinedCore] = objectSize;
                                    objectSize += virtualCoreTriangle;
                                }
                            } else {
                                delete[] this->coulombPointer[coulombCore][firstIrrepCore][firstIrrepVirtualCore];
                            }
                        } else {
                            if (allocate) {
                                const long long coreSquare = this->coreOrbitalsPerIrrep[firstIrrepCore] * this->coreOrbitalsPerIrrep[secondIrrepCore];
                                const long long virtualCoreSquare = this->totalOrbitalsPerIrrep[firstIrrepVirtualCore] * this->totalOrbitalsPerIrrep[secondIrrepVirtualCore];

                                this->coulombPointer[coulombCore][firstIrrepCore][firstIrrepVirtualCore] = new long long[coreSquare];

                                for (int combinedCore = 0; combinedCore < coreSquare; combinedCore++) {
                                    this->coulombPointer[coulombCore][firstIrrepCore][firstIrrepVirtualCore][combinedCore] = objectSize;
                                    objectSize += virtualCoreSquare;
                                }
                            } else {
                                delete[] this->coulombPointer[coulombCore][firstIrrepCore][firstIrrepVirtualCore];
                            }
                        }
                    }
                }
            }

            if (!allocate) {
                delete[] this->coulombPointer[coulombCore][firstIrrepCore];
            }
        }

        if (!allocate) {
            delete[] this->coulombPointer[coulombCore];
        }
    }

    if (!allocate) {
        delete[] this->coulombPointer;
    }

    return objectSize;
}

long long Integrals::getCoulombPointer(const int firstIrrepCore,
                                       const int secondIrrepCore,
                                       const int firstIrrepVirtualCore,
                                       const int secondIrrepVirtualCore,
                                       const int firstCore,
                                       const int secondCore,
                                       const int firstVirtualCore,
                                       const int secondVirtualCore) const {
    const int coulombCore = Irreps::directProduct(firstIrrepCore, secondIrrepCore);
    assert(coulombCore == Irreps::directProduct(firstIrrepVirtualCore, secondIrrepVirtualCore));

    if (coulombCore == 0) {
        const int IndexCore = (firstCore <= secondCore) ? firstCore + (secondCore * (secondCore + 1)) / 2 : secondCore + (firstCore * (firstCore + 1)) / 2;
        const int IndexVirtualCore = (firstVirtualCore <= secondVirtualCore) ? firstVirtualCore + (secondVirtualCore * (secondVirtualCore + 1)) / 2 : secondVirtualCore + (firstVirtualCore * (firstVirtualCore + 1)) / 2;
        return this->coulombPointer[coulombCore][firstIrrepCore][firstIrrepVirtualCore][IndexCore] + IndexVirtualCore;
    }

    const int irrepCore = (firstIrrepCore < secondIrrepCore) ? firstIrrepCore : secondIrrepCore;
    const int irrepVirtualCore = (firstIrrepVirtualCore < secondIrrepVirtualCore) ? firstIrrepVirtualCore : secondIrrepVirtualCore;
    const int indexCore = (firstIrrepCore < secondIrrepCore) ? firstCore + this->coreOrbitalsPerIrrep[firstIrrepCore] * secondCore : secondCore + this->coreOrbitalsPerIrrep[secondIrrepCore] * firstCore;
    const int indexVirtualCore = (firstIrrepVirtualCore < secondIrrepVirtualCore) ? firstVirtualCore + this->totalOrbitalsPerIrrep[firstIrrepVirtualCore] * secondVirtualCore : secondVirtualCore + this->totalOrbitalsPerIrrep[secondIrrepVirtualCore] * firstVirtualCore;
    return this->coulombPointer[coulombCore][irrepCore][irrepVirtualCore][indexCore] + indexVirtualCore;
}

long long Integrals::calcNumExchange(const bool allocate) {
    long long objectSize = 0;

    this->exchangePointer = allocate ? new long long***[this->irrepsNumber] : this->exchangePointer;

    for (int coulombCore = 0; coulombCore < this->irrepsNumber; coulombCore++) {
        this->exchangePointer[coulombCore] = allocate ? new long long**[this->irrepsNumber] : this->exchangePointer[coulombCore];

        for (int firstIrrepCore = 0; firstIrrepCore < this->irrepsNumber; firstIrrepCore++) {
            const int secondirrepCore = Irreps::directProduct(coulombCore, firstIrrepCore);

            if ((this->coreOrbitalsPerIrrep[firstIrrepCore] > 0) && (this->coreOrbitalsPerIrrep[secondirrepCore] > 0) && (firstIrrepCore <= secondirrepCore)) {
                this->exchangePointer[coulombCore][firstIrrepCore] = allocate ? new long long*[this->irrepsNumber] : this->exchangePointer[coulombCore][firstIrrepCore];

                for (int firstVirtual = 0; firstVirtual < this->irrepsNumber; firstVirtual++) {
                    const int secondVirtualCore = Irreps::directProduct(coulombCore, firstVirtual);

                    if ((this->coreOrbitalsPerIrrep[firstVirtual] > 0) && (this->coreOrbitalsPerIrrep[secondVirtualCore] > 0)) {
                        const long long virtualsquare = this->virtualOrbitalsPerIrrep[firstVirtual] * this->virtualOrbitalsPerIrrep[secondVirtualCore];

                        if (coulombCore == 0) {
                            if (allocate) {
                                const long long coreTriangle = (this->coreOrbitalsPerIrrep[firstIrrepCore] * (this->coreOrbitalsPerIrrep[firstIrrepCore] + 1)) / 2;
                                this->exchangePointer[coulombCore][firstIrrepCore][firstVirtual] = new long long[coreTriangle];

                                for (int combinedCore = 0; combinedCore < coreTriangle; combinedCore++) {
                                    this->exchangePointer[coulombCore][firstIrrepCore][firstVirtual][combinedCore] = objectSize;
                                    objectSize += virtualsquare;
                                }
                            } else {
                                delete[] this->exchangePointer[coulombCore][firstIrrepCore][firstVirtual];
                            }
                        } else {
                            if (allocate) {
                                const long long coreSquare = this->coreOrbitalsPerIrrep[firstIrrepCore] * this->coreOrbitalsPerIrrep[secondirrepCore];
                                this->exchangePointer[coulombCore][firstIrrepCore][firstVirtual] = new long long[coreSquare];

                                for (int combinedCore = 0; combinedCore < coreSquare; combinedCore++) {
                                    this->exchangePointer[coulombCore][firstIrrepCore][firstVirtual][combinedCore] = objectSize;
                                    objectSize += virtualsquare;
                                }
                            } else {
                                delete[] this->exchangePointer[coulombCore][firstIrrepCore][firstVirtual];
                            }
                        }
                    }
                }

                if (!allocate) {
                    delete[] this->exchangePointer[coulombCore][firstIrrepCore];
                }
            }
        }

        if (!allocate) {
            delete[] this->exchangePointer[coulombCore];
        }
    }

    if (!allocate) {
        delete[] this->exchangePointer;
    }

    return objectSize;
}

long long Integrals::getExchangePointer(const int firstIrrepCore,
                                        const int secondIrrepCore,
                                        const int firstIrrepVirtual,
                                        const int secondIrrepVirtual,
                                        const int firstCore,
                                        const int secondCore,
                                        const int firstVirtual,
                                        const int secondVirtual) const {
    const int coulombCore = Irreps::directProduct(firstIrrepCore, secondIrrepCore);
    assert(coulombCore == Irreps::directProduct(firstIrrepVirtual, secondIrrepVirtual));

    if (coulombCore == 0) {
        if (firstCore <= secondCore) {
            return this->exchangePointer[coulombCore][firstIrrepCore][firstIrrepVirtual][firstCore + (secondCore * (secondCore + 1)) / 2] +
                   firstVirtual - this->coreOrbitalsPerIrrep[firstIrrepVirtual] + this->virtualOrbitalsPerIrrep[firstIrrepVirtual] * (secondVirtual - this->coreOrbitalsPerIrrep[secondIrrepVirtual]);
        } else {
            return this->exchangePointer[coulombCore][secondIrrepCore][secondIrrepVirtual][secondCore + (firstCore * (firstCore + 1)) / 2] +
                   secondVirtual - this->coreOrbitalsPerIrrep[secondIrrepVirtual] + this->virtualOrbitalsPerIrrep[secondIrrepVirtual] * (firstVirtual - this->coreOrbitalsPerIrrep[firstIrrepVirtual]);
        }

    } else {
        if (firstIrrepCore < secondIrrepCore) {
            return this->exchangePointer[coulombCore][firstIrrepCore][firstIrrepVirtual][firstCore + this->coreOrbitalsPerIrrep[firstIrrepCore] * secondCore] +
                   firstVirtual - this->coreOrbitalsPerIrrep[firstIrrepVirtual] + this->virtualOrbitalsPerIrrep[firstIrrepVirtual] * (secondVirtual - this->coreOrbitalsPerIrrep[secondIrrepVirtual]);
        } else {
            return this->exchangePointer[coulombCore][secondIrrepCore][secondIrrepVirtual][secondCore + this->coreOrbitalsPerIrrep[secondIrrepCore] * firstCore] +
                   secondVirtual - this->coreOrbitalsPerIrrep[secondIrrepVirtual] + this->virtualOrbitalsPerIrrep[secondIrrepVirtual] * (firstVirtual - this->coreOrbitalsPerIrrep[firstIrrepVirtual]);
        }
    }

    return -1;
}
