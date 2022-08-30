#pragma once

#include "Indices.hpp"

class Integrals {
   private:
    int irrepsNumber;
    int* coreOrbitalsPerIrrep;
    int* virtualOrbitalsPerIrrep;
    int* totalOrbitalsPerIrrep;
    long long**** coulombPointer;
    long long coulombSize;
    double* coulombArray;
    long long calcNumCoulomb(const bool allocate);
    long long getCoulombPointer(const int firstIrrepCore,
                                const int secondIrrepCore,
                                const int firstIrrepVirtualCore,
                                const int secondIrrepVirtualCore,
                                const int firstCore,
                                const int secondCore,
                                const int firstVirtualCore,
                                const int secondVirtualCore) const;
    long long**** exchangePointer;
    long long exchangeSize;
    double* exchangeArray;
    long long calcNumExchange(const bool allocate);
    long long getExchangePointer(const int firstIrrepCore,
                                 const int secondIrrepCore,
                                 const int firstIrrepVirtual,
                                 const int secondIrrepVirtual,
                                 const int firstCore,
                                 const int secondCore,
                                 const int firstVirtual,
                                 const int secondVirtual) const;

   public:
    Integrals(Indices* indicesHandler);
    virtual ~Integrals();
    void clear();
    void setCoulomb(const int firstIrrepCore,
                    const int secondIrrepCore,
                    const int firstIrrepVirtualCore,
                    const int secondIrrepVirtualCore,
                    const int firstCore,
                    const int secondCore,
                    const int firstVirtualCore,
                    const int secondVirtualCore,
                    const double value);
    void addCoulomb(const int firstIrrepCore,
                     const int secondIrrepCore,
                     const int firstIrrepVirtualCore,
                     const int secondIrrepVirtualCore,
                     const int firstCore,
                     const int secondCore,
                     const int firstVirtualCore,
                     const int secondVirtualCore,
                     const double value);
    double getCoulomb(const int firstIrrepCore,
                       const int secondIrrepCore,
                       const int firstIrrepVirtualCore,
                       const int secondIrrepVirtualCore,
                       const int firstCore,
                       const int secondCore,
                       const int firstVirtualCore,
                       const int secondVirtualCore) const;
    void setExchange(const int firstIrrepCore,
                      const int secondIrrepCore,
                      const int firstIrrepVirtual,
                      const int secondIrrepVirtual,
                      const int firstCore,
                      const int secondCore,
                      const int firstVirtual,
                      const int secondVirtual,
                      const double value);
    void addExchange(const int firstIrrepCore,
                      const int secondIrrepCore,
                      const int firstIrrepVirtual,
                      const int secondIrrepVirtual,
                      const int firstCore,
                      const int secondCore,
                      const int firstVirtual,
                      const int secondVirtual,
                      const double value);
    double getExchange(const int firstIrrepCore,
                        const int secondIrrepCore,
                        const int firstIrrepVirtual,
                        const int secondIrrepVirtual,
                        const int firstCore,
                        const int secondCore,
                        const int firstVirtual,
                        const int secondVirtual) const;
    double FourIndexAPI(const int firstIrrep,
                        const int secondIrrep,
                        const int thirdIrrep,
                        const int fourthIrrep,
                        const int firstIndex,
                        const int secondIndex,
                        const int thirdIndex,
                        const int fourthIndex) const;
};