#pragma once

#include "Irreps.hpp"

class Indices {
   private:
    int orbitals;
    Irreps symmInfo;
    int irrepsNumber;
    int* orbitalsPerIrrep;
    int* occupiedOrbitalsPerIrrep;
    int* activeOrbitalsPerIrrep;
    int* virtualOrbitalsPerIrrep;
    int* cumulativeOrbitalsPerIrrep;
    int* cumulativeActiveOrbitalsPerIrrep;
    int* irrepPerDMRG;
    int* irrepPerOrbital;
   public:
    Indices(const int orbitals, const int group, int* const occupiedOrbitalsPerIrrep, int* const activeOrbitalsPerIrrep, int* const virtualOrbitalsPerIrrep);
    virtual ~Indices();
    int getOrbitals() const;
    int groupNumber() const;
    int getIrrepsNumber() const;
    int getOrbitalsPerIrrep(const int irrep) const;
    int getOccupiedOrbitalsPerIrrep(const int irrep) const;
    int getActiveOrbitalsPerIrrep(const int irrep) const;
    int getVirtualOrbitalsPerIrrep(const int irrep) const;
    int getCumulativeActiveOrbitalsPerIrrep(const int irrep) const;
    int getOriginalOccupiedOrbitalsPerIrrep(const int irrep) const;
    int getOriginalActiveOrbitalsPerIrrep(const int irrep) const;
    int getOriginalVirtualOrbitalsPerIrrep(const int irrep) const;
    int* getIrrepPerDMRG();
    int getIrrepPerOrbital(const int index) const;
    int getOccupiedOrbitalsPerIrrepSum() const;
    int getActiveOrbitalsPerIrrepSum() const;
    int getVirtualOrbitalsPerIrrepSum() const;
    void prnint() const;
};
