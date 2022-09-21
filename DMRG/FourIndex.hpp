#pragma once

#include <hdf5.h>

#include <cassert>
#include <climits>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

#include "Irreps.hpp"
#include "Lapack.hpp"
#include "Mpi.hpp"

class FourIndex {
   private:
    Irreps symmInfo;
    int* irrepsNumber;
    long long***** storage;
    long long calcNumberOfUniqueElements(const bool allocateStorage);
    long long arrayLength;
    double* elements;
    long long getPointer(const int irrep_i, const int irrep_j, const int irrep_k, const int irrep_l, const int i, const int j, const int k, const int l) const;
    long long getPtrIrrepOrderOK(const int irrep_i, const int irrep_j, const int irrep_k, const int irrep_l, const int i, const int j, const int k, const int l) const;
    long long getPtrAllOK1(const int Icent, const int irrep_i, const int irrep_k, const int i, const int j, const int k, const int l) const;
    long long getPtrAllOK2(const int Icent, const int irrep_i, const int irrep_k, const int i, const int j, const int k, const int l) const;
    long long getPtrAllOK5(const int Icent, const int irrep_i, const int irrep_k, const int i, const int j, const int k, const int l) const;

   public:
    FourIndex(const int group, const int* irrepsNumber);
    virtual ~FourIndex();
    void clear();
    void set(const int irrep_i, const int irrep_j, const int irrep_k, const int irrep_l, const int i, const int j, const int k, const int l, const double val);
    void add(const int irrep_i, const int irrep_j, const int irrep_k, const int irrep_l, const int i, const int j, const int k, const int l, const double val);
    double get(const int irrep_i, const int irrep_j, const int irrep_k, const int irrep_l, const int i, const int j, const int k, const int l) const;
    int get_irrep_size(const int irrep) const;
    void save(const std::string name) const;
    void read(const std::string name);
#ifdef CHEMPS2_MPI_COMPILATION
    void broadcast(const int ROOT);
#endif
};
