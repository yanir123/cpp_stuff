#pragma once

#include "Irreps.hpp"

#include <string>

class TwoIndex
{
private:
    Irreps symmInfo;
    int* irrepsNumber;
    double** storage;
public:
    TwoIndex(const int group, const int* irrepsNumber);
    virtual ~TwoIndex();
    void clear();
    void set(const int irrep, const int i, const int j, const int value);
    double get(const int irrep, const int i, const int j) const;
    void save(std::string fileName) const;
    void red(std::string fileName);
};