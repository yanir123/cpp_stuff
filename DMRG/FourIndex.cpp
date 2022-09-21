#include "FourIndex.hpp"

FourIndex::FourIndex(const int group, const int* irrepsNumber) {
    this->symmInfo.setGroup(group);

    this->irrepsNumber = new int[this->symmInfo.numberOfIrreps()];
    for (int irrep = 0; irrep < this->symmInfo.numberOfIrreps(); irrep++) {
        this->irrepsNumber[irrep] = irrepsNumber[irrep];
    }

    this->arrayLength = this->calcNumberOfUniqueElements(true);
    this->elements = new double[this->arrayLength];

    clear();
}

FourIndex::~FourIndex() {
    this->arrayLength = this->calcNumberOfUniqueElements(false);
    delete[] this->elements;
    delete[] this->irrepsNumber;
}

long long FourIndex::calcNumberOfUniqueElements(const bool allocateStorage) {
    long long theTotalSize = 0;

    if (allocateStorage) {
        this->storage = new long long****[this->symmInfo.numberOfIrreps()];
    }
    for (int Icenter = 0; Icenter < this->symmInfo.numberOfIrreps(); Icenter++) {
        if (allocateStorage) {
            this->storage[Icenter] = new long long***[this->symmInfo.numberOfIrreps()];
        }
        for (int I_i = 0; I_i < this->symmInfo.numberOfIrreps(); I_i++) {
            int I_j = Irreps::directProduct(Icenter, I_i);
            if ((this->irrepsNumber[I_i] > 0) && (this->irrepsNumber[I_j] > 0)) {
                if (allocateStorage) {
                    this->storage[Icenter][I_i] = new long long**[this->symmInfo.numberOfIrreps()];
                }
                for (int I_k = I_i; I_k < this->symmInfo.numberOfIrreps(); I_k++) {
                    int I_l = Irreps::directProduct(Icenter, I_k);
                    if ((this->irrepsNumber[I_k] > 0) && (this->irrepsNumber[I_l] > 0)) {
                        if ((I_i <= I_j) && (I_j <= I_l)) {
                            if (Icenter == 0) {  // I_i = I_j and I_k = I_l
                                if (I_i == I_k) {
                                    if (allocateStorage) {
                                        this->storage[Icenter][I_i][I_k] = new long long*[this->irrepsNumber[I_i] * (this->irrepsNumber[I_i] + 1) / 2];
                                    }
                                    for (int i = 0; i < this->irrepsNumber[I_i]; i++) {
                                        for (int k = i; k < this->irrepsNumber[I_k]; k++) {
                                            if (allocateStorage) {
                                                this->storage[Icenter][I_i][I_k][i + k * (k + 1) / 2] = new long long[this->irrepsNumber[I_j] - i];
                                                for (int j = i; j < this->irrepsNumber[I_j]; j++) {
                                                    this->storage[Icenter][I_i][I_k][i + k * (k + 1) / 2][j - i] = theTotalSize;
                                                    theTotalSize += this->irrepsNumber[I_l] - ((i == j) ? k : j);
                                                }
                                            } else {
                                                delete[] this->storage[Icenter][I_i][I_k][i + k * (k + 1) / 2];
                                            }
                                        }
                                    }
                                    if (!allocateStorage) {
                                        delete[] this->storage[Icenter][I_i][I_k];
                                    }
                                } else {  // I_i < I_k
                                    if (allocateStorage) {
                                        this->storage[Icenter][I_i][I_k] = new long long*[this->irrepsNumber[I_i] * this->irrepsNumber[I_k]];
                                    }
                                    for (int i = 0; i < this->irrepsNumber[I_i]; i++) {
                                        for (int k = 0; k < this->irrepsNumber[I_k]; k++) {
                                            if (allocateStorage) {
                                                this->storage[Icenter][I_i][I_k][i + k * this->irrepsNumber[I_i]] = new long long[this->irrepsNumber[I_j] - i];
                                                for (int j = i; j < this->irrepsNumber[I_j]; j++) {
                                                    this->storage[Icenter][I_i][I_k][i + k * this->irrepsNumber[I_i]][j - i] = theTotalSize;
                                                    theTotalSize += this->irrepsNumber[I_l] - ((i == j) ? k : 0);
                                                }
                                            } else {
                                                delete[] this->storage[Icenter][I_i][I_k][i + k * this->irrepsNumber[I_i]];
                                            }
                                        }
                                    }
                                    if (!allocateStorage) {
                                        delete[] this->storage[Icenter][I_i][I_k];
                                    }
                                }
                            } else {  //Icenter !=0 ; I_i < I_j and I_k != I_l
                                if (I_i == I_k) {
                                    if (allocateStorage) {
                                        this->storage[Icenter][I_i][I_k] = new long long*[this->irrepsNumber[I_i] * (this->irrepsNumber[I_i] + 1) / 2];
                                    }
                                    for (int i = 0; i < this->irrepsNumber[I_i]; i++) {
                                        for (int k = i; k < this->irrepsNumber[I_k]; k++) {
                                            if (allocateStorage) {
                                                this->storage[Icenter][I_i][I_k][i + k * (k + 1) / 2] = new long long[this->irrepsNumber[I_j]];
                                                for (int j = 0; j < this->irrepsNumber[I_j]; j++) {
                                                    this->storage[Icenter][I_i][I_k][i + k * (k + 1) / 2][j] = theTotalSize;
                                                    theTotalSize += this->irrepsNumber[I_l] - j;
                                                }
                                            } else {
                                                delete[] this->storage[Icenter][I_i][I_k][i + k * (k + 1) / 2];
                                            }
                                        }
                                    }
                                    if (!allocateStorage) {
                                        delete[] this->storage[Icenter][I_i][I_k];
                                    }
                                } else {  // I_i < I_k
                                    if (allocateStorage) {
                                        this->storage[Icenter][I_i][I_k] = new long long*[this->irrepsNumber[I_i] * this->irrepsNumber[I_k]];
                                    }
                                    for (int i = 0; i < this->irrepsNumber[I_i]; i++) {
                                        for (int k = 0; k < this->irrepsNumber[I_k]; k++) {
                                            if (allocateStorage) {
                                                this->storage[Icenter][I_i][I_k][i + k * this->irrepsNumber[I_i]] = new long long[this->irrepsNumber[I_j]];
                                                for (int j = 0; j < this->irrepsNumber[I_j]; j++) {
                                                    this->storage[Icenter][I_i][I_k][i + k * this->irrepsNumber[I_i]][j] = theTotalSize;
                                                    theTotalSize += this->irrepsNumber[I_l];
                                                }
                                            } else {
                                                delete[] this->storage[Icenter][I_i][I_k][i + k * this->irrepsNumber[I_i]];
                                            }
                                        }
                                    }
                                    if (!allocateStorage) {
                                        delete[] this->storage[Icenter][I_i][I_k];
                                    }
                                }
                            }
                        }
                    }
                }
                if (!allocateStorage) {
                    delete[] this->storage[Icenter][I_i];
                }
            }
        }
        if (!allocateStorage) {
            delete[] this->storage[Icenter];
        }
    }
    if (!allocateStorage) {
        delete[] this->storage;
    }

    return theTotalSize;
}

long long FourIndex::getPointer(const int irrep_i, const int irrep_j, const int irrep_k, const int irrep_l, const int i, const int j, const int k, const int l) const {
    assert(Irreps::directProduct(irrep_i, irrep_j) == Irreps::directProduct(irrep_k, irrep_l));

    if ((irrep_i <= irrep_j) && (irrep_i <= irrep_k) && (irrep_j <= irrep_l)) {
        return getPtrIrrepOrderOK(irrep_i, irrep_j, irrep_k, irrep_l, i, j, k, l);
    }

    if ((irrep_j <= irrep_i) && (irrep_j <= irrep_l) && (irrep_i <= irrep_k)) {
        return getPtrIrrepOrderOK(irrep_j, irrep_i, irrep_l, irrep_k, j, i, l, k);
    }

    if ((irrep_k <= irrep_j) && (irrep_k <= irrep_i) && (irrep_j <= irrep_l)) {
        return getPtrIrrepOrderOK(irrep_k, irrep_j, irrep_i, irrep_l, k, j, i, l);
    }

    if ((irrep_j <= irrep_k) && (irrep_j <= irrep_l) && (irrep_k <= irrep_i)) {
        return getPtrIrrepOrderOK(irrep_j, irrep_k, irrep_l, irrep_i, j, k, l, i);
    }

    if ((irrep_i <= irrep_l) && (irrep_i <= irrep_k) && (irrep_l <= irrep_j)) {
        return getPtrIrrepOrderOK(irrep_i, irrep_l, irrep_k, irrep_j, i, l, k, j);
    }

    if ((irrep_l <= irrep_i) && (irrep_l <= irrep_j) && (irrep_i <= irrep_k)) {
        return getPtrIrrepOrderOK(irrep_l, irrep_i, irrep_j, irrep_k, l, i, j, k);
    }

    if ((irrep_k <= irrep_l) && (irrep_k <= irrep_i) && (irrep_l <= irrep_j)) {
        return getPtrIrrepOrderOK(irrep_k, irrep_l, irrep_i, irrep_j, k, l, i, j);
    }

    if ((irrep_l <= irrep_k) && (irrep_l <= irrep_j) && (irrep_k <= irrep_i)) {
        return getPtrIrrepOrderOK(irrep_l, irrep_k, irrep_j, irrep_i, l, k, j, i);
    }

    return -1;
}

long long FourIndex::getPtrIrrepOrderOK(const int irrep_i, const int irrep_j, const int irrep_k, const int irrep_l, const int i, const int j, const int k, const int l) const {
    int Icenter = Irreps::directProduct(irrep_i, irrep_j);

    if (Icenter > 0) {
        if (irrep_i == irrep_k) {
            if (k >= i) {
                if (l >= j)
                    return getPtrAllOK5(Icenter, irrep_i, irrep_k, i, j, k, l);
                else
                    return getPtrAllOK5(Icenter, irrep_i, irrep_k, i, l, k, j);
            } else {
                if (l >= j)
                    return getPtrAllOK5(Icenter, irrep_k, irrep_i, k, j, i, l);
                else
                    return getPtrAllOK5(Icenter, irrep_k, irrep_i, k, l, i, j);
            }

        } else {
            return this->storage[Icenter][irrep_i][irrep_k][i + this->irrepsNumber[irrep_i] * k][j] + l;
        }

    } else {
        if (irrep_i == irrep_k) {
            if ((i < j) && (i <= k) && (j <= l)) return getPtrAllOK2(Icenter, irrep_i, irrep_k, i, j, k, l);
            if ((i == j) && (i <= k) && (j <= l)) {
                if (l >= k)
                    return getPtrAllOK1(Icenter, irrep_i, irrep_k, i, j, k, l);
                else
                    return getPtrAllOK1(Icenter, irrep_j, irrep_l, j, i, l, k);
            }
            if ((j < i) && (j <= l) && (i <= k)) return getPtrAllOK2(Icenter, irrep_j, irrep_l, j, i, l, k);

            if ((k < l) && (k <= i) && (l <= j)) return getPtrAllOK2(Icenter, irrep_k, irrep_i, k, l, i, j);
            if ((k == l) && (k <= i) && (l <= j)) {
                if (j >= i)
                    return getPtrAllOK1(Icenter, irrep_k, irrep_i, k, l, i, j);
                else
                    return getPtrAllOK1(Icenter, irrep_l, irrep_j, l, k, j, i);
            }
            if ((l < k) && (l <= j) && (k <= i)) return getPtrAllOK2(Icenter, irrep_l, irrep_j, l, k, j, i);

            if ((k < j) && (k <= i) && (j <= l)) return getPtrAllOK2(Icenter, irrep_k, irrep_i, k, j, i, l);
            if ((k == j) && (k <= i) && (j <= l)) {
                if (l >= i)
                    return getPtrAllOK1(Icenter, irrep_k, irrep_i, k, j, i, l);
                else
                    return getPtrAllOK1(Icenter, irrep_j, irrep_l, j, k, l, i);
            }
            if ((j < k) && (j <= l) && (k <= i)) return getPtrAllOK2(Icenter, irrep_j, irrep_l, j, k, l, i);

            if ((i < l) && (i <= k) && (l <= j)) return getPtrAllOK2(Icenter, irrep_i, irrep_k, i, l, k, j);
            if ((i == l) && (i <= k) && (l <= j)) {
                if (j >= k)
                    return getPtrAllOK1(Icenter, irrep_i, irrep_k, i, l, k, j);
                else
                    return getPtrAllOK1(Icenter, irrep_l, irrep_j, l, i, j, k);
            }
            if ((l < i) && (l <= j) && (i <= k)) return getPtrAllOK2(Icenter, irrep_l, irrep_j, l, i, j, k);

        } else {
            if (j == i) {
                if (l >= k) {
                    return this->storage[Icenter][irrep_i][irrep_k][i + this->irrepsNumber[irrep_i] * k][j - i] + l - k;
                } else {
                    return this->storage[Icenter][irrep_j][irrep_l][j + this->irrepsNumber[irrep_j] * l][i - j] + k - l;
                }
            } else {
                if (j > i) {
                    return this->storage[Icenter][irrep_i][irrep_k][i + this->irrepsNumber[irrep_i] * k][j - i] + l;
                } else {
                    return this->storage[Icenter][irrep_j][irrep_l][j + this->irrepsNumber[irrep_j] * l][i - j] + k;
                }
            }
        }
    }

    return -1;
}

long long FourIndex::getPtrAllOK1(const int Icent, const int irrep_i, const int irrep_k, const int i, const int j, const int k, const int l) const {
    return this->storage[Icent][irrep_i][irrep_k][i + k * (k + 1) / 2][j - i] + l - k;
}

long long FourIndex::getPtrAllOK2(const int Icent, const int irrep_i, const int irrep_k, const int i, const int j, const int k, const int l) const {
    return this->storage[Icent][irrep_i][irrep_k][i + k * (k + 1) / 2][j - i] + l - j;
}

long long FourIndex::getPtrAllOK5(const int Icent, const int irrep_i, const int irrep_k, const int i, const int j, const int k, const int l) const {
    return this->storage[Icent][irrep_i][irrep_k][i + k * (k + 1) / 2][j] + l - j;
}

void FourIndex::clear() {
    for (long long i = 0; i < this->arrayLength; i++) {
        this->elements[i] = 0.0;
    }
}

void FourIndex::set(const int irrep_i, const int irrep_j, const int irrep_k, const int irrep_l, const int i, const int j, const int k, const int l, const double val) {
    this->elements[this->getPointer(irrep_i, irrep_j, irrep_k, irrep_l, i, j, k, l)] = val;
}

void FourIndex::add(const int irrep_i, const int irrep_j, const int irrep_k, const int irrep_l, const int i, const int j, const int k, const int l, const double val) {
    this->elements[this->getPointer(irrep_i, irrep_j, irrep_k, irrep_l, i, j, k, l)] += val;
}

double FourIndex::get(const int irrep_i, const int irrep_j, const int irrep_k, const int irrep_l, const int i, const int j, const int k, const int l) const {
    this->elements[this->getPointer(irrep_i, irrep_j, irrep_k, irrep_l, i, j, k, l)];
}

int FourIndex::get_irrep_size(const int irrep) const {
    return this->irrepsNumber[irrep];
}

void FourIndex::save(const std::string name) const {
    hid_t file_id = H5Fcreate(name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    hid_t group_id = H5Gcreate(file_id, "/MetaData", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    hsize_t dimarray = this->symmInfo.numberOfIrreps();
    hid_t dataspace_id = H5Screate_simple(1, &dimarray, NULL);
    hid_t dataset_id = H5Dcreate(group_id, "IrrepSizes", H5T_STD_I32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, this->irrepsNumber);

    hid_t attribute_space_id1 = H5Screate(H5S_SCALAR);
    hid_t attribute_id1 = H5Acreate(dataset_id, "nGroup", H5T_STD_I32LE, attribute_space_id1, H5P_DEFAULT, H5P_DEFAULT);
    int nGroup = this->symmInfo.getGroupNumber();
    H5Awrite(attribute_id1, H5T_NATIVE_INT, &nGroup);

    hid_t attribute_space_id2 = H5Screate(H5S_SCALAR);
    hid_t attribute_id2 = H5Acreate(dataset_id, "nIrreps", H5T_STD_I32LE, attribute_space_id2, H5P_DEFAULT, H5P_DEFAULT);
    int nIrreps = this->symmInfo.numberOfIrreps();
    H5Awrite(attribute_id2, H5T_NATIVE_INT, &nIrreps);

    hid_t attribute_space_id3 = H5Screate(H5S_SCALAR);
    hid_t attribute_id3 = H5Acreate(dataset_id, "theTotalSize", H5T_STD_I64LE, attribute_space_id3, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attribute_id3, H5T_NATIVE_LLONG, &arrayLength);

    H5Aclose(attribute_id1);
    H5Aclose(attribute_id2);
    H5Aclose(attribute_id3);
    H5Sclose(attribute_space_id1);
    H5Sclose(attribute_space_id2);
    H5Sclose(attribute_space_id3);

    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    H5Gclose(group_id);

    hid_t group_id7 = H5Gcreate(file_id, "/FourIndexObject", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    hsize_t dimarray7 = arrayLength;
    hid_t dataspace_id7 = H5Screate_simple(1, &dimarray7, NULL);
    hid_t dataset_id7 = H5Dcreate(group_id7, "Matrix elements", H5T_IEEE_F64LE, dataspace_id7, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id7, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, this->elements);

    H5Dclose(dataset_id7);
    H5Sclose(dataspace_id7);

    H5Gclose(group_id7);

    H5Fclose(file_id);
}

void FourIndex::read(const std::string name) {
    hid_t file_id = H5Fopen(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    hid_t group_id = H5Gopen(file_id, "/MetaData", H5P_DEFAULT);

    hid_t dataset_id = H5Dopen(group_id, "IrrepSizes", H5P_DEFAULT);

    hid_t attribute_id1 = H5Aopen_by_name(group_id, "IrrepSizes", "nGroup", H5P_DEFAULT, H5P_DEFAULT);
    int nGroup;
    H5Aread(attribute_id1, H5T_NATIVE_INT, &nGroup);
    assert(nGroup == this->symmInfo.getGroupNumber());

    hid_t attribute_id2 = H5Aopen_by_name(group_id, "IrrepSizes", "nIrreps", H5P_DEFAULT, H5P_DEFAULT);
    int nIrreps;
    H5Aread(attribute_id2, H5T_NATIVE_INT, &nIrreps);
    assert(nIrreps == this->symmInfo.numberOfIrreps());

    hid_t attribute_id3 = H5Aopen_by_name(group_id, "IrrepSizes", "theTotalSize", H5P_DEFAULT, H5P_DEFAULT);
    long long theTotalSize;
    H5Aread(attribute_id3, H5T_NATIVE_LLONG, &theTotalSize);
    assert(theTotalSize == arrayLength);

    H5Aclose(attribute_id1);
    H5Aclose(attribute_id2);
    H5Aclose(attribute_id3);

    int* IsizesAgain = new int[this->symmInfo.numberOfIrreps()];
    H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, IsizesAgain);
    for (int cnt = 0; cnt < this->symmInfo.numberOfIrreps(); cnt++) {
        assert(IsizesAgain[cnt] == this->irrepsNumber[cnt]);
    }
    delete[] IsizesAgain;
    H5Dclose(dataset_id);

    H5Gclose(group_id);

    std::cout << "FourIndex::read : loading " << arrayLength << " doubles." << std::endl;

    hid_t group_id7 = H5Gopen(file_id, "/FourIndexObject", H5P_DEFAULT);

    hid_t dataset_id7 = H5Dopen(group_id7, "Matrix elements", H5P_DEFAULT);
    H5Dread(dataset_id7, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, this->elements);
    H5Dclose(dataset_id7);

    H5Gclose(group_id7);

    H5Fclose(file_id);

    std::cout << "FourIndex::read : everything loaded!" << std::endl;
}

#ifdef CHEMPS2_MPI_COMPILATION
void FourIndex::broadcast(const int ROOT) {
    assert(this->arrayLength <= INT_MAX);
    const int size = (int)(this->arrayLength);
    if (size > 0) {
        Mpi::broadcast_array_double(this->elements, size, ROOT);
    }
}
#endif
