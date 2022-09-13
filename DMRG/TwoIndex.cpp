#include "TwoIndex.hpp"

#include <hdf5.h>

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

TwoIndex::TwoIndex(const int group, const int* irrepsNumber) {
    this->symmInfo.setGroup(group);

    this->irrepsNumber = new int[this->symmInfo.numberOfIrreps()];
    this->storage = new double*[this->symmInfo.numberOfIrreps()];

    for (int i = 0; i < this->symmInfo.numberOfIrreps(); i++) {
        this->irrepsNumber[i] = irrepsNumber[i];
        if (this->irrepsNumber[i] > 0) {
            this->storage[i] = new double[this->irrepsNumber[i] * (this->irrepsNumber[i] + 1) / 2];
        }
    }

    this->clear();
}

TwoIndex::~TwoIndex() {
    for (int i = 0; i < this->symmInfo.numberOfIrreps(); i++) {
        if (this->irrepsNumber[i] > 0) {
            delete[] this->storage[i];
        }
    }

    delete[] this->irrepsNumber;
    delete[] this->storage;
}

void TwoIndex::clear() {
    for (int i = 0; i < this->symmInfo.numberOfIrreps(); i++) {
        const int loopSize = this->irrepsNumber[i] * (this->irrepsNumber[i] + 1) / 2;

        for (int j = 0; j < loopSize; j++) {
            this->storage[i][j] = 0.0;
        }
    }
}

void TwoIndex::set(const int irrep, const int i, const int j, const int value) {
    if (i > j) {
        this->storage[irrep][j + i * (i + 1) / 2] = value;
    } else {
        this->storage[irrep][i + j * (j + 1) / 2] = value;
    }
}

double TwoIndex::get(const int irrep, const int i, const int j) const {
    if (i > j) {
        return this->storage[irrep][j + i * (i + 1) / 2];
    }

    return this->storage[irrep][i + j * (j + 1) / 2];
}

void TwoIndex::save(std::string fileName) const {
    hid_t file_id = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

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

    H5Aclose(attribute_id1);
    H5Aclose(attribute_id2);
    H5Sclose(attribute_space_id1);
    H5Sclose(attribute_space_id2);

    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    H5Gclose(group_id);

    for (int cnt = 0; cnt < this->symmInfo.numberOfIrreps(); cnt++) {
        if (this->irrepsNumber[cnt] > 0) {
            std::stringstream sstream;
            sstream << "/TwoIndex" << cnt;
            hid_t group_id3 = H5Gcreate(file_id, sstream.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

            hsize_t dimarray3 = this->irrepsNumber[cnt] * (this->irrepsNumber[cnt] + 1) / 2;
            hid_t dataspace_id3 = H5Screate_simple(1, &dimarray3, NULL);
            hid_t dataset_id3 = H5Dcreate(group_id3, "Matrix elements", H5T_IEEE_F64LE, dataspace_id3, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Dwrite(dataset_id3, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, this->storage[cnt]);

            H5Dclose(dataset_id3);
            H5Sclose(dataspace_id3);

            H5Gclose(group_id3);
        }
    }

    H5Fclose(file_id);
}

void TwoIndex::red(std::string fileName) {
    hid_t file_id = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

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

    H5Aclose(attribute_id1);
    H5Aclose(attribute_id2);

    int* IsizesAgain = new int[this->symmInfo.numberOfIrreps()];
    H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, IsizesAgain);
    for (int cnt = 0; cnt < this->symmInfo.numberOfIrreps(); cnt++) {
        assert(IsizesAgain[cnt] == this->irrepsNumber[cnt]);
    }
    delete[] IsizesAgain;
    H5Dclose(dataset_id);

    H5Gclose(group_id);

    for (int cnt = 0; cnt < this->symmInfo.numberOfIrreps(); cnt++) {
        if (this->irrepsNumber[cnt] > 0) {
            std::stringstream sstream;
            sstream << "/TwoIndex" << cnt;
            hid_t group_id3 = H5Gopen(file_id, sstream.str().c_str(), H5P_DEFAULT);

            hid_t dataset_id3 = H5Dopen(group_id3, "Matrix elements", H5P_DEFAULT);
            H5Dread(dataset_id3, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, this->storage[cnt]);
            H5Dclose(dataset_id3);

            H5Gclose(group_id3);
        }
    }

    H5Fclose(file_id);
}
