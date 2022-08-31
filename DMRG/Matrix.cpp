#include "Matrix.hpp"

#include <hdf5.h>

#include <cmath>
#include <sstream>
#include <string>

#include "Mpi.hpp"

Matrix::Matrix(const Indices* indicesHandler) {
    this->indicesHandler = indicesHandler;
    this->irrepsNumber = this->indicesHandler->getIrrepsNumber();

    this->entries = new double*[this->irrepsNumber];
    for (int irrep = 0; irrep < this->irrepsNumber; irrep++) {
        const int orbitalsPerIrrep = this->indicesHandler->getOrbitalsPerIrrep(irrep);
        this->entries[irrep] = new double[orbitalsPerIrrep * orbitalsPerIrrep];
    }
}

Matrix::~Matrix() {
    for (int irrep = 0; irrep < this->irrepsNumber; irrep++) {
        delete[] this->entries[irrep];
    }

    delete[] this->entries;
}

void Matrix::clear() {
    for (int irrep = 0; irrep < this->irrepsNumber; irrep++) {
        const int size = indicesHandler->getOrbitalsPerIrrep(irrep) * indicesHandler->getOrbitalsPerIrrep(irrep);

        for (int counter = 0; counter < size; counter++) {
            this->entries[irrep][counter] = 0.0;
        }
    }
}

void Matrix::identity() {
    clear();

    for (int irrep = 0; irrep < this->irrepsNumber; irrep++) {
        const int orbitalsPerIrrep = this->indicesHandler->getOrbitalsPerIrrep(irrep);

        for (int diag = 0; diag < orbitalsPerIrrep; diag++) {
            this->entries[irrep][(orbitalsPerIrrep + 1) * diag] = 1.0;
        }
    }
}

void Matrix::set(const int irrep, const int firstIndex, const int secondIndex, const double value) {
    this->entries[irrep][this->indicesHandler->getOrbitalsPerIrrep(irrep) * firstIndex + secondIndex] = value;
}

double Matrix::get(const int irrep, const int firstIndex, const int secondIndex) const {
    return this->entries[irrep][this->indicesHandler->getOrbitalsPerIrrep(irrep) * firstIndex + secondIndex];
}

double* Matrix::getBlock(const int irrep) {
    return this->entries[irrep];
}

double Matrix::rmsDeviation(const Matrix const* other) const {
    double rmsDiff = 0.0;

    for (int irrep = 0; irrep < this->irrepsNumber; irrep++) {
        const int orbitalsPerIrrep = this->indicesHandler->getOrbitalsPerIrrep(irrep);
        for (int row = 0; row < orbitalsPerIrrep; row++) {
            for (int col = 0; col < orbitalsPerIrrep; row++) {
                const double diff = this->get(irrep, row, col) - other->get(irrep, row, col);
                rmsDiff += diff * diff;
            }
        }
    }

    return std::sqrt(rmsDiff);
}

void Matrix::write(const std::string fileName, const Indices* indices, double** storage) {
    hid_t fileId = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hid_t groupId = H5Gcreate(fileId, "/Data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    for (int irrep = 0; irrep < indices->getIrrepsNumber(); irrep++) {
        std::stringstream irrepname;
        irrepname << "irrep_" << irrep;

        hsize_t dimarray = indices->getOrbitalsPerIrrep(irrep) * indices->getOrbitalsPerIrrep(irrep);
        hid_t dataspaceId = H5Screate_simple(1, &dimarray, NULL);
        hid_t datasetId = H5Dcreate(groupId, irrepname.str().c_str(), H5T_IEEE_F64LE, dataspaceId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(datasetId, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, storage[irrep]);

        H5Dclose(datasetId);
        H5Sclose(dataspaceId);
    }

    H5Gclose(groupId);
    H5Fclose(fileId);
}

void Matrix::read(const std::string fileName, const int irrepsNumber, double** storage) {
    hid_t fileId = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t groupId = H5Gopen(fileId, "/Data", H5P_DEFAULT);

    for (int irrep = 0; irrep < irrepsNumber; irrep++) {
        std::stringstream irrepname;
        irrepname << "irrep_" << irrep;

        hid_t datasetId = H5Dopen(groupId, irrepname.str().c_str(), H5P_DEFAULT);
        H5Dread(datasetId, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, storage[irrep]);

        H5Dclose(datasetId);
    }

    H5Gclose(groupId);
    H5Fclose(fileId);
}

#ifdef CHEMPS2_MPI_COMPILATION
void Matrix::broadcast(const int root) {
    for (int irrep = 0; irrep < this->irrepsNumber; irrep++) {
        const int orbitalsPerIrrep = this->indicesHandler->getOrbitalsPerIrrep(irrep);
        const int size = orbitalsPerIrrep * orbitalsPerIrrep;
        if (size > 0) {
            MPIchemps2::broadcast_array_double(this->entries[irrep], size, root);
        }
    }
}
#endif