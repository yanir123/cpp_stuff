#pragma once

#ifdef CHEMPS2_MPI_COMPILATION

#include <assert.h>
#include <mpi.h>

#include "Tensor.h"

#define MPI_CHEMPS2_MASTER 0

#define MPI_CHEMPS2_4D1AB 1
#define MPI_CHEMPS2_4D2AB 2
#define MPI_CHEMPS2_4I1AB 3
#define MPI_CHEMPS2_4I2AB 4
#define MPI_CHEMPS2_4F1AB 5
#define MPI_CHEMPS2_4F2AB 6
#define MPI_CHEMPS2_4G1AB 7
#define MPI_CHEMPS2_4G2AB 8

#define MPI_CHEMPS2_4D3ABCD 9
#define MPI_CHEMPS2_4D4ABCD 10
#define MPI_CHEMPS2_4I3ABCD 11
#define MPI_CHEMPS2_4I4ABCD 12
#define MPI_CHEMPS2_4F3ABCD 13
#define MPI_CHEMPS2_4F4ABCD 14
#define MPI_CHEMPS2_4G3ABCD 15
#define MPI_CHEMPS2_4G4ABCD 16

#define MPI_CHEMPS2_4E1 17
#define MPI_CHEMPS2_4E2 18
#define MPI_CHEMPS2_4H1 19
#define MPI_CHEMPS2_4H2 20

#define MPI_CHEMPS2_4E3A 21
#define MPI_CHEMPS2_4E3B 22
#define MPI_CHEMPS2_4E4A 23
#define MPI_CHEMPS2_4E4B 24
#define MPI_CHEMPS2_4H3A 25
#define MPI_CHEMPS2_4H3B 26
#define MPI_CHEMPS2_4H4A 27
#define MPI_CHEMPS2_4H4B 28

#define MPI_CHEMPS2_5A1 29
#define MPI_CHEMPS2_5A2 30
#define MPI_CHEMPS2_5A3 31
#define MPI_CHEMPS2_5A4 32

#define MPI_CHEMPS2_5B1 33
#define MPI_CHEMPS2_5B2 34
#define MPI_CHEMPS2_5B3 35
#define MPI_CHEMPS2_5B4 36

#define MPI_CHEMPS2_5C1 37
#define MPI_CHEMPS2_5C2 38
#define MPI_CHEMPS2_5C3 39
#define MPI_CHEMPS2_5C4 40

#define MPI_CHEMPS2_5D1 41
#define MPI_CHEMPS2_5D2 42
#define MPI_CHEMPS2_5D3 43
#define MPI_CHEMPS2_5D4 44

#define MPI_CHEMPS2_5E1 45
#define MPI_CHEMPS2_5E2 46
#define MPI_CHEMPS2_5E3 47
#define MPI_CHEMPS2_5E4 48

#define MPI_CHEMPS2_5F1 49
#define MPI_CHEMPS2_5F2 50
#define MPI_CHEMPS2_5F3 51
#define MPI_CHEMPS2_5F4 52

#define MPI_CHEMPS2_OFFSET 53

#endif

class MPI {
   public:
    MPI() {}
    virtual ~MPI() {}
    static int mpi_size() {
#ifdef CHEMPS2_MPI_COMPILATION
        int size;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        return size;
#else
        return 1;
#endif
    }
    static int mpi_rank() {
#ifdef CHEMPS2_MPI_COMPILATION
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        return rank;
#else
        return 0;
#endif
    }

#ifdef CHEMPS2_MPI_COMPILATION
    static void mpi_init() {
        int zero = 0;
        MPI_Init(&zero, NULL);
    }
#endif

#ifdef CHEMPS2_MPI_COMPILATION
    static void mpi_finalize() {
        MPI_Finalize();
    }
#endif

#ifdef CHEMPS2_MPI_COMPILATION
    static int owner_x() { return MPI_CHEMPS2_MASTER; }
#endif

#ifdef CHEMPS2_MPI_COMPILATION
    static int owner_absigma(const int index1, const int index2) {
        assert(index1 <= index2);
        return (1 + index1 + (index2 * (index2 + 1)) / 2) % mpi_size();
    }
#endif

#ifdef CHEMPS2_MPI_COMPILATION
    static int owner_cdf(const int L, const int index1, const int index2) { 
        assert(index1 <= index2);
        return (1 + (L * (L + 1)) / 2 + index1 + (index2 * (index2 + 1)) / 2) % mpi_size();
    }
#endif

#ifdef CHEMPS2_MPI_COMPILATION
    static int owner_q(const int L, const int index) { 
        return (1 + L * (L + 1) + index) % mpi_size();
    }
#endif

#ifdef CHEMPS2_MPI_COMPILATION
    static int owner_3rdm_diagram(const int L, const int index1, const int index2, const int index3) {  
        assert(index1 <= index2);
        assert(index2 <= index3);
        return (1 + L * (L + 1) + index1 + (index2 * (index2 + 1)) / 2 + (index3 * (index3 + 1) * (index3 + 2)) / 6) % mpi_size();
    }
#endif

#ifdef CHEMPS2_MPI_COMPILATION
    static int owner_1cd2d3eh() { return MPI_CHEMPS2_MASTER; }
#endif

#ifdef CHEMPS2_MPI_COMPILATION
    static int owner_specific_diagram(const int L, const int macro) {
        return (macro + L * (L + 2)) % mpi_size();
    }
#endif

#ifdef CHEMPS2_MPI_COMPILATION
    static int owner_specific_excitation(const int L, const int excitation) {
        return (MPI_CHEMPS2_OFFSET + L * (L + 2) + excitation) % mpi_size();
    }
#endif

#ifdef CHEMPS2_MPI_COMPILATION
    static void broadcast_tensor(Tensor* object, int ROOT) {
        int arraysize = object->gKappa2index(object->gNKappa());
        MPI_Bcast(object->gStorage(), arraysize, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    }
#endif

#ifdef CHEMPS2_MPI_COMPILATION
    static void broadcast_array_double(double* array, int length, int ROOT) {
        MPI_Bcast(array, length, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
    }
#endif

#ifdef CHEMPS2_MPI_COMPILATION
    static void broadcast_array_int(int* array, int length, int ROOT) {
        MPI_Bcast(array, length, MPI_INT, ROOT, MPI_COMM_WORLD);
    }
#endif

#ifdef CHEMPS2_MPI_COMPILATION
    static bool all_booleans_equal(const bool mybool) {
        int my_value = (mybool) ? 1 : 0;
        int tot_value;
        MPI_Allreduce(&my_value, &tot_value, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        return (my_value * MPIchemps2::mpi_size() == tot_value);
    }
#endif

#ifdef CHEMPS2_MPI_COMPILATION
    static void sendreceive_tensor(Tensor* object, int SENDER, int RECEIVER, int tag) {
        if (SENDER != RECEIVER) {
            const int MPIRANK = mpi_rank();
            if (SENDER == MPIRANK) {
                int arraysize = object->gKappa2index(object->gNKappa());
                MPI_Send(object->gStorage(), arraysize, MPI_DOUBLE, RECEIVER, tag, MPI_COMM_WORLD);
            }
            if (RECEIVER == MPIRANK) {
                int arraysize = object->gKappa2index(object->gNKappa());
                MPI_Recv(object->gStorage(), arraysize, MPI_DOUBLE, SENDER, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }
#endif

#ifdef CHEMPS2_MPI_COMPILATION
    static void reduce_array_double(double* vec_in, double* vec_out, int size, int ROOT) {
        MPI_Reduce(vec_in, vec_out, size, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
    }
#endif

#ifdef CHEMPS2_MPI_COMPILATION
    static void allreduce_array_double(double* vec_in, double* vec_out, int size) {
        MPI_Allreduce(vec_in, vec_out, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
#endif
};