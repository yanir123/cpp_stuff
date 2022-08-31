#include "Indices.hpp"

class Matrix {
   protected:
    const Indices* indicesHandler;
    double** entries;
    int irrepsNumber;

   public:
    Matrix(const Indices* indicesHandler);
    virtual ~Matrix();
    void clear();
    void identity();
    void set(const int irrep, const int firstIndex, const int secondIndex, const double value);
    double get(const int irrep, const int firstIndex, const int secondIndex) const;
    double* getBlock(const int irrep);
    double rmsDeviation(const Matrix const* other) const;
    static void write(const std::string fileName, const Indices* indices, double** storage);
    static void read(const std::string fileName, const int irrepsNumber, double** storage);
#ifdef CHEMPS2_MPI_COMPILATION
    void broadcast(const int root);
#endif
};