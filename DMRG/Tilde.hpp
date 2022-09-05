#include "Indices.hpp"

class Tilde {
   private:
    Indices* handler;
    int* occupiedElements;
    double**** tildeMatrix;

   public:
    Tilde(Indices* handler);
    virtual ~Tilde();
    void clear();
    void set(const int firstIrrep, const int secondIrrep, const int first, const int second, const int third, const int fourth, const double value);
    double get(const int firstIrrep, const int secondIrrep, const int first, const int second, const int third, const int fourth) const;
    double * getBlock(const int firstIrrep, const int secondIrrep, const int first, const int third);
};