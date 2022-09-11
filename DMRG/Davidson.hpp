#pragma once

class Davidson {
   private:
    int vectorLength;
    int multipicationsNumber;
    char state;
    bool debugPrint;
    char problemType;

    int maxVectorsNumber;
    int keepVectorsNumber;
    double diagCutoff;
    double residualTolerence;

    int vectorsNumber;
    double** vectors;
    double** hVectors;
    int numAllocated;

    double* mXMatrix;
    double* mXMatrixEigens;
    double* mXMatrixWork;
    double* mXMatrixVectors;
    double* mXMatrixRhs;
    int mXMatrixWorkLength;

    double* tVector;
    double* uVector;
    double* workVector;
    double* diag;
    double* rhs;

    double* ReorthoLowdin;
    double* ReorthoOverlapEigens;
    double* ReorthoOverlap;
    double* ReorthoEigensVectors;

    double FrobeniusNorm(double* currentVector);
    void SafetyCheckGuess();
    void AddNewVector();
    double DiagonalizeSmallMatrixAndCalcResidual();
    void CalculateNewVector();
    void Deflation();
    void MxMafterDeflation();
    void SolveLinearSystemDeflation(const int SOLUTIONS_NUMBER);

   public:
    Davidson(const int vectorLength, const int maxVectorsNumber, const int keepVectorsNumber, const double residualTolerence, const double diagCutoff, const bool debugPrint, const char problemType = 'E');
    virtual ~Davidson();
    char fetchInstruction(double** const pointers);
    int getNumMultiplications() const;
};