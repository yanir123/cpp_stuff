#ifndef ANN
#define ANN

#include <eigen3/Eigen/Eigen>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

typedef float Scalar;
typedef Eigen::MatrixXf Matrix;
typedef Eigen::RowVectorXf RowVector;
typedef Eigen::VectorXf ColVector;

class NeuralNetwork {
public:
    NeuralNetwork(std::vector<uint> topology, Scalar learningRate = Scalar(0.005));
    ~NeuralNetwork();
    void forwardProp(RowVector& input);
    void backProp(RowVector& output);
    void loss(RowVector& output);
    void updateWeights();
    void train(std::vector<RowVector*> input_data, std::vector<RowVector*> output_data);
    static Scalar activationFunction(Scalar x);
    static Scalar activationFunctionDerivative(Scalar x);
    
    std::vector<RowVector*> neuronLayers; 
    std::vector<RowVector*> cacheLayers; 
    std::vector<RowVector*> deltas; 
    std::vector<Matrix*> weights; 
    std::vector<uint> topology;
    Scalar learningRate;
};

void genData(std::string filename);
void ReadCSV(std::string filename, std::vector<RowVector*>& data);

#endif