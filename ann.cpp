#include "ann.hpp"

NeuralNetwork::NeuralNetwork(std::vector<uint> topology, Scalar learningRate) {
    this->topology = topology;
    this->learningRate = learningRate;

    for (uint i = 0; i < this->topology.size(); i++) {
        if (i == this->topology.size() - 1) {
            this->neuronLayers.push_back(new RowVector(this->topology[i]));
        } else {
            this->neuronLayers.push_back(new RowVector(this->topology[i + 1]));
        }

        this->cacheLayers.push_back(new RowVector(this->neuronLayers.size()));
        this->deltas.push_back(new RowVector(this->neuronLayers.size()));

        if (i != this->topology.size() - 1) {
            this->neuronLayers.back()->coeffRef(this->topology[i]) = 1.0;
            this->cacheLayers.back()->coeffRef(this->topology[i]) = 1.0;
        }

        if (i > 0) {
            if (i != this->topology.size() - 1) {
                this->weights.push_back(new Matrix(this->topology[i - 1] + 1, this->topology[i] + 1));
                this->weights.back()->setRandom();
                this->weights.back()->col(this->topology[i]).setZero();
                this->weights.back()->coeffRef(this->topology[i - 1], this->topology[i]) = 1.0;
            } else {
                this->weights.push_back(new Matrix(this->topology[i - 1] + 1, topology[i]));
                this->weights.back()->setRandom();
            }
        }
    }
}

NeuralNetwork::~NeuralNetwork() {
    for (uint i = 0; i < this->neuronLayers.size(); i++) {
        delete this->neuronLayers[i];
    }

    for (uint i = 0; i < this->cacheLayers.size(); i++) {
        delete this->cacheLayers[i];
    }

    for (uint i = 0; i < this->deltas.size(); i++) {
        delete this->deltas[i];
    }

    for (uint i = 0; i < this->weights.size(); i++) {
        delete this->weights[i];
    }
}

void NeuralNetwork::forwardProp(RowVector& input) {
    this->neuronLayers.front()->block(0, 0, 1, this->neuronLayers.front()->size() - 1) = input;

    for (uint i = 0; i < this->topology.size(); i++) {
        (*this->neuronLayers[i]) = (*this->neuronLayers[i - 1]) * (*this->weights[i - 1]);
        this->neuronLayers[i]->block(0, 0, 1, this->topology[i]).unaryExpr(std::function<Scalar(Scalar)>(NeuralNetwork::activationFunction));
    }
}

void NeuralNetwork::loss(RowVector& output) {
    (*this->deltas.back()) = output - (*this->neuronLayers.back());

    for (uint i = this->topology.size() - 2; i > 0; i--) {
        (*this->deltas[i]) = (*this->deltas[i + 1]) * (this->weights[i]->transpose());
    }
}

void NeuralNetwork::updateWeights() {
    for (uint i = 0; i < this->topology.size() - 1; i++) {
        for (uint c = 0; c < this->weights[i]->cols() - (i != this->topology.size() - 2); c++) {
            for (uint r = 0; r < this->weights[i]->rows(); r++) {
                this->weights[i]->coeffRef(r, c) +=
                    this->learningRate *
                    this->deltas[i + 1]->coeffRef(c) *
                    this->activationFunctionDerivative(this->cacheLayers[i + 1]->coeffRef(c) * this->neuronLayers[i]->coeffRef(r));
            }
        }
    }
}

void NeuralNetwork::backProp(RowVector& output) {
    this->loss(output);
    this->updateWeights();
}

Scalar NeuralNetwork::activationFunction(Scalar x) {
    return tanhf(x);
}

Scalar NeuralNetwork::activationFunctionDerivative(Scalar x) {
    return 1 - tanhf(x) * tanhf(x);
}

void NeuralNetwork::train(std::vector<RowVector*> input_data, std::vector<RowVector*> output_data) {
    for(uint i = 0; i < input_data.size(); i++) {
        std::cout << "Input to neural network is : " << *input_data[i] << std::endl;
        this->forwardProp(*input_data[i]);
        std::cout << "Expected output is : " << *output_data[i] << std::endl;
        std::cout << "Output produced is : " << *neuronLayers.back() << std::endl;
        this->backProp(*output_data[i]);
        std::cout << "MSE : " << std::sqrt((*deltas.back()).dot((*deltas.back())) / deltas.back()->size()) << std::endl;
    }
}

void ReadCSV(std::string filename, std::vector<RowVector*>& data)
{
    data.clear();
    std::ifstream file(filename);
    std::string line, word;
    // determine number of columns in file
    getline(file, line, '\n');
    std::istringstream ss(line);
    std::vector<Scalar> parsed_vec;
    while (getline(ss, word, ',')) {
        parsed_vec.push_back(Scalar(std::stof(&word[0])));
    }
    uint cols = parsed_vec.size();
    data.push_back(new RowVector(cols));
    for (uint i = 0; i < cols; i++) {
        data.back()->coeffRef(1, i) = parsed_vec[i];
    }
 
    // read the file
    if (file.is_open()) {
        while (getline(file, line, '\n')) {
            std::stringstream ss(line);
            data.push_back(new RowVector(1, cols));
            uint i = 0;
            while (getline(ss, word, ',')) {
                data.back()->coeffRef(i) = Scalar(std::stof(&word[0]));
                i++;
            }
        }
    }
}

void genData(std::string filename)
{
    std::ofstream file1(filename + "-in");
    std::ofstream file2(filename + "-out");
    for (uint r = 0; r < 1000; r++) {
        Scalar x = rand() / Scalar(RAND_MAX);
        Scalar y = rand() / Scalar(RAND_MAX);
        file1 << x << "," << y << std::endl;
        file2 << 2 * x + 10 + y << std::endl;
    }
    file1.close();
    file2.close();
}