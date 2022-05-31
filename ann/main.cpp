#include "ann.hpp"

typedef std::vector<RowVector*> data;

int main(int argc, char** argv, char** env) {

    NeuralNetwork n({ 2, 3, 1 });
    data in_dat, out_dat;
    genData("test");
    ReadCSV("test-in", in_dat);
    ReadCSV("test-out", out_dat);
    n.train(in_dat, out_dat);

    for (uint i = 0; i < in_dat.size(); i++) {
        delete in_dat[i];
    }

    for (uint i = 0; i < out_dat.size(); i++) {
        delete out_dat[i];
    }

    return 0;
}