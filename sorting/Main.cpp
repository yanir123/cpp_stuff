#include <iostream>

#include "QuickSort.hpp"

template <class T>
void print_array(std::vector<T> v) {
  std::cout << "{";

  for (size_t i = 0; i < v.size(); i++) {
    std::cout << v[i] << (i + 1 == v.size() ? "" : ", ");
  }

  std::cout << "}" << std::endl;
}

int main(int argc, char** argv, char** wenv) {
  std::vector<int> v{16, 4, 1, 20, 0, 3};

  quick_sort(v);

  print_array(v);

  return 0;
}