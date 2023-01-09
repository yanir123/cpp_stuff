#include <algorithm>
#include <vector>

template <class T>
T select_pivot_median(std::vector<T>& vec, const int left, const int right) {
  int middle = (left + right) / 2;
  if (vec[left] > vec[middle]) {
    std::swap(vec[left], vec[middle]);
  }
  if (vec[middle] > vec[right]) {
    std::swap(vec[middle], vec[right]);
  }
  if (vec[left] > vec[right]) {
    std::swap(vec[left], vec[right]);
  }

  std::swap(vec[middle], vec[right]);
  return vec[right];
}

template <class T>
int partition(std::vector<T>& vec, const T pivot, const int left,
              const int right) {
  int j = left - 1;
  for (int i = left; i <= right; i++) {
    if (vec[i] <= pivot) {
      j++;
      std::swap(vec[i], vec[j]);
    }
  }
  return j;
}

template <class T>
void quick_sort(std::vector<T>& vec, const int left, const int right) {
  if (left >= right) {
    return;
  }

  T pivot = select_pivot_median(vec, left, right);
  int q = partition(vec, pivot, left, right);

  quick_sort(vec, left, q - 1);
  quick_sort(vec, q + 1, right);
}

template <class T>
void quick_sort(std::vector<T>& vec) {
  int left = 0;
  int right = vec.size() - 1;
  quick_sort(vec, left, right);
}