#include <vector>
using namespace std;

template<typename T>
void alloc2D (vector < T * > & array2, int r1, int r2) {
  for (int i = 0 ; i < r1 ; ++i) {
    T * temp1 = new T[r2];
    array2 . push_back(temp1);
  }
}

template<typename T>
void dealloc2D (vector < T * > & array2) {
  for (int i = 0 ; i < array2.size() ; ++i) {
    delete[] array2[i];
  }
  array2.clear();
}

template<typename T>
void alloc3D (vector < vector < T * > > & array3,
              int r1, int r2, int r3) {
  for (int i1 = 0 ; i1 < r1 ; ++i1) {
    vector<T *> temp2;
    for (int i2 = 0 ; i2 < r2 ; ++i2) {
      T * temp1 = new T[r3];
      temp2. push_back(temp1);
    }
    array3 . push_back(temp2);
  }
}

template<typename T>
void dealloc3D (vector < vector < T * > > & array3) {
  for (int i1 = 0 ; i1 < array3.size() ; ++i1) {
    for (int i2 = 0 ; i2 < array3[i1].size(); ++i2) {
      delete[] array3[i1][i2];
    }
    array3[i1].clear();
  }
  array3.clear();
}

template<typename T>
void alloc4D (vector < vector <vector < T * > > > & array4,
              int r1, int r2, int r3, int r4) {
  for (int i1 = 0 ; i1 < r1 ; ++i1) {
    vector < vector < T * > > temp3;
    for (int i2 = 0 ; i2 < r2 ; ++i2) {
      vector<T *> temp2;
      for (int i3 = 0 ; i3 < r3 ; ++i3) {
        T * temp1 = new T[r4];
        temp2. push_back(temp1);
      }
      temp3 . push_back(temp2);
    }
    array4 . push_back(temp3);
  }
}

template<typename T>
void dealloc4D (vector < vector < T * > > & array4) {
  for (int i1 = 0 ; i1 < array4.size() ; ++i1) {
    for (int i2 = 0 ; i2 < array4[i1].size(); ++i2) {
      for (int i3 = 0 ; i3 < array4[i1][i2].size(); ++i3) {
	delete[] array4[i1][i2][i3];
      }
      array4[i1][i2].clear();
    }
    array4[i1].clear();
  }
  array4.clear();
}
