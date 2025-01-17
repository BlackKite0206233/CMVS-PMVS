
#include "mylapack.h"
#include <cstdlib>
#include <iostream>

// Use Eigen library or LAPACK
//#define PMVS_USE_LAPACK
#include <Eigen/Dense>

#if defined(PMVS_USE_LAPACK)
extern "C" {
#include <clapack/clapack.h>
#include <clapack/f2c.h>
};
#endif

using namespace std;

/*
// Solve Ax = 0.
// Values contain singular values stored in an increasing order
void Cmylapack::hlls(const std::vector<std::vector<float> >& A,
                     std::vector<float>& vec,
                     std::vector<float>& values) {
  char jobu = 'N';    char jobvt = 'S';
  integer M = (int)A.size();
  if (M == 0) {
    cerr << "Error in hlls" << endl;    exit (1);
  }
  integer N = (int)A[0].size();
  

  float *C = new float[M * N];
  int count = 0;
  for (int x = 0; x < N; ++x)
    for (int y = 0; y < M; ++y) {
      C[count++] = A[y][x];
    }
  integer lda = M;
  float *S = new float[M];
  integer LDU = 1;
  float *VT = new float[N * N];
  integer LDVT = N;
  integer lwork = 5 * max(N, M);
  float *work = new float[lwork];
  integer info;
  

  sgesvd_(&jobu, &jobvt, &M, &N, C, &lda, S, NULL, &LDU,
          VT, &LDVT, work, &lwork, &info);

  vec.resize(N);
  values.resize(N);
  for (int i = 0; i < N; ++i) {
    vec[i] = VT[N * (i + 1) - 1];
    values[i] = S[N - 1 - i];
  }
  delete [] work;
  delete [] VT;
  delete [] S;
  delete []C;
}

// Solve Ax = 0.
// Values contain singular values stored in an increasing order
void Cmylapack::hlls(const std::vector<std::vector<double> >& A,
                     std::vector<double>& vec,
                     std::vector<double>& values) {
  

  char jobu = 'N';    char jobvt = 'S';
  integer M = (int)A.size();
  if (M == 0) {
    cerr << "Error in hlls" << endl;    exit (1);
  }
  integer N = (int)A[0].size();
  

  double *C = new double[M * N];
  int count = 0;
  for (int x = 0; x < N; ++x)
    for (int y = 0; y < M; ++y) {
      C[count++] = A[y][x];
    }
  integer lda = M;
  double *S = new double[M];
  integer LDU = 1;
  double *VT = new double[N * N];
  integer LDVT = N;
  integer lwork = 5 * max(N, M);
  double *work = new double[lwork];
  integer info;
  

  dgesvd_(&jobu, &jobvt, &M, &N, C, &lda, S, NULL, &LDU,
          VT, &LDVT, work, &lwork, &info);

  vec.resize(N);
  values.resize(N);
  for (int i = 0; i < N; ++i) {
    vec[i] = VT[N * (i + 1) - 1];
    values[i] = S[N - 1 - i];
  }
  delete [] work;
  delete [] VT;
  delete [] S;
  delete []C;
}
*/

void Cmylapack::lls(const std::vector<std::vector<float>> &A, const std::vector<float> &b, std::vector<float> &ans) {
  int m = static_cast<int>(A.size());
  int n = static_cast<int>(A[0].size());

#if defined(PMVS_USE_LAPACK)
  char trans = 'N';
  integer nrhs = 1;
  vector<float> a;
  a.resize(m * n);

  int count = 0;
  for (int x = 0; x < n; ++x)
    for (int y = 0; y < m; ++y)
      a[count++] = A[y][x];
  integer lda = m;
  vector<float> b2;
  b2.resize(m);
  for (int i = 0; i < m; ++i)
    b2[i] = b[i];

  integer ldb = m;
  integer lwork = n + m;
  vector<float> work;
  work.resize(lwork);
  integer info;
  sgels_(&trans, &m, &n, &nrhs, &a[0], &lda, &b2[0], &ldb, &work[0], &lwork,
         &info);

  ans.resize(n);
  for (int i = 0; i < n; ++i)
    ans[i] = b2[i];
#else
  // Eigen implementation
  Eigen::MatrixXd matA(m, n);
  Eigen::VectorXd vecb(m);
  for (int x = 0; x < n; ++x)
    for (int y = 0; y < m; ++y)
      matA(y, x) = static_cast<double>(A[y][x]);
  for (int i = 0; i < m; ++i)
    vecb(i) = static_cast<double>(b[i]);
  Eigen::VectorXd vecx = matA.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(vecb);

  for (int i = 0; i < n; ++i)
    ans[i] = static_cast<float>(vecx(i));
#endif
}
/*
void Cmylapack::lls(const std::vector<std::vector<double> >& A,
                    const std::vector<double>& b,
                    std::vector<double>& ans) {
  char trans = 'N';
  integer m = (int)A.size();
  integer n = (int)A[0].size();
  integer nrhs = 1;
  vector<double> a;
  a.resize(m * n);

  int count = 0;
  for (int x = 0; x < n; ++x)
    for (int y = 0; y < m; ++y)
      a[count++] = A[y][x];
  integer lda = m;
  vector<double> b2;
  b2.resize(m);
  for (int i = 0; i < m; ++i)
    b2[i] = b[i];
  

  integer ldb = m;
  integer lwork = n + m;
  vector<double> work;
  work.resize(lwork);
  integer info;
  dgels_(&trans, &m, &n, &nrhs, &a[0], &lda, &b2[0], &ldb, &work[0],
         &lwork, &info);
  

  ans.resize(n);
  for (int i = 0; i < n; ++i)
    ans[i] = b2[i];
}

void Cmylapack::lls(std::vector<float>& A,
                    std::vector<float>& b,
                    integer width, integer height) {
  char trans = 'N';
  integer nrhs = 1;
  integer lwork = width * height;
  integer info;
  vector<float> work(width * height);
  

  sgels_(&trans, &width, &height, &nrhs, &A[0], &width, &b[0], &width, &work[0],
         &lwork, &info);
}

void Cmylapack::lls(std::vector<double>& A,
                    std::vector<double>& b,
                    integer width, integer height) {
  char trans = 'N';
  integer nrhs = 1;
  integer lwork = width * height;
  integer info;
  vector<double> work(width * height);
  

  dgels_(&trans, &width, &height, &nrhs, &A[0], &width, &b[0], &width, &work[0],
         &lwork, &info);
}

// SVD
void Cmylapack::svd(const std::vector<std::vector<float> >& A,
                    std::vector<std::vector<float> >& U,
                    std::vector<std::vector<float> >& VT,
                    std::vector<float>& S) {
  char jobu = 'A';    char jobvt = 'A';
  integer M = (int)A.size();
  if (M == 0) {
    cerr << "Error in hlls" << endl;    exit (1);
  }
  integer N = (int)A[0].size();
  

  float *C = new float[M * N];
  int count = 0;
  for (int x = 0; x < N; ++x)
    for (int y = 0; y < M; ++y) {
      C[count++] = A[y][x];
    }
  const integer minMN = min(M, N);
  integer lda = M;
  float *S2= new float[minMN];
  float *U2= new float[M * M];
  integer LDU = M;
  float *VT2= new float[N * N];
  integer LDVT = N;
  integer lwork = 5 * max(N, M);
  float *work= new float[lwork];
  integer info;
  

  sgesvd_(&jobu, &jobvt, &M, &N, C, &lda, S2, U2, &LDU,
          VT2, &LDVT, work, &lwork, &info);

  U.resize(M);
  for (int y = 0; y < M; ++y)
    U[y].resize(M);
  

  VT.resize(N);
  for (int y = 0; y < N; ++y)
    VT[y].resize(N);
  

  count = 0;
  for (int x = 0; x < M; ++x)
    for (int y = 0; y < M; ++y)
      U[y][x] = U2[count++];
  

  count = 0;
  for (int x = 0; x < N; ++x)
    for (int y = 0; y < N; ++y)
      VT[y][x] = VT2[count++];

  S.resize(minMN);
  for (int i = 0; i < minMN; ++i)
    S[i] = S2[i];

  delete [] work;
  delete [] VT2;
  delete [] U2;
  delete [] S2;
  delete [] C;
}
*/