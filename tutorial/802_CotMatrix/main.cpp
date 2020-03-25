#include <iostream>
#include <igl/cotmatrix_numeric.h>
#include <Eigen/Dense>
#include <fstream>
#include <igl/readOBJ.h>
#include <igl/cotmatrix.h>
#include <igl/doublearea.h>
#include <igl/Timer.h>
#include <igl/squared_edge_lengths.h>

std::string mesh_name = "";
Eigen::MatrixXd V, U;
Eigen::MatrixXi F;
// ========================================================================================
// ========================================================================================
// ========================================================================================
template <class ValueType>
ValueType read_option(const char *option, int argc, char **argv, const char *default_value = nullptr);

template <>
std::string read_option<std::string>(const char *option, int argc, char **argv, const char *default_value)
{
  for (int i = 0; i < argc - 1; i++)
  {
    if (!strcmp(argv[i], option))
    {
      return std::string(argv[i + 1]);
    }
  }
  if (default_value)
    return std::string(default_value);
  std::cerr << "Option " << option << " was not provided. Exiting...\n";
  exit(1);
}
template <>
int read_option<int>(const char *option, int argc, char **argv, const char *default_value)
{
  return strtol(read_option<std::string>(option, argc, argv, default_value).c_str(), NULL, 10);
}
template <>
long read_option<long>(const char *option, int argc, char **argv, const char *default_value)
{
  return strtol(read_option<std::string>(option, argc, argv, default_value).c_str(), NULL, 10);
}
template <>
float read_option<float>(const char *option, int argc, char **argv, const char *default_value)
{
  return strtod(read_option<std::string>(option, argc, argv, default_value).c_str(), NULL);
}
template <>
double read_option<double>(const char *option, int argc, char **argv, const char *default_value)
{
  return strtof(read_option<std::string>(option, argc, argv, default_value).c_str(), NULL);
}

bool hasEnding(std::string const &fullString, std::string const &ending)
{
  if (fullString.length() >= ending.length())
  {
    return (0 == fullString.compare(fullString.length() - ending.length(), ending.length(), ending));
  }
  else
  {
    return false;
  }
}
// ========================================================================================
// ========================================================================================
// ========================================================================================

int main(int argc, char *argv[])
{
  Eigen::SparseMatrix<double> L;
  mesh_name = read_option<std::string>("-f", argc, argv, "");
  if (mesh_name == "")
  {
    mesh_name = TUTORIAL_SHARED_PATH "/armadillo.obj";
  }
  if (hasEnding(mesh_name, "obj") || hasEnding(mesh_name, "OBJ"))
  {
    igl::readOBJ(mesh_name, V, F);
  }
  igl::cotmatrix(V, F, L);
  igl::Timer t;
  t.start();
  for (int i = 0; i < 100; i++)
  {
    igl::cotmatrix(V, F, L);
  }
  t.stop();
  std::cout << "Standard Eigen method took: " << t.getElapsedTimeInMicroSec() / 100 << " microseconds on average\n";
  // std::cout << L.rows() << " " << L.cols() << " " << L.nonZeros() << "\n";

  igl::COTMATRIXData cd;
  cd.V = V;
  cd.F = F;
  igl::cotmatrix_numeric_intermediate(cd);

  // Eigen::MatrixXd l2;
  // igl::squared_edge_lengths(V, F, l2);
  // for (int i=0; i<l2.cols(); i++){
  //   for (int j = 0; j<l2.rows(); j++){
  //     std::cout<<l2(j, i)<<" "<<cd.datas[0][i*l2.rows()+j]<<"\n";
  //   }
  // }


  t.start();
  for (int i = 0; i < 100; i++)
  {
    igl::cotmatrix_numeric_intermediate(cd);
  }
  t.stop();
  std::cout << "Our method took: " << t.getElapsedTimeInMicroSec() / 100 << " microseconds on average\n";

  double err = 0;
  for (int i = 0; i < L.nonZeros(); i++)
  {
    err += std::pow(L.valuePtr()[i] - cd.result[i], 2);
  }

  std::cout<<"Error: "<<err<<"\n";
}