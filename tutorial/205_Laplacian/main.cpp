#include <igl/barycenter.h>
#include <igl/cotmatrix.h>
#include <igl/doublearea.h>
#include <igl/grad.h>
#include <igl/jet.h>
#include <igl/massmatrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/readDMAT.h>
#include <igl/readOFF.h>
#include <igl/readSTL.h>
#include <igl/repdiag.h>
// #include <igl/opengl/glfw/Viewer.h>
#include <igl/readOBJ.h>
#include <iostream>
#include "tutorial_shared_path.h"
#include <igl/writeOBJ.h>

#include <igl/massmatrix.h>
#include <igl/cot_smoothing.h>

Eigen::MatrixXd V, U;
Eigen::MatrixXi F;
Eigen::SparseMatrix<double> L;
// igl::opengl::glfw::Viewer viewer;

igl::COTSMOOTHData c;
// variables responsible for tests
std::string mesh_name = "";

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
  using namespace Eigen;
  using namespace std;

  c.method = read_option<int>("-m", argc, argv, "0");
  mesh_name = read_option<std::string>("-f", argc, argv, "");
  c.w = read_option<int>("-w", argc, argv, "10000000");
  if (mesh_name == "")
  {
    mesh_name = TUTORIAL_SHARED_PATH "/armadillo.obj";
  }
  else
  {
    mesh_name = TUTORIAL_SHARED_PATH "/" + mesh_name;
  }
  std::cout << mesh_name << "\n";
  // Load a mesh in OFF format
  if (hasEnding(mesh_name, "obj") || hasEnding(mesh_name, "OBJ"))
  {
    igl::readOBJ(mesh_name, V, F);
  }
  else if (hasEnding(mesh_name, "stl") || hasEnding(mesh_name, "STL"))
  {
    Eigen::MatrixXd N;
    igl::readSTL(mesh_name, V, F, N);
  }
  else if (hasEnding(mesh_name, "off") || hasEnding(mesh_name, "OFF"))
  {
    igl::readOFF(mesh_name, V, F);
  }

  c.V = V;
  c.F = F;

  const auto &key_down = [](unsigned char key, int mod) -> bool {
    switch (key)
    {
    case 'r':
    case 'R':
      U = V;
      break;
    case ' ':
    {
      igl::cot_smooth_solve(c);
      // U = c.smoothedV;
      break;
    }
    default:
      return false;
    }
    // Send new positions, update normals, recenter
    // viewer.data().set_vertices(U);
    // viewer.data().compute_normals();
    // viewer.core().align_camera_center(U, F);
    return true;
  };

  // Use original normals as pseudo-colors
  // MatrixXd N;
  // igl::per_vertex_normals(V, F, N);
  // MatrixXd C = N.rowwise().normalized().array() * 0.5 + 0.5;

  // Initialize smoothing with base mesh
  U = V;
  // viewer.data().set_mesh(U, F);

  // viewer.data().show_lines = false;

  // viewer.data().set_colors(C);
  // viewer.callback_key_down = key_down;

  key_down(' ', 0); // first iteration
  for (int i = 0; i < 100; i++)
  {
    key_down(' ', 0); // the following 100 iterations
    if ((i + 1) % 25 == 0)
      igl::writeOBJ("cot_smoothed_"+to_string(i+1)+".obj", c.smoothedV, c.F);
  }

  // cout << "Press [space] to smooth." << endl;
  // return viewer.launch();
}
