#include <iostream>
#include <igl/optical_flow.h>
#include <Eigen/Dense>
#include <fstream>

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
// ========================================================================================
// ========================================================================================
// ========================================================================================

int main(int argc, char *argv[])
{

  igl::OPTICALData o;
  o.method = read_option<int>("-m", argc, argv, "0"); // what method to use
  int demo_picture = read_option<int>("-d", argc, argv, "0");
  std::cout << "load image\n";
  if (demo_picture == 0)
  {
    igl::load_image1(o, TUTORIAL_SHARED_PATH "/box.0.bmp");
    igl::load_image2(o, TUTORIAL_SHARED_PATH "/box.1.bmp");
  }else{
    igl::load_image1(o, TUTORIAL_SHARED_PATH "/frame10.png");
    igl::load_image2(o, TUTORIAL_SHARED_PATH "/frame11.png");
  }

  o.alpha = 1.0;

  igl::solve_flow(o);
  for (int i = 0; i < 100; i++)
  {
    std::cout << "iteration " << i << "\n";
    igl::solve_flow(o);
  }

  std::ofstream u_file;
  u_file.open("u.txt");
  u_file << o.u;
  u_file.close();

  std::ofstream v_file;
  v_file.open("v.txt");
  v_file << o.v;
  v_file.close();
}