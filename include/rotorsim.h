//  rotorsim.h
//  IXP Code

#ifndef ROTORSIM_H_
#define ROTORSIM_H_

#include <cstddef>
#include <filesystem>
#include <vector>

class Rotorsim {
public:
  Rotorsim();
  ~Rotorsim();
  
  // Executes rotorsims
  virtual std::vector<double> RunSolver(const std::vector<double>& parameters,
                                        const std::filesystem::path& mesh);

private:
  // Patb to directory containing the ROTORSIM flow solver
  static const std::filesystem::path rotorsim_dir_;
  
  // Path to UOBMSB mesher
  static const std::filesystem::path mesher_app_;
    
  // Path to rotorsim preprocessor and flow solver
  static const std::filesystem::path preproc_app_;
  static const std::filesystem::path rotorsim_app_;
  
  // Path to template configuration files
  static const std::filesystem::path mesher_template_;
  static const std::filesystem::path rotorsim_template_;
  
  // Set up configuration file to use the aerofoil
  void Configure(const std::vector<double>& parameters, std::size_t thread_num);

  // Run pre-processor and flow-solver on a new thread
  void Solve(std::size_t thread_num, bool preprocess = true);
  
  // Creates run directory for the current thread and copies the mesh file
  void CopyMesh(const std::filesystem::path& mesh, std::size_t thread_num);
  
  // Returns path to solver execution folder for specified thread
  inline std::filesystem::path RunDir(std::size_t thread_num);
  
  // Reads the loads file output by Rotorsim and returns the values as a vector
  std::vector<double> ReadLoads(std::size_t thread_num);
};

#endif // ROTORSIM_H_
