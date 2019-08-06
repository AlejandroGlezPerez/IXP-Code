//  rotorsim.cc
//  IXP Code

#include "rotorsim.h"

#include <fstream>
#include <sstream>
#include <string>
#include <sys/wait.h>
#include <unistd.h>

#include <omp.h>

// Use noah Hypervisor to un Linux binaries on macOS
#if defined(__APPLE__) && defined(__MACH__)
  #define HYPERVISOR "/usr/local/bin/noah "
#elif defined(__gnu_linux__)
  #define HYPERVISOR ""
#endif

// Construct paths to file/folder locations
const std::filesystem::path Rotorsim::rotorsim_dir_{
  std::filesystem::current_path() / "rotorsim"};
const std::filesystem::path Rotorsim::mesher_app_{rotorsim_dir_.string()
                                                  + "/UOBSMB.x"};
const std::filesystem::path Rotorsim::preproc_app_{rotorsim_dir_.string()
                                                   + "/preproc.x"};
const std::filesystem::path Rotorsim::rotorsim_app_{rotorsim_dir_.string()
                                                    + "/rotorsim.x"};
const std::filesystem::path Rotorsim::mesher_template_{rotorsim_dir_.string()
                                                       + "/meshgen2D.conf"};
const std::filesystem::path Rotorsim::rotorsim_template_{rotorsim_dir_.string()
                                                         + "/rotorsim.conf"};

Rotorsim::Rotorsim() {
  std::filesystem::create_directory(rotorsim_dir_ / "run");
}

Rotorsim::~Rotorsim() {
  std::filesystem::remove_all(rotorsim_dir_ / "run");
}

std::vector<double> Rotorsim::RunSolver(const std::vector<double>& parameters,
                                         const std::filesystem::path& mesh) {
  // Get thread to run solver
  std::size_t thread_num = omp_get_thread_num();
  
  // Create directory to run solver
  std::filesystem::create_directory(RunDir(thread_num));
  
  // Run solver
  CopyMesh(mesh, thread_num);
  Configure(parameters, thread_num);
  Solve(thread_num);
  
  return ReadLoads(thread_num);
}

void Rotorsim::Configure(const std::vector<double>& parameters,
                         std::size_t thread_num) {
  // Open config file template for reading
  std::ifstream template_file{rotorsim_template_};
  
  // Open destination config file to save settings
  std::filesystem::path config_file_path{RunDir(thread_num).string()
                                         + "/rotorsim.conf"};
  std::ofstream config_file{config_file_path};
  
  // Substitute parameters in template file and copy to destination file
  std::string line;
  while (getline(template_file, line)) {
    // Locate Mach number placeholder and substitute
    std::size_t pos = line.find("_MACH_");
    if (pos != std::string::npos)
      line.replace(pos, strlen("_MACH_"), std::to_string(parameters[0]));
    
    // Locate angle of attack placeholder and substitute
    pos = line.find("_AOA_");
    if (pos != std::string::npos)
      line.replace(pos, strlen("_AOA_"), std::to_string(parameters[1]));
    config_file << line << std::endl;
  }
}

void Rotorsim::Solve(std::size_t thread, bool preprocess) {
  // Run as a different process to allow changing the current directory - this
  // is necessary for rotorsims to be able to read the configuration file
  pid_t child_pid;
  if ((child_pid = fork()) == 0) {
    std::ostringstream preproc_command_stream;
    preproc_command_stream << "cd " << RunDir(thread) << " && "
                           << HYPERVISOR << preproc_app_
                           << " >> preproc.out";
    std::string preproc_command = preproc_command_stream.str();
    std::system(preproc_command.c_str());
    
    // Run flow solver as a different process
    std::ostringstream rotorsim_command_stream;
    rotorsim_command_stream << "cd " << RunDir(thread) << " && "
                            << HYPERVISOR << rotorsim_app_
                            << " >> rotorsim.out";
    std::string rotorsim_command = rotorsim_command_stream.str();
    std::system(rotorsim_command.c_str());
    
    exit(0);
  }
  
  // Wait for child process to terminate
  int returnStatus;
  waitpid(child_pid, &returnStatus, 0);
}

void Rotorsim::CopyMesh(const std::filesystem::path& mesh,
                        std::size_t thread_num) {
  std::filesystem::path dest{RunDir(thread_num).string() + "/mesh.blk"};
  std::filesystem::copy_file(mesh, dest, std::filesystem::
                             copy_options::overwrite_existing);
}

inline std::filesystem::path Rotorsim::RunDir(std::size_t thread) {
  std::ostringstream run_dir_stream;
  run_dir_stream << rotorsim_dir_.string() << "/run/thread-" << thread;
  std::filesystem::path run_dir{run_dir_stream.str()};
  return run_dir;
}

std::vector<double> inline Rotorsim::ReadLoads(std::size_t thread) {
  // Throw exception if no results file found
  std::filesystem::path loads_file{RunDir(thread).string() + "/finalloads.dat"};
  
  // Store results
  double Cl, Cd, Cmy, dummy;
  std::ifstream loads_data{loads_file};
  loads_data >> Cl >> Cd >> dummy >> Cmy;
  std::vector<double> loads{Cl, Cd, Cmy};
  
  return loads;
}
