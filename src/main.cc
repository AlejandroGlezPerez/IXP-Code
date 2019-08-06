//  main.cc
//  IXP Code

#include "aerofoil_library.h"
#include "rotorsim.h"

#include <chrono>
#include <filesystem>
#include <iostream>
#include <omp.h>

// Sample calculation - Korn Technology Factor of NACA 0012
int main() {
  omp_set_num_threads(6); // Deactivate hyperthreading
  
  // Start timer
  auto start = std::chrono::high_resolution_clock::now();
  
  // Load library, get NACA 0012 (#601 in library) and initialise flow solver
  AerofoilLibrary library{"library", 301};
  Aerofoil* naca0012 = library.Get(601);
  Rotorsim rotorsim;
  
  // Output aerofoil characteristics
  std::cout << "Name: " << naca0012->Name() << std::endl;
  std::cout << "Thickness: " << naca0012->Thickness() << std::endl;
  std::cout << "Lift: " << 0.74 << std::endl;
  
  // Calculate drag divergence and korn technology factor
  double korn;
  #pragma omp parallel shared(korn)
  #pragma omp master
    korn =  naca0012->Korn(&rotorsim, 0.74);
  std::cout << "Drag divergence: "
            << korn - 0.74 / 10.0 - naca0012->Thickness() << std::endl;
  std::cout << "Korn factor: " << korn << std::endl;
  
  // Stop timer and output time elapsed since start
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "Elapsed time: " << elapsed.count() << " s" << std::endl;
  
  return 0;
}
