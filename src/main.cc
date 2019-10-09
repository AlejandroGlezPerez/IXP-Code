//  main.cc
//  IXP Code

#include "aerofoil_library.h"
#include "logger.h"
#include "rotorsim.h"

#include <chrono>
#include <filesystem>
#include <iostream>
#include <omp.h>

// Sample calculation - Korn Technology Factor of NACA 0012
int main() {
  
  LOG(info) << "*********************************************************\n";
  LOG(info) << "*              Starting SVD Optimiser v0.1              *\n";
  LOG(info) << "*          Written by Alejandro Gonzalez Perez          *\n";
  LOG(info) << "*           UoB Final Year Project 2019/2020            *\n";
  LOG(info) << "*********************************************************\n\n";
  
  omp_set_num_threads(6); // Deactivate hyperthreading
  LOG(info) << "Number of threads set to 6 (1 task per thread).";
  
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
  std::cout << "Number: " << naca0012->Number() << std::endl;
  
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
