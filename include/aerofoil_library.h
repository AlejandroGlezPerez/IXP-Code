//  aerofoil_library.h
//  IXP Code

#ifndef AEROFOIL_LIBRARY_H_
#define AEROFOIL_LIBRARY_H_

#include "aerofoil.h"

#include <cstddef>
#include <filesystem>
#include <vector>

// Linear algebra and matrix class
#include <boost/numeric/ublas/matrix.hpp>

class AerofoilLibrary {
 public:
  AerofoilLibrary(std::filesystem::path library_directory,
                  std::size_t num_coords);

  // Adds a new aerofoil to the library specified by its data file. Throws a
  // runtime error if number of aerofoil surface points not the same as library
  // number of surface points
  void Add(std::filesystem::path& datafile);
  
  // Copies an existing aerofoil into the library. Aerofoil ownership is now
  // assumed by the library. Throws a runtime error if number of aerofoil
  // surface points not the same as number of surface points in library
  void Add(Aerofoil*& aerofoil);
  
  // Returns a pointer to the aerofoil with the specified number
  Aerofoil* Get(std::size_t aerofoil_num) const;
  
 private:
  // Maximum number of aerofoils in library
  static constexpr std::size_t kMaxAerofoils = 3000;
  
  // Paths to library directories
  const std::filesystem::path directory_;
  const std::filesystem::path surface_directory_;
  const std::filesystem::path mesh_directory_;
  
  // Storage containers for aerofoils in library
  std::vector<std::unique_ptr<Aerofoil>> aerofoils_;
  boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major>
      coordinates_;
  
  // Number of aerofoil surface points - constant for all aerofoils in library
  size_t num_coords_;
  
  // Constructs the library path to the datafile of the specified aerofoil
  std::filesystem::path Datafile(size_t aerofoil_num) const;
  
  // Constructs the library path to the mesh of the specified aerofoil
  std::filesystem::path Meshfile(size_t aerofoil_num) const;
  
};

#endif // AEROFOIL_LIBRARY_H_
