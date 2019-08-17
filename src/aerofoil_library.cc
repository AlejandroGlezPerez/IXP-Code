//  aerofoil_library.cc
//  IXP Code

#include "aerofoil_library.h"

#include <regex>
#include <sstream>

// Matrices and linear algebra
namespace {
  void append_coordinates(boost::numeric::ublas::matrix<double,
                          boost::numeric::ublas::column_major>& matrix,
                          std::vector<double>& coordinates) {
    matrix.resize(matrix.size1(), matrix.size2() + 1);
    auto last_column = matrix.end2() - 1;
    std::copy(coordinates.begin(), coordinates.end(), last_column.begin());
  }
}

AerofoilLibrary::AerofoilLibrary(std::filesystem::path library_directory,
                                 std::size_t num_coords)
    : directory_{std::filesystem::canonical(library_directory)},
      surface_directory_{directory_.string() + "/surfaces"},
      mesh_directory_{directory_.string() + "/meshes"},
      num_coords_{num_coords} {
        
  coordinates_.resize(num_coords * 2, 0);
        
  // Read all surface files in library, stop when next expected file not found
  for (std::size_t i = 0; i < kMaxAerofoils; i++) {
    std::filesystem::path surface_file = Datafile(i);
    if (!std::filesystem::exists(surface_file))
      break;
    Add(surface_file);
  }
}

void AerofoilLibrary::Add(std::filesystem::path& datafile) {
  // Disable addition if aerofoil count exceeds limit
  std::size_t aerofoil_num = aerofoils_.size();
  if (aerofoil_num > kMaxAerofoils)
    return;
  
  // Copy aerofoil to library if not already there
  std::filesystem::path dest = Datafile(aerofoil_num);
  if (datafile != Datafile(aerofoil_num))
    std::filesystem::copy(datafile, dest, std::filesystem::copy_options::
                          overwrite_existing);
  
  // Construct aerofoil with mesh if created
  std::filesystem::path meshfile = Meshfile(aerofoil_num);
  Aerofoil* aerofoil;
  if (std::filesystem::exists(meshfile))
      aerofoil = new Aerofoil{dest, meshfile};
  else
      aerofoil = new Aerofoil{dest};
  
  Add(aerofoil);
}

void AerofoilLibrary::Add(Aerofoil*& aerofoil) {
  // Disable addition if aerofoil count exceeds limit
  if (aerofoils_.size() > kMaxAerofoils)
    return;
  
  // Store surface points
  auto surface_points = aerofoil->SurfacePoints();
  
  // Fail if adding aerofoil with a difference number of surface points
  if (num_coords_ != surface_points.size() / 2) {
    std::ostringstream err_stream;
    err_stream << "Error adding aerofoil " << aerofoil->Name() << " to "
    "library. Number of aerofoil surface points does not match number of "
    "surface points expected by library." << std::endl
    << "Expected: " << num_coords_ << "."
    << "Found: "<< surface_points.size() / 2 << ".";
    throw(std::runtime_error{err_stream.str()});
  }
  append_coordinates(coordinates_, surface_points);
  aerofoils_.push_back(std::unique_ptr<Aerofoil>{aerofoil});
}

Aerofoil* AerofoilLibrary::Get(std::size_t aerofoil_num) const {
  return aerofoils_[aerofoil_num].get();
}

std::filesystem::path AerofoilLibrary::Datafile(std::size_t aerofoil_num) const {
  std::filesystem::path surface_file{surface_directory_.string() + "/"
                                     + std::to_string(aerofoil_num) + ".dat"};
  return surface_file;
};

std::filesystem::path AerofoilLibrary::Meshfile(std::size_t aerofoil_num)
    const {
  std::filesystem::path mesh_file{mesh_directory_.string() + "/"
    + std::to_string(aerofoil_num) + ".blk"};
  return mesh_file;
};
