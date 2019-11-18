//  aerofoil.h
//  IXP Code

#ifndef AEROFOIL_H_
#define AEROFOIL_H_

#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

#include "rotorsim.h"

class Aerofoil {
 public:
  Aerofoil(std::filesystem::path& datafile, std::size_t number = 0);
  Aerofoil(std::filesystem::path& datafile, std::filesystem::path meshfile,
           std::size_t number = 0);

  // Returns the difference between the highest and lowest Y coordinates
  double Thickness() const;
  
  // Returns the Korn technology factor of the aerofoil
  double Korn(Rotorsim* rotorsim, double lift) const;

  // Returns the drag divergence Mach number for a given lift coefficient
  double MachDivergence(Rotorsim* rotorsim, double Cl);
  
  // Returns the mesh file of the current aerofoil
  inline std::filesystem::path MeshFile() const;
  
  // Updated the mesh file associated to the current aerofoil
  void MeshFile(std::filesystem::path file);
  
  // Returns the number of coordinate points
  std::size_t SurfacePointCount() const { return SurfacePoints().size() / 2; };

  // Returns the vector of surface point coordinates
  std::vector<double> SurfacePoints() const;
  
  // Updates the vector of surface point coordinates to new values
  void SurfacePoints(std::vector<double> surface_points);
  
  // Copies the specified aerofoil parameters into the aerofoil memory
  void Load(bool surface_points = true, bool name = true, bool mesh = false);
  
  // Outputs the formatted aerofoil surface points into the specified file
  void Print(std::filesystem::path& file) const;

  // Returns the name of the aerofoil
  std::string Name() const;
  
  // Updates the name of the aerofoil
  void Name(std::string name);
  
  // Returns the numeric identifier of the aerofoil
  std::size_t Number() const;
  
  // Sets the numeric identifier of the aerofoil to a given value
  void Number(std::size_t number);
  
  // Calculates the required incidence to achieve the specified lift at the
  // specified Mach number. Uses the secant method to converge to a solution
  double Incidence(Rotorsim* rotorsim, double Mach, double lift) const;
  
  // Returns the value of the dCd/dMa gradient at a given incidence and Mach
  double DragToMach(Rotorsim* rotorsim, double Mach, double incidence) const;
  
  // Returns the value of the drag divergence Mach number for constant lift
  double DivergenceMach(Rotorsim* rotorsim, double lift) const;
  
 private:
  // Aerofoil properties
  std::string name_;
  std::size_t number_;
  std::filesystem::path datafile_;
  std::filesystem::path meshfile_;
  
  // Vector of surface point coordinates consisting of the vector of X
  // coordinates and the vector of Y coordinates stacked one after the other
  // in a single vector of lenght 2*N, where N is the number of surface points.
  // The first N points correspond to the X coordinates and the following N
  // points, to the corresponding Y coordinates (in order)
  std::vector<double> surface_points_;
  
  // Vector of mesh point coordinates consisting of the vector of X coordinates
  // and the vector of Y coordinates stacked one after the other in a single
  // vector of lenght 2*N, where N is the number of mesh points.
  // The first N points correspond to the X coordinates and the following N
  // points, to the corresponding Y coordinates (in order)
  std::vector<double> mesh_;
  
  // Reads the aerofoil coordinates datafile and returns the vector of
  // coordinates of the surface points
  std::vector<double> ReadSurfacePoints() const;
  
  // Reads the aerofoil coordinates datafile and returns the aerofoil name
  std::string ReadName() const;
  
  // Reads the aerofoil meshfile and returns the vector of coordinates of the
  // mesh points
  std::vector<double> ReadMesh() const;
};


#endif // AEROFOIL_H_
