//  aerofoil.cc
//  IXP Code

#include "aerofoil.h"

#include <cmath>
#include <fstream>
#include <sstream>

// OpenMP
#include <omp.h>

// Mathematics implementation functions
namespace  {
  // Return sign of a value
  template <typename T> int sign(T value) {
    return (T(0) < value) - (value < T(0));
  }
}

Aerofoil::Aerofoil(std::filesystem::path& datafile) : datafile_{datafile} {}

Aerofoil::Aerofoil(std::filesystem::path& datafile,
                   std::filesystem::path meshfile)
    : datafile_{datafile}, meshfile_(meshfile) {}

double Aerofoil::Thickness() const {
  // Initialise iterator and set to beggining of Y-coordinate in array
  std::vector<double> surface_points = SurfacePoints();
  auto it = surface_points.begin();
  std::advance(it, surface_points.size() / 2);

  double ymax = 0.0;
  double ymin = 0.0;

  while(it != surface_points.end()) {
    if (*it > ymax)
      ymax = *it;
    else if (*it < ymin)
      ymin = *it;
    it++;
  }

  return ymax-ymin;
}

double Aerofoil::Korn(Rotorsim* rotorsim, double lift) const {
  return DivergenceMach(rotorsim, lift) + lift / 10.0 + Thickness();
}

std::filesystem::path Aerofoil::MeshFile() const {
  return meshfile_;
}

std::vector<double> Aerofoil::SurfacePoints() const {
  return surface_points_.empty() ? ReadSurfacePoints() : surface_points_;
}

std::string Aerofoil::Name() const {
  return name_.empty() ? ReadName() : name_;
}

std::vector<double> Aerofoil::ReadSurfacePoints() const {
  // Open aerofoil datafile for reading
  std::ifstream data{datafile_};

  // Skip aerofoil name
  data.ignore(std::numeric_limits<std::streamsize>::max(), data.widen('\n'));

  // Store number of datapoints specified in file
  std::size_t num_coords;
  data >> num_coords;

  // Read and insert surface coordinates at corresponding positions
  std::vector<double> surface_points;
  double x_coord, y_coord;
  for (std::size_t i = 0; data >> x_coord >> y_coord; i++) {
    surface_points.insert(surface_points.begin() + i, x_coord);
    surface_points.push_back(y_coord);
  }

  // Throw runtime error if number of surface points specified does not match
  // the number of surface points read
  if (num_coords != surface_points.size() / 2) {
    std::ostringstream err_stream;
    err_stream << "Error parsing coordinate file for aerofoil " << Name()
               << ". Number of surface points in file does not match "
               << "specified number." << std::endl
               << "Expected: " << num_coords << ". Found: "
               << surface_points.size() / 2 << ".";
    throw(std::runtime_error{err_stream.str()});
  }

  return surface_points;
}

std::string Aerofoil::ReadName() const {
  // Open aerofoil datafile for reading
  std::ifstream data{datafile_};

  // Read aerofoil name
  std::string line;
  std::getline(data, line);
  std::size_t name_start = line.find_first_of('\"') + 1;
  std::size_t name_end = line.find_last_of('\"');
  std::string name{line.substr(name_start, name_end - name_start)};

  return name;
}

std::vector<double> Aerofoil::ReadMesh() const {
  // Open aerofoil datafile for reading
  std::ifstream data{datafile_};

  // Skip aerofoil name
  data.ignore(std::numeric_limits<std::streamsize>::max(), data.widen('\n'));

  // Store number of datapoints specified in file
  std::size_t num_coords;
  data >> num_coords;

  double xcoord;
  double ycoord;

  // Read and insert surface coordinates at corresponding positions
  std::vector<double> surface_points;
  for (std::size_t i = 0; data >> xcoord >> ycoord; i++) {
    auto it = surface_points_.begin();
    advance (it, i);
    surface_points.insert(it, xcoord);
    surface_points.push_back(ycoord);
  }

  // Throw runtime error if number of surface points specified does not match
  // the number of surface points read
  if (num_coords != SurfacePointCount()) {
    std::ostringstream err_stream;
    err_stream << "Error parsing coordinate file for aerofoil " << Name()
    << ". Number of surface points in file does not match "
    << "specified number." << std::endl
    << "Expected: " << num_coords << ". Found: "
    << SurfacePointCount() << ".";
    throw(std::runtime_error{err_stream.str()});
  }


  return surface_points_;
}

void Aerofoil::Load(bool surface_points, bool name, bool mesh) {
  if (surface_points)
    surface_points_ = ReadSurfacePoints();

  if (name)
    name_ = ReadName();

  if (mesh)
    mesh_ = ReadMesh();
}

double Aerofoil::Incidence(Rotorsim* rotorsim, double Mach, double lift) const {
  // Set search parameters
  double constexpr kInitialIncidenceStep = 2.0;
  double constexpr kLiftTolerance = 0.01;
  std::size_t constexpr kMaxIterations = 5;

  // Use thin aerofoil theory for initial incidence guess
  std::vector<double> params_in{Mach, 0.0};
  double lift_in = rotorsim->RunSolver(params_in, MeshFile())[0];

  double dCl_in = lift - lift_in;
  double incidence_in = dCl_in / (2.0 * M_PI * M_PI / 180.0);

  // Calculate side solution points for gradient calculation
  double incidence_up = incidence_in + kInitialIncidenceStep;
  double dCl_up;
  #pragma omp task shared(dCl_up)
  {
    std::vector<double> params_up{Mach, incidence_up};
    double lift_up = rotorsim->RunSolver(params_up, MeshFile())[0];
    dCl_up = lift - lift_up;
  }
  double incidence_low = incidence_in - kInitialIncidenceStep;
  double dCl_low;
  #pragma omp task shared(dCl_low)
  {
    std::vector<double> params_low{Mach, incidence_low};
    double lift_low = rotorsim->RunSolver(params_low, MeshFile())[0];
    dCl_low = lift - lift_low;
  }
  #pragma omp taskwait

  // Use secant search to find incidence_in that solves dCl_in = 0
  for (std::size_t i = 0;
       fabs(dCl_in) > kLiftTolerance
       && i < kMaxIterations;
       i++) {
    // Calculate gradient and update central value
    double dCl_dAlpha = (dCl_up - dCl_low) / (incidence_up - incidence_low);
    incidence_in = incidence_up - dCl_up/dCl_dAlpha;
    std::vector<double> params_in{Mach, incidence_in};
    lift_in = rotorsim->RunSolver(params_in, MeshFile())[0];
    dCl_in = lift - lift_in;

    // Update search boundaries
    incidence_low = incidence_up;
    incidence_up = incidence_in;
    dCl_low = dCl_up;
    dCl_up = dCl_in;
  }

  return incidence_in;
};

double Aerofoil::DragToMach(Rotorsim* rotorsim, double Mach,
                            double incidence) const {
  // Set step for central difference
  constexpr double kMachStep = 1e-4;

  // Run flow solver asynchronously for specified parameters
  double drag_up;
  #pragma omp task shared(drag_up)
  {
    std::vector<double> params_up{Mach + kMachStep, incidence};
    drag_up = rotorsim->RunSolver(params_up, MeshFile())[1];
  }
  double drag_low;
  #pragma omp task shared(drag_low)
  {
    std::vector<double> params_low{Mach - kMachStep, incidence};
    drag_low = rotorsim->RunSolver(params_low, MeshFile())[1];
  }
  #pragma omp taskwait

  // Return central difference approximation
  return (drag_up - drag_low) / (2.0 * kMachStep);
}

double Aerofoil::DivergenceMach(Rotorsim* rotorsim, double lift) const {
  // Set search parameters
  double constexpr kGradientTolerance = 0.005;
  double constexpr kMachTolerance = 0.01;
  std::size_t constexpr kMaxIterations = 5;

  // Set Mach boundaries for initial estimate
  double Mach_up = 0.8;
  double Mach_low = 0.4;
  double Mach_in = (Mach_up + Mach_low) / 2.0;

  // Calculate initial incidence
  double incidence_up, incidence_low, incidence_in;
  #pragma omp task shared(incidence_up)
  incidence_up = Incidence(rotorsim, Mach_up, lift);
  #pragma omp task shared(incidence_low)
  incidence_low = Incidence(rotorsim, Mach_low, lift);
  #pragma omp task shared(incidence_in)
  incidence_in = Incidence(rotorsim, Mach_in, lift);
  #pragma omp taskwait

  // Calculate initial gradients - take away 0.1 to find zeroes of function
  // at dCd_dMach = 0.1
  double dCd_dMach_up, dCd_dMach_low, dCd_dMach_in;
  #pragma omp task shared(dCd_dMach_up)
  dCd_dMach_up = DragToMach(rotorsim, Mach_up, incidence_up) - 0.1;
  #pragma omp task shared(dCd_dMach_low)
  dCd_dMach_low = DragToMach(rotorsim, Mach_low, incidence_low) - 0.1;
  #pragma omp task shared(dCd_dMach_in)
  dCd_dMach_in = DragToMach(rotorsim, Mach_in, incidence_in) - 0.1;
  #pragma omp taskwait

  // Use bisection search to find Mach_in that solves dCd_dMach = 0.1
  for (std::size_t i = 0;
       (std::abs(dCd_dMach_in) > kGradientTolerance ||
           dCd_dMach_in != dCd_dMach_in) // continue when dCd_dMach_in = NaN
         && i < kMaxIterations
         && (Mach_low - Mach_up) < kMachTolerance;
       i++) {
    // Update search boundaries
    if (sign(dCd_dMach_in) == sign(dCd_dMach_low)) {
      Mach_low = Mach_in;
      dCd_dMach_low = dCd_dMach_in;
    } else {
      Mach_up = Mach_in;
      dCd_dMach_up = dCd_dMach_in;
    }

    // Update central value and repeat function evaluations
    Mach_in = (Mach_up + Mach_low) / 2.0;
    incidence_in = Incidence(rotorsim, Mach_in, lift);
    dCd_dMach_in = DragToMach(rotorsim, Mach_in, incidence_in) - 0.1;
  }

  return Mach_in;
};
