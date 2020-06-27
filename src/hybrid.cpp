#include "openmc/hybrid.h"


#include <string>
#include <utility> // for move

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/hdf5_interface.h"
#include "openmc/endf.h"
#include "openmc/nuclide.h"
#include "openmc/cross_sections.h"
#include "openmc/simulation.h"

namespace openmc {

namespace simulation{
xt::xtensor<double, 3> hybrid_tallies_;
}

//==============================================================================
//Groupr Reaction implementation
//==============================================================================

hybrid::hybrid(hid_t group) 
{
   // Read reactions
   // Need them in this order to match the score strings on the python/c api
  std::string depletion_reactions[7] = {"fission","(n,2n)","(n,3n)","(n,4n)","(n,p)","(n,a)","(n,gamma)"};
  std::string name_ = object_name(group).substr(1);
  hid_t rxs_group = open_group(group, "reactions");
  for (auto irx : depletion_reactions) {
      hid_t rx_group = open_group(rxs_group, irx.c_str());
      HybridReactions_.push_back(std::make_unique<HybridReaction>(rx_group));
      close_group(rx_group);
}
  close_group(rxs_group);
}


//========================================================================
// Non-member functions
//========================================================================

void check_hybrid_version(hid_t file_id)

{
  if (attribute_exists(file_id, "version")) {
    std::vector<int> version;
    read_attribute(file_id, "version", version);
    if (version[0] != HDF5_VERSION[0]) {
      fatal_error("HDF5 data format uses version " + std::to_string(version[0])
        + "." + std::to_string(version[1]) + " whereas your installation of "
        "OpenMC expects version " + std::to_string(HDF5_VERSION[0])
        + ".x data.");
    }
  } else {
    fatal_error("HDF5 data does not indicate a version. Your installation of "
      "OpenMC expects version " + std::to_string(HDF5_VERSION[0]) +
      ".x data.");
  }
}

void read_hybrid_data(int i_nuclide)
{
  // Look for hybrid data in cross_sections.xml
  const auto& nuc {data::nuclides[i_nuclide]};
  auto it = data::library_map.find({Library::Type::hybrid, nuc->name_});

  // If no hybrid library for this nuclide, just return
  if (it == data::library_map.end()) return;

  // Check if hybrid library exists
  int idx = it->second;
  std::string& filename = data::libraries[idx].path_;

  // Display message
  write_message("Reading " + nuc->name_ + " hybrid data from " + filename, 6);

  // Open file and make sure version is sufficient
  hid_t file = file_open(filename, 'r');
  check_hybrid_version(file);

  // Read nuclide data from HDF5
  hid_t group = open_group(file, nuc->name_.c_str());
  nuc->hybrid_ = std::make_unique<hybrid>(group);
  close_group(group);
  file_close(file);

}

//==============================================================================
// C-API function
//==============================================================================

//! \brief Returns a pointer to the hybrid tally results array along with its shape. This
//! allows a user to obtain in-memory results for hybri tallies directly in Python.
extern "C" int
openmc_hybrid_tally_results(double** h_results, size_t* h_shape)
{
  // Set pointer to results and copy shape
  *h_results = simulation::hybrid_tallies_.data();
  auto s = simulation::hybrid_tallies_.shape();
  h_shape[0] = s[0];
  h_shape[1] = s[1];
  h_shape[2] = s[2];
  return 0;
}

} // end namespace