//! \file groupr.h
//! MG Data for an incident neutron reaction

#ifndef OPENMC_GROUPR_REACTION_H
#define OPENMC_GROUPR_REACTION_H
#include <string>
#include <vector>
#include "hdf5.h"


namespace openmc {

//==============================================================================
//! MG Data for a single reaction including cross sections
//==============================================================================

class HybridReaction {
public:
  //! Construct reaction from HDF5 data
  //! \param[in] group HDF5 group containing reaction data
  explicit HybridReaction(hid_t rx_groupr);

  //! Cross section at a single temperature
  struct HybridXS {
    int threshold;
    std::vector<double> value;
  };

  HybridXS xs_;
};

} // namespace openmc

#endif // OPENMC_GROUPR_H
