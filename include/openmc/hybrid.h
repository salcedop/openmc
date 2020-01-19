//! \file nuclide.h
//! \brief Nuclide type and other associated types/data

#ifndef OPENMC_HYBRID_H
#define OPENMC_HYBRID_H

#include <hdf5.h>

#include "openmc/hybrid_reaction.h"

namespace openmc {

//==============================================================================
// Data for a nuclide
//==============================================================================

class hybrid {
public:
  // Constructors
  hybrid(hid_t group);


  std::vector<std::unique_ptr<HybridReaction>> HybrydReactions_; //!< Reactions
}

//! \param[in] file  HDF5 file object
void check_hybrid_version(hid_t file);

//! \brief Checks for the existence of a multipole library in the directory and
//! loads it
//!
//! \param[in] i_nuclide  Index in global nuclides array
void read_hybrid_data(int i_nuclide);



} //namespace

#endif // OPENMC_HYBRID_H
