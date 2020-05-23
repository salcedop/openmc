#include "openmc/hybrid_reaction.h"
#include "openmc/hybrid.h"
#include <string>
#include <utility> // for move

#include "openmc/constants.h"
#include "openmc/hdf5_interface.h"
#include "openmc/endf.h"

namespace openmc {

//==============================================================================
//Groupr Reaction implementation
//==============================================================================

HybridReaction::HybridReaction(hid_t rx_group)
{
    hid_t xs_dataset = open_dataset(rx_group,"groupr");
    read_attribute(rx_group,"start_point",xs_.threshold);
    
    // Read cross section values
    read_dataset(xs_dataset, xs_.value);
    close_dataset(xs_dataset);

}

} // namespace openmc
