/**
 * @file record_fillers.h
 * @brief Declaration of various functions for filling the SPINE reconstruction
 * output objects for use within the CAF files.
 * @details This file contains the declaration of the functions that are
 * responsible for copying the SPINE reconstruction output data products from
 * the HDF5 files into the CAF objects. The functions are designed to be called
 * from either the standalone CAF generator or the merged CAF generator.
 * @author mueller@fnal.gov
 */
#ifndef RECORD_FILLERS_H
#define RECORD_FILLERS_H
#include <vector>
#include <map>
#include <ctype.h>

#include "event.h"
#include "reco_interaction.h"
#include "reco_particle.h"
#include "true_interaction.h"
#include "true_particle.h"

#include "duneanaobj/StandardRecord/StandardRecord.h"
#include "duneanaobj/StandardRecord/SRInteraction.h"
#include "duneanaobj/StandardRecord/SRRecoParticle.h"
#include "duneanaobj/StandardRecord/SRTrueInteraction.h"
#include "duneanaobj/StandardRecord/SRTrueParticle.h"

#define COPY(x,y) std::copy(y.begin(), y.end(), x)

// Forward declare.
template <class T>
std::vector<T> get_product(H5::H5File & file, dlp::types::Event & evt);

/**
 * @brief Constructs an instance of the caf::SRTrueParticle class from the
 * data in the dlp::types::TruthParticle data product.
 */
caf::SRTrueParticle fill_truth_particle(dlp::types::TruthParticle &p, uint64_t offset=0);

/**
 * @brief Constructs an instance of the caf::SRRecoParticle class from the data
 * in the dlp::types::RecoParticle data product.
 */
caf::SRRecoParticle fill_particle(dlp::types::RecoParticle &p, const std::map<int64_t, caf::TrueParticleID>& match_map, uint64_t offset=0);

/**
 * @brief Constructs an instance of the caf::SRTrueInteraction class from
 * the data in the dlp::types::TruthInteraction data product.
 */
caf::SRTrueInteraction fill_truth_interaction(dlp::types::TruthInteraction &in, std::vector<caf::SRTrueParticle> &particles, uint64_t offset=0);

/**
 * @brief Constructs an instance of the caf::SRInteraction class from the
 * data in the dlp::types::RecoInteraction data product.
 */
caf::SRInteraction fill_interaction(dlp::types::RecoInteraction &in, std::vector<caf::SRRecoParticle> &particles, const std::map<int64_t, int>& match_map, uint64_t offset=0);

/**
 * @brief Populates the StandardRecord object with the ML reconstruction
 * outputs.
 */
void package_event(caf::StandardRecord * rec, H5::H5File & file, dlp::types::Event & evt, uint64_t offset=0);

#endif
