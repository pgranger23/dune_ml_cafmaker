/**
 * @file record_fillers.cc
 * @brief Implementation of various functions for filling the SPINE
 * reconstruction output objects for use within the DUNE CAF files.
 */
#include <vector>
#include <ctype.h>
#include "H5Cpp.h"

#include "record_fillers.h"
#include "event.h"
#include "reco_interaction.h"
#include "reco_particle.h"
#include "true_interaction.h"
#include "true_particle.h"

#include "duneanaobj/StandardRecord/StandardRecord.h"
#include "logger.h"

caf::SRTrueParticle fill_truth_particle(dlp::types::TruthParticle &p, uint64_t offset)
{
    caf::SRTrueParticle part;
    part.pdg = p.pdg_code;
    // SRLorentzVector has E, px, py, pz members
    part.p.E = p.energy_init;
    part.p.px = p.momentum[0];
    part.p.py = p.momentum[1];
    part.p.pz = p.momentum[2];
    
    // SRTrueParticle has start_pos and end_pos
    part.start_pos = caf::SRVector3D(p.position[0], p.position[1], p.position[2]);
    part.end_pos = caf::SRVector3D(p.end_position[0], p.end_position[1], p.end_position[2]);
    
    part.interaction_id = p.interaction_id;
    part.G4ID = p.track_id;
    return part;
}

caf::SRRecoParticle fill_particle(dlp::types::RecoParticle &p, const std::map<int64_t, caf::TrueParticleID>& match_map, uint64_t offset)
{
    p.match_ids.reset(&p.match_ids_handle);
    p.match_overlaps.reset(&p.match_overlaps_handle);

    caf::SRRecoParticle part;
    
    // Core attributes
    part.primary = p.is_primary;
    part.pdg = p.pdg_code;
    part.E = p.ke / 1000.0; // MeV -> GeV
    part.p = caf::SRVector3D(p.momentum[0] / 1000.0, p.momentum[1] / 1000.0, p.momentum[2] / 1000.0); // MeV/c -> GeV/c
    part.start = caf::SRVector3D(p.start_point[0], p.start_point[1], p.start_point[2]);
    part.end = caf::SRVector3D(p.end_point[0], p.end_point[1], p.end_point[2]);
    part.contained = p.is_contained;

    // Generic extensions
    part.start_dir = caf::SRVector3D(p.start_dir[0], p.start_dir[1], p.start_dir[2]);
    part.end_dir = caf::SRVector3D(p.end_dir[0], p.end_dir[1], p.end_dir[2]);
    part.length = p.length;

    // Energy estimators (MeV -> GeV)
    part.energy_est.calo_ke = p.calo_ke / 1000.0;
    part.energy_est.csda_ke = p.csda_ke / 1000.0;
    part.energy_est.mcs_ke = p.mcs_ke / 1000.0;
    part.energy_est.depositions = p.depositions_sum / 1000.0;
    part.energy_est.mass = p.mass / 1000.0;
    part.energy_est.start_dedx = p.start_dedx; // kept in MeV/cm

    const std::vector<int> hypothesis_pdgs = {22, 11, 13, 211, 2212, 321};
    part.energy_est.ke_pdg = hypothesis_pdgs;
    for (int i = 0; i < 6; ++i) {
        part.energy_est.csda_ke_per_pid.push_back(p.csda_ke_per_pid[i] / 1000.0);
        part.energy_est.mcs_ke_per_pid.push_back(p.mcs_ke_per_pid[i] / 1000.0);
    }

    // PID probabilities & score arrays
    part.score = p.pid_scores[static_cast<int>(p.pid)]; // Best PID probability as score
    part.pid_prob.pid_pdg = hypothesis_pdgs;
    part.pid_prob.chi2_pdg = hypothesis_pdgs;
    for (int i = 0; i < 6; ++i) {
        part.pid_prob.pid_scores.push_back(p.pid_scores[i]);
        part.pid_prob.chi2_per_pid.push_back(p.chi2_per_pid[i]);
    }
    part.pid_prob.primary_scores.push_back(p.primary_scores[0]);
    part.pid_prob.primary_scores.push_back(p.primary_scores[1]);

    // Clustering
    part.cluster.size = p.size;
    part.cluster.num_fragments = p.num_fragments;

    // LArTPC-specific metadata
    part.lartpc.pid = static_cast<caf::LArTPCPid>(p.pid);
    part.lartpc.shape = static_cast<caf::LArTPCShape>(p.shape);
    part.lartpc.is_valid = p.is_valid;
    part.lartpc.is_cathode_crosser = p.is_cathode_crosser;
    part.lartpc.is_time_contained = p.is_time_contained;
    part.lartpc.start_straightness = p.start_straightness;
    part.lartpc.axial_spread = p.axial_spread;
    part.lartpc.directional_spread = p.directional_spread;
    part.lartpc.cathode_offset = p.cathode_offset;
    part.lartpc.vertex_distance = p.vertex_distance;

    // Truth matching
    for(size_t i = 0; i < p.match_ids.size(); ++i) {
        caf::TrueParticleID tid;
        int64_t match_id = p.match_ids[i];
        if (match_map.count(match_id) && match_map.at(match_id).ixn != -1) {
            tid = match_map.at(match_id);
        } else {
            tid.type = caf::TrueParticleID::kUnknown;
            tid.ixn = -1;
            tid.part = -1;
            LOG_WARN("Reco particle match_id " << match_id << " could not be resolved to a valid CAF TrueParticleID!" << std::endl);
        }
        part.truth.push_back(tid);
        part.truthOverlap.push_back(p.match_overlaps[i]);
        LOG_DEBUG("Reco particle matched to truth particle ID " << match_id << " (mapped to ixn " << tid.ixn << ", part " << tid.part << ") with overlap " << p.match_overlaps[i] << std::endl);
    }

    return part;
}

caf::SRTrueInteraction fill_truth_interaction(dlp::types::TruthInteraction &in, std::vector<caf::SRTrueParticle> &particles, uint64_t offset)
{
    in.particle_ids.reset(&in.particle_ids_handle);

    LOG_DEBUG("Filling truth interaction ID: " << in.id << ", PDG: " << in.pdg_code << ", E: " << in.energy_init << " MeV" << std::endl);

    caf::SRTrueInteraction ret;
    ret.id = in.id;
    ret.pdg = in.pdg_code;
    ret.E = in.energy_init;
    ret.vtx = caf::SRVector3D(in.vertex[0], in.vertex[1], in.vertex[2]);
    ret.momentum = caf::SRVector3D(in.momentum[0], in.momentum[1], in.momentum[2]);
    
    LOG_DEBUG("Interaction " << in.id << " has " << in.particle_ids.size() << " particles." << std::endl);
    for(int64_t id : in.particle_ids) {
        if (id < (int64_t)particles.size())
        {
            ret.prim.push_back(particles.at(id));
        }
        else
        {
            LOG_WARN("Particle ID " << id << " out of range (max " << particles.size() << ")" << std::endl);
        }
    }
    ret.nprim = ret.prim.size();

    return ret;
}

caf::SRInteraction fill_interaction(dlp::types::RecoInteraction &in, std::vector<caf::SRRecoParticle> &particles, const std::map<int64_t, int>& match_map, uint64_t offset)
{
    in.particle_ids.reset(&in.particle_ids_handle);
    in.match_ids.reset(&in.match_ids_handle);
    in.match_overlaps.reset(&in.match_overlaps_handle);

    LOG_DEBUG("Filling reco interaction ID: " << in.id << " with " << in.particle_ids.size() << " particles." << std::endl);

    caf::SRInteraction ret;
    ret.id = in.id;
    ret.vtx = caf::SRVector3D(in.vertex[0], in.vertex[1], in.vertex[2]);
    
    for(int64_t id : in.particle_ids) {
        if (id < (int64_t)particles.size())
            ret.part.dlp.push_back(particles.at(id));
    }
    ret.part.ndlp = ret.part.dlp.size();

    // Truth matching
    for(size_t i = 0; i < in.match_ids.size(); ++i) {
        int ixn_id = -1;
        int64_t match_id = in.match_ids[i];
        if (match_map.count(match_id)) {
            ixn_id = match_map.at(match_id);
        }
        
        if (ixn_id == -1) {
            LOG_WARN("Reco interaction match_id " << match_id << " could not be resolved to a valid CAF True Interaction index!" << std::endl);
        }
        
        ret.truth.push_back(ixn_id);
        ret.truthOverlap.push_back(in.match_overlaps[i]);
    }

    return ret;
}

void package_event(caf::StandardRecord * rec, H5::H5File & file, dlp::types::Event & evt, uint64_t offset)
{
    std::map<int64_t, caf::TrueParticleID> hdf5_part_to_caf_id;
    std::map<int64_t, int> hdf5_ixn_to_caf_idx;

    #ifdef MC_NOT_DATA
    LOG_DEBUG("Getting truth particles..." << std::endl);
    std::vector<dlp::types::TruthParticle> true_particles(get_product<dlp::types::TruthParticle>(file, evt));
    LOG_DEBUG("Found " << true_particles.size() << " truth particles." << std::endl);
    std::vector<caf::SRTrueParticle> caf_true_particles;
    for(dlp::types::TruthParticle &p : true_particles)
    {
        LOG_DEBUG("TruthParticle PDG: " << p.pdg_code << ", Energy: " << p.energy_init << " MeV, ID: " << p.id << std::endl);
        caf_true_particles.push_back(fill_truth_particle(p, offset));
    }

    LOG_DEBUG("Getting truth interactions..." << std::endl);
    std::vector<dlp::types::TruthInteraction> true_interactions(get_product<dlp::types::TruthInteraction>(file, evt));
    LOG_DEBUG("Found " << true_interactions.size() << " truth interactions." << std::endl);
    
    if (rec->mc.nu.empty()) {
        for(dlp::types::TruthInteraction &i : true_interactions)
            rec->mc.nu.push_back(fill_truth_interaction(i, caf_true_particles, offset));
        rec->mc.nnu = rec->mc.nu.size();
    }

    // MAP BUILDING: Create lookup map from G4ID to TrueParticleID for existing (or newly added) CAF truth
    std::map<int, caf::TrueParticleID> g4_to_caf_id;
    for (size_t ixn_idx = 0; ixn_idx < rec->mc.nu.size(); ++ixn_idx) {
        const auto& ixn = rec->mc.nu[ixn_idx];
        for (size_t p_idx = 0; p_idx < ixn.prim.size(); ++p_idx) {
            if (ixn.prim[p_idx].G4ID >= 0) {
                caf::TrueParticleID tid;
                tid.ixn = ixn_idx;
                tid.type = caf::TrueParticleID::kPrimary;
                tid.part = p_idx;
                g4_to_caf_id[ixn.prim[p_idx].G4ID] = tid;
            }
        }
        for (size_t p_idx = 0; p_idx < ixn.sec.size(); ++p_idx) {
            if (ixn.sec[p_idx].G4ID >= 0) {
                caf::TrueParticleID tid;
                tid.ixn = ixn_idx;
                tid.type = caf::TrueParticleID::kSecondary;
                tid.part = p_idx;
                g4_to_caf_id[ixn.sec[p_idx].G4ID] = tid;
            }
        }
    }

    if (get_global_log_level() <= LogLevel::DEBUG) {
        std::string g4_list = "";
        for (const auto& pair : g4_to_caf_id) {
            g4_list += std::to_string(pair.first) + " ";
        }
        LOG_DEBUG("Available G4IDs in CAF truth: " << g4_list << std::endl);
    }

    // Map HDF5 TruthParticle ID to valid TrueParticleID
    for (const auto& p : true_particles) {
        int g4id = p.track_id; // Using track_id as G4ID
        if (g4_to_caf_id.find(g4id) != g4_to_caf_id.end()) {
            hdf5_part_to_caf_id[p.id] = g4_to_caf_id[g4id];
        } else {
            caf::TrueParticleID tid;
            tid.ixn = -1;
            tid.type = caf::TrueParticleID::kUnknown;
            tid.part = -1;
            LOG_WARN("HDF5 TruthParticle ID " << p.id << " (G4ID: " << g4id << ") not found in existing CAF truth records!" << std::endl);
            hdf5_part_to_caf_id[p.id] = tid;
        }
    }

    // Map HDF5 TruthInteraction ID to valid CAF index
    for (auto& ixn : true_interactions) {
        int matched_caf_ixn = -1;
        ixn.particle_ids.reset(&ixn.particle_ids_handle);
        for(size_t i = 0; i < ixn.particle_ids.size(); i++) {
            int64_t part_id = ixn.particle_ids[i];
            if (hdf5_part_to_caf_id.count(part_id) && hdf5_part_to_caf_id[part_id].ixn != -1) {
                matched_caf_ixn = hdf5_part_to_caf_id[part_id].ixn;
                break;
            }
        }
        if (matched_caf_ixn != -1) {
            hdf5_ixn_to_caf_idx[ixn.id] = matched_caf_ixn;
        }
    }
    #endif

    LOG_DEBUG("Getting reco particles..." << std::endl);
    std::vector<dlp::types::RecoParticle> reco_particles(get_product<dlp::types::RecoParticle>(file, evt));
    LOG_DEBUG("Found " << reco_particles.size() << " reco particles." << std::endl);
    std::vector<caf::SRRecoParticle> caf_reco_particles;
    for(dlp::types::RecoParticle &p : reco_particles) {
        caf_reco_particles.push_back(fill_particle(p, hdf5_part_to_caf_id, offset));
    }

    LOG_DEBUG("Getting reco interactions..." << std::endl);
    std::vector<dlp::types::RecoInteraction> reco_interactions(get_product<dlp::types::RecoInteraction>(file, evt));
    
    rec->common.ixn.dlp.clear();
    for(dlp::types::RecoInteraction &i : reco_interactions)
        rec->common.ixn.dlp.push_back(fill_interaction(i, caf_reco_particles, hdf5_ixn_to_caf_idx, offset));
    rec->common.ixn.ndlp = rec->common.ixn.dlp.size();
}
