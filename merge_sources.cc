/**
 * @file merge_sources.cc
 * @brief This file contains the main function for merging the outputs of the
 * standard reconstruction pathway with the SPINE reconstruction.
 */
#include <iostream>
#include <vector>
#include <tuple>
#include <ctype.h>
#include "H5Cpp.h"

#include "include/products.h"
#include "include/event.h"
#include "include/reco_interaction.h"
#include "include/reco_particle.h"
#include "include/true_interaction.h"
#include "include/true_particle.h"
#include "include/record_fillers.h"
#include "include/logger.h"

#include "duneanaobj/StandardRecord/StandardRecord.h"
#include "duneanaobj/StandardRecord/SRInteraction.h"
#include "duneanaobj/StandardRecord/SRRecoParticle.h"
#include "duneanaobj/StandardRecord/SRTrueInteraction.h"
#include "duneanaobj/StandardRecord/SRTrueParticle.h"

#include "TFile.h"
#include "TDirectoryFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TKey.h"

typedef std::tuple<size_t, size_t, size_t> index_t;

/**
 * @brief Helper to get run info from the first enabled detector in meta.
 */
index_t get_run_info(const caf::StandardRecord * rec) {
    if (rec->meta.fd_hd.enabled) return {rec->meta.fd_hd.run, rec->meta.fd_hd.subrun, rec->meta.fd_hd.event};
    if (rec->meta.fd_vd.enabled) return {rec->meta.fd_vd.run, rec->meta.fd_vd.subrun, rec->meta.fd_vd.event};
    if (rec->meta.nd_lar.enabled) return {rec->meta.nd_lar.run, rec->meta.nd_lar.subrun, rec->meta.nd_lar.event};
    if (rec->meta.nd_gar.enabled) return {rec->meta.nd_gar.run, rec->meta.nd_gar.subrun, rec->meta.nd_gar.event};
    if (rec->meta.tms.enabled) return {rec->meta.tms.run, rec->meta.tms.subrun, rec->meta.tms.event};
    if (rec->meta.sand.enabled) return {rec->meta.sand.run, rec->meta.sand.subrun, rec->meta.sand.event};
    if (rec->meta.lar2x2.enabled) return {rec->meta.lar2x2.run, rec->meta.lar2x2.subrun, rec->meta.lar2x2.event};
    if (rec->meta.minerva.enabled) return {rec->meta.minerva.run, rec->meta.minerva.subrun, rec->meta.minerva.event};
    if (rec->meta.pd_hd.enabled) return {rec->meta.pd_hd.run, rec->meta.pd_hd.subrun, rec->meta.pd_hd.event};
    
    // Default fallback (though this shouldn't happen in a valid CAF)
    return {rec->meta.fd_hd.run, rec->meta.fd_hd.subrun, rec->meta.fd_hd.event};
}

void copy_keyval_tree(TDirectory * output, TDirectory * input, const char * name)
{
    output->cd();
    if(input)
    {
        TDirectory * dir = output->mkdir(input->GetName());
        dir->cd();
        TTree * tree = static_cast<TTree *>(input->Get(name));
        TTree * new_tree = new TTree(name, tree->GetTitle());
        
        std::string key, value;
        std::string * pkey = &key;
        std::string * pvalue = &value;
        new_tree->Branch("key", &pkey);
        new_tree->Branch("value", &pvalue);
        tree->SetBranchAddress("key", &pkey);
        tree->SetBranchAddress("value", &pvalue);
        for(Long64_t i(0); i < tree->GetEntries(); ++i)
        {
            tree->GetEntry(i);
            new_tree->Fill();
        }
        new_tree->Write();
        output->cd();
    }
}

int main(int argc, char const * argv[])
{
    if(argc < 4)
    {
        LOG_CRITICAL("Usage: ./merge_sources <output_file> <input_caf_file> <input_h5_file>" << std::endl);
        return 0;
    }

    H5::H5File input_h5(argv[3], H5F_ACC_RDONLY);
    std::map<index_t, size_t> event_map;
    std::vector<dlp::types::Event> events(get_all_events(input_h5));
    for(size_t e(0); e < events.size(); ++e)
    {
        try
        {
            std::vector<dlp::types::RunInfo> run_info(get_product<dlp::types::RunInfo>(input_h5, events.at(e)));
            index_t index(run_info.back().run, run_info.back().subrun, run_info.back().event);
            event_map.insert(std::make_pair(index, e));
        }
        catch(const H5::ReferenceException & e)
        {
            LOG_CRITICAL("Found incomplete entry for event." << std::endl);
        }
    }

    TFile *input_caf = TFile::Open(argv[2], "read");
    if (!input_caf || input_caf->IsZombie())
    {
        LOG_CRITICAL("CAF file " << argv[2] << " could not be opened or does not exist." << std::endl);
        return 0;
    }
    TTree *input_tree = (TTree*)input_caf->Get("cafmaker/cafTree");
    if(!input_tree)
    {
        LOG_CRITICAL("CAF file does not contain recTree." << std::endl);
        return 0;
    }
    caf::StandardRecord *rec = new caf::StandardRecord;
    input_tree->SetBranchAddress("rec", &rec);

    TFile output_caf(argv[1], "recreate");
    TTree *output_tree = new TTree("recTree", "records");
    output_tree->Branch("rec", &rec);

    for(size_t n(0); n < input_tree->GetEntries(); ++n)
    {
        input_tree->GetEntry(n);
        rec->common.ixn.dlp.clear();
        rec->common.ixn.ndlp = 0;

        index_t index = get_run_info(rec);
        if(event_map.find(index) != event_map.end())
        {
            LOG_INFO("Matched Event " << std::get<2>(index) << " in (Run, Subrun) = (" << std::get<0>(index) << ", " << std::get<1>(index) << ") of CAF input to HDF5 event." << std::endl);
            try
            {
                package_event(rec, input_h5, events[event_map[index]]);
            }
            catch(const H5::ReferenceException & e)
            {
                LOG_CRITICAL("Found incomplete entry for event." << std::endl);
            }
        }
        output_tree->Fill();
    }

    TDirectoryFile * env = (TDirectoryFile *)input_caf->Get("env");
    copy_keyval_tree(&output_caf, env, "envtree");
    TH1D * total_pot = (TH1D*)input_caf->Get("TotalPOT");
    if (total_pot) total_pot->Write();
    TH1D * total_events = (TH1D*)input_caf->Get("TotalEvents");
    if (total_events) total_events->Write();
    output_tree->Write();
    TDirectoryFile * metadata = (TDirectoryFile *)input_caf->Get("metadata");
    copy_keyval_tree(&output_caf, metadata, "metatree");
    
    #ifdef MC_NOT_DATA
    TTree * global_tree = (TTree*)input_caf->Get("globalTree");
    if(global_tree)
    {
        output_caf.cd();
        TTree * clone = global_tree->CloneTree(-1, "fast");
        clone->Write();
        delete clone;
    }
    TTree * genie_tree = (TTree*)input_caf->Get("GenieEvtRecTree");
    if(genie_tree)
    {
        output_caf.cd();
        TTree * clone = genie_tree->CloneTree(-1, "fast");
        clone->Write();
        delete clone;
    }
    #endif
    
    input_h5.close();
    input_caf->Close();
    output_caf.Close();

    return 0;
}
