/**
 * @file make_standalone.cc
 * @brief This file contains the main function for the standalone CAF generator.
 */
#include <iostream>
#include <vector>
#include <ctype.h>
#include "H5Cpp.h"

#include "include/record_fillers.h"
#include "include/products.h"
#include "include/event.h"
#include "include/reco_interaction.h"
#include "include/reco_particle.h"
#include "include/true_interaction.h"
#include "include/true_particle.h"

#include "duneanaobj/StandardRecord/StandardRecord.h"
#include "duneanaobj/StandardRecord/SRInteraction.h"
#include "duneanaobj/StandardRecord/SRRecoParticle.h"
#include "duneanaobj/StandardRecord/SRTrueInteraction.h"
#include "duneanaobj/StandardRecord/SRTrueParticle.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"

int main(int argc, char const *argv[])
{
    if(argc < 4)
    {
        std::cerr << "Usage: ./make_standalone <output file> <event offset> <input file(s)>" << std::endl;
        return 0;
    }

    TFile caf(argv[1], "recreate");
    TTree rec_tree("recTree", "Event tree for ML reconstruction");
    caf::StandardRecord *rec = new caf::StandardRecord;
    rec_tree.Branch("rec", &rec);

    TH1F * pot = new TH1F("TotalPOT", "TotalPOT", 1, 0, 1);
    TH1F * nevt = new TH1F("TotalEvents", "TotalEvents", 1, 0, 1);

    for(size_t n(3); n < argc; ++n)
    {
        H5::H5File file(argv[n], H5F_ACC_RDONLY);
        std::cout << "Opened file: " << argv[n] << std::endl;
        std::cout << "Calling get_all_events..." << std::endl;
        std::vector<dlp::types::Event> events(get_all_events(file));
        std::cout << "Found " << events.size() << " events. Starting event loop..." << std::endl;

        int i_evt = 0;
        for(dlp::types::Event &evt : events)
        {
            std::cout << "Processing event " << i_evt++ << "..." << std::endl;
            rec->common.ixn.spine.clear();
            rec->common.ixn.nspine = 0;
            rec->mc.nu.clear();
            rec->mc.nnu = 0;
            try
            {
                std::cout << "  Calling package_event..." << std::endl;
                package_event(rec, file, evt, std::atoi(argv[2]));
                std::cout << "  Calling get_product<RunInfo>..." << std::endl;
                std::vector<dlp::types::RunInfo> run_info(get_product<dlp::types::RunInfo>(file, evt));
                
                // For standalone FD CAFs, we default to FD HD
                rec->meta.fd_hd.run = run_info[0].run;
                rec->meta.fd_hd.subrun = run_info[0].subrun;
                rec->meta.fd_hd.event = run_info[0].event;
                rec->meta.fd_hd.enabled = true;
                
                pot->Fill(1);
                nevt->Fill(1);
                rec_tree.Fill();
            }
            catch(const H5::ReferenceException & e)
            {
                std::cerr << "Found incomplete entry for event." << std::endl;
            }
        }
        file.close();
    }
    pot->Write();
    nevt->Write();
    rec_tree.Write();
    caf.Close();

    return 0;
}
