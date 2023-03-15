
#include <Acts/Surfaces/DiscSurface.hpp>
#include <Acts/Surfaces/RadialBounds.hpp>

#include "TrackPropagationTest_processor.h"

#include <JANA/JApplication.h>
#include <JANA/JEvent.h>

#include <Math/GenVector/PxPyPzM4D.h>

#include <spdlog/spdlog.h>

#include <algorithms/tracking/ParticlesFromTrackFitResult.h>
#include <services/rootfile/RootFile_service.h>
#include <services/geometry/acts/ACTSGeo_service.h>

#include <edm4hep/MCParticle.h>
#include <edm4hep/ReconstructedParticle.h>
#include <edm4hep/SimTrackerHit.h>

#include <TVector3.h>


//------------------
// OccupancyAnalysis (Constructor)
//------------------
TrackPropagationTest_processor::TrackPropagationTest_processor(JApplication *app) :
	JEventProcessor(app)
{
}

//------------------
// Init
//------------------
void TrackPropagationTest_processor::Init()
{
    std::string plugin_name=("track_propagation_test");

    // Get JANA application
    auto app = GetApplication();

    // Ask service locator a file to write histograms to
    auto root_file_service = app->GetService<RootFile_service>();

    // Get TDirectory for histograms root file
    auto globalRootLock = app->GetService<JGlobalRootLock>();
    globalRootLock->acquire_write_lock();
    auto file = root_file_service->GetHistFile();
    globalRootLock->release_lock();

    // Create a directory for this plugin. And subdirectories for series of histograms
    m_dir_main = file->mkdir(plugin_name.c_str());

   //Define ntups
    mytup_part = new TNtuple("mytup_part","","mceta:mcp:mctheta:mcphi:receta:recp:rectheta:recphi");
    mytup_dirc_traj = new TNtuple("mytup_dirc_traj","","dirc_traj_x:dirc_traj_y:dirc_traj_z:dirc_traj_eta:dirc_traj_p:dirc_traj_theta:dirc_traj_phi");
    
    // Get log level from user parameter or default
    InitLogger(plugin_name);

    auto acts_service = GetApplication()->GetService<ACTSGeo_service>();

    m_propagation_algo.init(acts_service->actsGeoProvider(), logger());

    // Create HCal surface that will be used for propagation
    auto transform = Acts::Transform3::Identity();

    // make a reference disk for pfRICH
    const auto pfrichNZ = -1200.0;
    const auto pfrichNMinR = 46.0;
    const auto pfrichNMaxR = 630.0;
    auto pfrichNBounds = std::make_shared<Acts::RadialBounds>(pfrichNMinR, pfrichNMaxR);
    auto pfrichNTrf = transform * Acts::Translation3(Acts::Vector3(0, 0, pfrichNZ));
    m_pfrich_surface = Acts::Surface::makeShared<Acts::DiscSurface>(pfrichNTrf, pfrichNBounds);

  // make a reference cylinder to mimic DIRC
  const auto dirc_R = 700.0; //minr dirc
  const auto dirc_halfZ = 1700.0;
  // auto DIRC_center_Bounds =
      // std::make_shared<Acts::RadialBounds>(DIRC_center_MinR, DIRC_center_MaxR);
  auto dirc_center_Trf = transform * Acts::Translation3(Acts::Vector3(0, 0, 0));
  m_dirc_center_surface = Acts::Surface::makeShared<Acts::CylinderSurface>(dirc_center_Trf, dirc_R, dirc_halfZ);


}


//------------------
// Process
//------------------
// This function is called every event
void TrackPropagationTest_processor::Process(const std::shared_ptr<const JEvent>& event)
{
    m_log->trace("TrackPropagationTest_processor event");
    
    //Generated particles (only one particle generated)
    double mceta = 0; //particle eta
    double mctheta = 0; //particle phi
    double mcphi = 0; //particle phi
    double mcp = 0; //total momentum
    double mcpx = 0; //x momentum
    double mce = 0; //energy
    int num_primary = 0; //Number of primary particles

    auto mcParticles = event->Get<edm4hep::MCParticle>("MCParticles");
    for( size_t iParticle=0; iParticle<mcParticles.size(); iParticle++ )
      {
        auto mcparticle = mcParticles[iParticle];
        if(mcparticle->getGeneratorStatus() != 1) continue;
        auto& mom = mcparticle->getMomentum();
	TVector3 vec_mc_p(mom.x,mom.y,mom.z);
	mceta = vec_mc_p.Eta();
	mctheta = vec_mc_p.Theta();
	mcphi = vec_mc_p.Phi();
	mcp = vec_mc_p.Mag();
	mcpx = mom.x;
	num_primary++;
      }//End mcParticle loop

    //Reconstructed track momentum
    double recp = 0; //total momentum
    double receta = 0;//reconstructed eta
    double rectheta = 0;//reconstructed theta
    double recphi = 0;//reconstructed phi
    double recpx = 0;//x momentum
    int num_rec = 0; //Number of reconstructed particles

    auto RecParticles = event->Get<edm4eic::ReconstructedParticle>("ReconstructedChargedParticles");     
    for( size_t iParticle=0;iParticle<RecParticles.size(); iParticle++ )
      {
        auto Recparticle = RecParticles[iParticle];
        auto& mom = Recparticle->getMomentum();
	TVector3 vec_rec_p(mom.x,mom.y,mom.z);
	receta = vec_rec_p.Eta();
	rectheta = vec_rec_p.Theta();
	recphi = vec_rec_p.Phi();
	recp = vec_rec_p.Mag();
	recpx = mom.x;
        num_rec++;
      }
    // Get trajectories from tracking
    auto trajectories = event->Get<eicrecon::TrackingResultTrajectory>("CentralCKFTrajectories");
    //dirc points
    double dirc_traj_p = 0; //trajectory momentum
    double dirc_traj_eta = 0; //trajectory eta
    double dirc_traj_theta = 0; //trajectory theta
    double dirc_traj_phi = 0; //trajectory phi
    double dirc_traj_x = 0; //trajectory x
    double dirc_traj_y = 0; //trajectory y
    double dirc_traj_z = 0; //trajectory z

    // Iterate over trajectories
    m_log->debug("Propagating through {} trajectories", trajectories.size());
    for (size_t traj_index = 0; traj_index < trajectories.size(); traj_index++) {
        auto &trajectory = trajectories[traj_index];
        m_log->trace(" -- trajectory {} --", traj_index);

        edm4eic::TrackPoint* projection_point_dirc;
        try {
            // >>> try to propagate to surface <<<
            projection_point_dirc = m_propagation_algo.propagate(trajectory, m_dirc_center_surface);
        }
        catch(std::exception &e) {
            m_log->warn("Exception in underlying algorithm: {}. Trajectory is skipped", e.what());
        }

        if(!projection_point_dirc) {
            m_log->trace("   could not propagate!", traj_index);
            continue;
        }
        auto dirc_pos = projection_point_dirc->position;
        auto dirc_length =  projection_point_dirc->pathlength;
        auto dirc_traj_mom = projection_point_dirc->momentum;
	TVector3 vec_dirc_traj_p(dirc_traj_mom.x,dirc_traj_mom.y,dirc_traj_mom.z);
	TVector3 vec_dirc_traj_pos(dirc_pos.x,dirc_pos.y,dirc_pos.z);
        dirc_traj_p = vec_dirc_traj_p.Mag();
	dirc_traj_eta = vec_dirc_traj_p.Eta();
	dirc_traj_theta = vec_dirc_traj_p.Theta();
	dirc_traj_phi = vec_dirc_traj_p.Phi();
        dirc_traj_x = dirc_pos.x;
        dirc_traj_y = dirc_pos.y;
        dirc_traj_z = dirc_pos.z;

        m_log->trace("MC p,eta,th,ph: {:>10.2f} {:>10.2f} {:>10.2f} {:>10.2f}", mcp, mceta, mctheta,mcphi);
        m_log->trace("Rec p,eta,th,ph: {:>10.2f} {:>10.2f} {:>10.2f} {:>10.2f}", recp, receta, rectheta,recphi);
        m_log->trace("traj index, pos(x,y,z): {:>10} {:>10.2f} {:>10.2f} {:>10.2f}", traj_index, dirc_traj_x, dirc_traj_y, dirc_traj_z);
        m_log->trace("traj p,eta,th,ph: {:>10.2f} {:>10.2f} {:>10.2f} {:>10.2f}", dirc_traj_p, dirc_traj_eta, dirc_traj_theta,dirc_traj_phi);
	
	//Fill histograms
        if(num_primary==1 && num_rec==1)
          {
            mytup_part->Fill(mceta,mcp,mctheta,mcphi,receta,recp,rectheta,recphi);
            mytup_dirc_traj->Fill(dirc_traj_x,dirc_traj_y,dirc_traj_z,dirc_traj_eta,dirc_traj_p,dirc_traj_theta,dirc_traj_phi);
          }
    }
}


//------------------
// Finish
//------------------
void TrackPropagationTest_processor::Finish()
{
//    m_log->trace("TrackPropagationTest_processor finished\n");

}
