import ROOT

testFile = "/eos/experiment/fcc/ee/analyses/case-studies/bsm/LLPs/H_SS_4b/output_MadgraphPythiaDelphes/exoticHiggs_scalar_ms20GeV_sine-6.root"

#Mandatory: List of processes
processList = {

        #privately-produced signals
        'exoticHiggs_scalar_ms20GeV_sine-5':{},
        'exoticHiggs_scalar_ms20GeV_sine-6':{},
        'exoticHiggs_scalar_ms20GeV_sine-7':{},
        'exoticHiggs_scalar_ms60GeV_sine-5':{},
        'exoticHiggs_scalar_ms60GeV_sine-6':{},
        'exoticHiggs_scalar_ms60GeV_sine-7':{},

        # #centrally produced backgrounds
        # 'p8_ee_ZH_ecm240':{'fraction':0.01},
        # 'p8_ee_ZZ_ecm240':{'fraction':0.01},   
        # 'p8_ee_WW_ecm240':{'fraction':0.01},     
}

#Production tag. This points to the yaml files for getting sample statistics
#Mandatory when running over EDM4Hep centrally produced events
#Comment out when running over privately produced events
#prodTag     = "FCCee/spring2021/IDEA/"


#Input directory
#Comment out when running over centrally produced events
#Mandatory when running over privately produced events
inputDir = "/eos/experiment/fcc/ee/analyses/case-studies/bsm/LLPs/H_SS_4b/output_MadgraphPythiaDelphes"


#Optional: output directory, default is local dir
#outputDir = "/eos/experiment/fcc/ee/analyses/case-studies/bsm/LLPs/H_SS_4b/Reco_output_stage1/"
#outputDirEos = "/eos/experiment/fcc/ee/analyses/case-studies/bsm/LLPs/H_SS_4b/Reco_output_stage1/"
outputDir = "Vertex_output_stage1_seltracksOnlyd0_230612/"

#Optional: ncpus, default is 4
nCPUS       = 8

#Optional running on HTCondor, default is False
runBatch    = False
#runBatch    = True

#Optional batch queue name when running on HTCondor, default is workday
#batchQueue = "longlunch"

#Optional computing account when running on HTCondor, default is group_u_FCC.local_gen
#compGroup = "group_u_FCC.local_gen"

##-------- USER DEFINED CODE --------------------
# For costum displaced vertex selection, apply selection on the DVs with invariant mass higher than 1 GeV and distance from PV to DV less than 2000 mm, but longer than 4 mm
# and count the number of DVs that fulfill this selection
ROOT.gInterpreter.Declare("""
ROOT::VecOps::RVec<float> get_tracks_momentum(ROOT::VecOps::RVec<edm4hep::TrackState> tracks) {
    ROOT::VecOps::RVec<float> result;
    result.resize(tracks.size(), -1.);
    for (size_t i = 0; i < tracks.size(); ++i) {
        result[i] = FCCAnalyses::VertexingUtils::get_trackMom(tracks[i]);
    }
    return result;
}
""")

# method to get the corresponding MC indices from the track indices of the DVs
ROOT.gInterpreter.Declare("""
ROOT::VecOps::RVec<int> get_VertexTrack2MC(ROOT::VecOps::RVec<int> vertexRecoInd, ROOT::VecOps::RVec<int> RPMCindex) {
    ROOT::VecOps::RVec<int> result;
    result.resize(vertexRecoInd.size(),-1.);
    for (size_t i = 0; i < vertexRecoInd.size(); ++i) {
        result[i] = RPMCindex[vertexRecoInd[i]];
    }
    return result;
}
""")

# Method to truthmatch the DVs by comparing the DVs mcind to the mcind of the stable particles.
# Returns a score of matched tracks in DVs, i.e 0 = no matched DV tracks, 1 = one matched, 2 = two matched. 
# Important that the input indices are the MC indices, i.e the indices of the particles/tracks in the MC collection
ROOT.gInterpreter.Declare("""
double n_matchedDVtracks(ROOT::VecOps::RVec<int> MCtruthind, ROOT::VecOps::RVec<int> DVMCind) {
    double n_matchedDV_tracks = 0.;
    for (size_t i = 0; i < DVMCind.size(); ++i) {
        for (size_t j = 0; j < MCtruthind.size(); ++j) {
            if (DVMCind[i] == MCtruthind[j]) {
                n_matchedDV_tracks += 1.; 
            }
        }
    }
    return n_matchedDV_tracks;
}
""")

# Method to truthmatch the DVs by comparing the DVs mcind to the mcind of the stable particles.
# Returns a score of matched tracks in DVs, i.e 0 = no matched DV tracks, 1 = one matched, 2 = two matched. 
# Important that the input indices are the MC indices, i.e the indices of the particles/tracks in the MC collection
ROOT.gInterpreter.Declare("""
ROOT::VecOps::RVec<edm4hep::TrackState> matchedDVtracks(ROOT::VecOps::RVec<int> MCtruthind, ROOT::VecOps::RVec<int> DVsrecoInd, ROOT::VecOps::RVec<int> RPMCindex, const ROOT::VecOps::RVec<edm4hep::TrackState>& alltracks) {
    ROOT::VecOps::RVec<edm4hep::TrackState> matchedDV_tracks;
    ROOT::VecOps::RVec<int> DVMCind = get_VertexTrack2MC(DVsrecoInd, RPMCindex);
    for (size_t i = 0; i < DVMCind.size(); ++i) {
        //std::cout << "MC DV index: " << DVMCind[i] << std::endl;
        for (size_t j = 0; j < MCtruthind.size(); ++j) {
            if (DVMCind[i] == MCtruthind[j]) {
                //std::cout << "Matched!" << std::endl;
                matchedDV_tracks.push_back(alltracks[DVsrecoInd.at(i)]);
            }
        }
    }
    return matchedDV_tracks;
}
""")

# Get the absolute truthmatch score for all DVs compared to both of the true scalars
ROOT.gInterpreter.Declare("""
ROOT::VecOps::RVec<double> n_matchedDVtracks_all(ROOT::VecOps::RVec<int> MCtruthind_scalar1, ROOT::VecOps::RVec<int> MCtruthind_scalar2,
                                                ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> DVs,
                                                const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& reco,
                                                ROOT::VecOps::RVec<int> RPMCindex) {
    ROOT::VecOps::RVec<double> result;
    result.resize(DVs.size(),-1.);
    for (size_t i = 0; i < DVs.size(); ++i) {
        double n_matchedtracks = 0.;
        ROOT::VecOps::RVec<int> DVsRecoIndices;
        ROOT::VecOps::RVec<int> DVsMCIndices;
        DVsRecoIndices = FCCAnalyses::VertexingUtils::get_VertexRecoParticlesInd(DVs[i], reco);
        DVsMCIndices = get_VertexTrack2MC(DVsRecoIndices, RPMCindex);
        n_matchedtracks = n_matchedDVtracks(MCtruthind_scalar1, DVsMCIndices);
        if (n_matchedtracks == 0.) {
            n_matchedtracks = n_matchedDVtracks(MCtruthind_scalar2, DVsMCIndices);
        }
        result[i] = n_matchedtracks;
    }
    return result;
}
""")

# Get the purity of one DV
ROOT.gInterpreter.Declare("""
double truthmatchDV_purity(double n_matchedDVtracks, FCCAnalyses::VertexingUtils::FCCAnalysesVertex DV) {
    double result;
    double n_total_tracks_DV = FCCAnalyses::VertexingUtils::get_VertexNtrk(DV);
    result = n_matchedDVtracks/n_total_tracks_DV;
    return result;
}
""")

# Get the weighted purity of one DV
ROOT.gInterpreter.Declare("""
double truthmatchDV_purity_weighted(ROOT::VecOps::RVec<edm4hep::TrackState> truthmatched_DV_tracks,
                                    FCCAnalyses::VertexingUtils::FCCAnalysesVertex DV,
                                    const ROOT::VecOps::RVec<edm4hep::TrackState>& alltracks) {
    double result;
    ROOT::VecOps::RVec<edm4hep::TrackState> tracks_DV;

    ROOT::VecOps::RVec<float> tracks_truthmatched_DV_mom = get_tracks_momentum(truthmatched_DV_tracks);

    ROOT::VecOps::RVec<int> tracks_indices = DV.reco_ind;
    for (size_t k = 0; k < tracks_indices.size(); k++) {
        tracks_DV.push_back(alltracks[tracks_indices.at(k)]);
    }
    ROOT::VecOps::RVec<float> tracks_DV_mom = get_tracks_momentum(tracks_DV);

    double sum_truthmatched_track_mom = Sum(tracks_truthmatched_DV_mom, 0.);
    double sum_track_mom = Sum(tracks_DV_mom, 0.);
    //std::cout << "Sum momentum truthmatched: " << sum_truthmatched_track_mom << std::endl;
    //std::cout << "Sum momentum: " << sum_track_mom << std::endl;
    result = sum_truthmatched_track_mom/sum_track_mom;
    //std::cout << "Purity: " << result << std::endl;
    return result;
}
""")


# Get the purity of all DVs
ROOT.gInterpreter.Declare("""
ROOT::VecOps::RVec<double> truthmatchDV_purity_all(ROOT::VecOps::RVec<double> n_matchedDVtracks_all,
                                                ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> DVs) {
    ROOT::VecOps::RVec<double> result;
    result.resize(DVs.size(),-1.);
    for (size_t i = 0; i < DVs.size(); ++i) {
        double purity = 0.;
        double n_total_tracks_DV = FCCAnalyses::VertexingUtils::get_VertexNtrk(DVs[i]);
        double n_truthmatched_tracks_DV = n_matchedDVtracks_all[i];
        purity = n_truthmatched_tracks_DV/n_total_tracks_DV;
        result[i] = purity;
        }
    return result;
}
""")

# Get the purity of all DVs
ROOT.gInterpreter.Declare("""
ROOT::VecOps::RVec<double> truthmatchDV_purity_all_weighted(ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> DVs,
                                                            const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& reco,
                                                            ROOT::VecOps::RVec<int> MCtruthind1,
                                                            ROOT::VecOps::RVec<int> MCtruthind2,
                                                            ROOT::VecOps::RVec<int> RPMCindex,
                                                            const ROOT::VecOps::RVec<edm4hep::TrackState>& alltracks) {
    ROOT::VecOps::RVec<double> result;
    result.resize(DVs.size(),-1.);
    for (size_t i = 0; i < DVs.size(); ++i) {
        ROOT::VecOps::RVec<int> DVsRecoIndices = FCCAnalyses::VertexingUtils::get_VertexRecoParticlesInd(DVs[i], reco);

        ROOT::VecOps::RVec<edm4hep::TrackState> matchedDV_tracks = matchedDVtracks(MCtruthind1, DVsRecoIndices, RPMCindex, alltracks);
        if (matchedDV_tracks.size() == 0) {
            matchedDV_tracks = matchedDVtracks(MCtruthind2, DVsRecoIndices, RPMCindex, alltracks);
        }
        double purity = truthmatchDV_purity_weighted(matchedDV_tracks, DVs[i], alltracks);
        result[i] = purity;
    }
    //std::cout << "Purity in all method: " << result << std::endl;
    return result;
}
""")

# Get the acceptance of the true scalars, old way
ROOT.gInterpreter.Declare("""
float get_acceptance_trueparticles(ROOT::VecOps::RVec<float> decaylengths_scalar) {
    float result;
    float n_scalar_in_ID = 0.;
    for (size_t i = 0; i < decaylengths_scalar.size(); ++i) {
        if (decaylengths_scalar.at(i) < 2000) 
            n_scalar_in_ID += 1.;
    }
    result = n_scalar_in_ID/decaylengths_scalar.size();
    return result;
}
""")

# Get the acceptance of the true scalars, new way
ROOT.gInterpreter.Declare("""
ROOT::VecOps::RVec<float> get_acceptance_trueparticles_new(ROOT::VecOps::RVec<float> decaylengths_scalar) {
    ROOT::VecOps::RVec<float> result;
    for (size_t i = 0; i < decaylengths_scalar.size(); ++i) {
        if (decaylengths_scalar.at(i) < 2000) result[i] = 1.;
        else result[i] = 0.;
    }
    return result;
}
""")

# Count number of true decays in the decay volume
ROOT.gInterpreter.Declare("""
int get_n_trueparticles_vol(ROOT::VecOps::RVec<float> decaylengths_scalar) {
    int n_scalar_in_ID = 0;
    for (size_t i = 0; i < decaylengths_scalar.size(); ++i) {
        if (decaylengths_scalar.at(i) > 4 && decaylengths_scalar.at(i) < 2000) 
            n_scalar_in_ID += 1;
    }
    return n_scalar_in_ID;
}
""")

# Get the unique set of stable particles
ROOT.gInterpreter.Declare("""
ROOT::VecOps::RVec<int> get_list_of_FSparticles_indices(int index, ROOT::VecOps::RVec<edm4hep::MCParticleData> particle_collection, ROOT::VecOps::RVec<int> indices_collection) {
    ROOT::VecOps::RVec<int> result;
    std::vector<int> list_of_particles = FCCAnalyses::MCParticle::get_list_of_stable_particles_from_decay(index, particle_collection, indices_collection);
    std::set<int> unique_set_of_particles(list_of_particles.begin(), list_of_particles.end());
    for (auto & idx: unique_set_of_particles) {
        result.push_back(idx);
    }
    return result;
}
""")

# Truthmatch DVs by comparing the MC indices of the tracks from the DV to the MC indicies of the stable particles from a true scalar decay
ROOT.gInterpreter.Declare("""
ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> truthmatchDVs_vertex(ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> DVs,
                                                                                        const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& reco,
                                                                                        ROOT::VecOps::RVec<int> MCtruthind,
                                                                                        ROOT::VecOps::RVec<int> RPMCindex,
                                                                                        const ROOT::VecOps::RVec<edm4hep::TrackState>& alltracks) {
    ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> matchedDVs;
    for (size_t i = 0; i < DVs.size(); ++i) {
        ROOT::VecOps::RVec<int> DVsRecoIndices;
        ROOT::VecOps::RVec<int> DVsMCIndices;

        DVsRecoIndices = FCCAnalyses::VertexingUtils::get_VertexRecoParticlesInd(DVs[i], reco);

        ROOT::VecOps::RVec<edm4hep::TrackState> matchedDV_tracks = matchedDVtracks(MCtruthind, DVsRecoIndices, RPMCindex, alltracks);
        double purity = truthmatchDV_purity_weighted(matchedDV_tracks, DVs[i], alltracks);
        if (purity >= 0.5)   // "purity" threshold, i.e the fraction of truthmatched tracks of all tracks in the DV should be more than 0.5
            matchedDVs.push_back(DVs[i]);
    }
    return matchedDVs;
}
""")

# Get the fake rate of DVs in the event
ROOT.gInterpreter.Declare("""
float get_fake_rate(ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> truthmatchedDVs_scalar1, ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> truthmatchedDVs_scalar2, ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> DVs) {
    float n_DVs = DVs.size();
    
    if (n_DVs == 0) {
        return 0;
    }
    
    //std::cout << "number of DVs again: " << n_DVs << "   ";
    //std::cout << "number of truthmatched DVs 1: " << truthmatchedDVs_scalar1.size() << "   ";
    //std::cout << "number of truthmatched DVs 12: " << truthmatchedDVs_scalar2.size() << "   ";
    float nonmatchedDVs = n_DVs - truthmatchedDVs_scalar1.size() - truthmatchedDVs_scalar2.size();
    //std::cout << "number of nonmatched DVs: " << nonmatchedDVs << "   ";
    float fake_rate = nonmatchedDVs/n_DVs;
    return fake_rate;
}
""")

# Get the reconstruction effiency of one truthmatched DV to the true scalar, since the DVs are truthmatched to each scalar this method should be used two times
ROOT.gInterpreter.Declare("""
float get_reconstruction_efficiency(ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> truthmatchedDVs, 
                                    float decaylength_scalar) {
    float n_truthmatchedDVs = truthmatchedDVs.size();
    
    float truthdecays = 0.;
    if (decaylength_scalar < 2000) truthdecays += 1.;
    if (truthdecays == 0.) return 0.;

    float result = n_truthmatchedDVs/truthdecays;
    if (result > 1.) return 1.;     // since there could be several DVs matched to one scalar decay atm, this restricts the rec eff to not be higher than 1
    return result;
}
""")

ROOT.gInterpreter.Declare("""
std::vector<float> get_vector(float in1, float in2) {
    std::vector<float> result;
    result.push_back(in1);
    result.push_back(in2);
    return result;
}
""")

# get all DVs in a single vector
ROOT.gInterpreter.Declare("""
ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> get_all_DVs( ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> DV1,
							ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> DV2 ) {
  ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> result;
  for (auto &p:DV1){
    result.push_back(p);
  }
  for (auto &p:DV2){
    result.push_back(p);
  }
  return result;  
}
""")

# get the MC particles from the indices
ROOT.gInterpreter.Declare("""
ROOT::VecOps::RVec<edm4hep::MCParticleData> get_MCparticles(ROOT::VecOps::RVec<int> list_of_indices,  ROOT::VecOps::RVec<edm4hep::MCParticleData> in) {
    ROOT::VecOps::RVec<edm4hep::MCParticleData> result;
    for(size_t i = 0; i < list_of_indices.size(); ++i) {
        result.push_back(FCCAnalyses::MCParticle::sel_byIndex(list_of_indices[i], in));
    }
    return result;
}
""")

# get the sum of the invariant mass from the charged MC particles
ROOT.gInterpreter.Declare("""
double sum_charged_InvMass(ROOT::VecOps::RVec<edm4hep::MCParticleData> mc_in) {
    double result = 0.;
    ROOT::VecOps::RVec<float> all_charges = FCCAnalyses::MCParticle::get_charge(mc_in);
    ROOT::VecOps::RVec<TLorentzVector> all_InvMasses = FCCAnalyses::MCParticle::get_tlv(mc_in);

    int i = 0;
    for (auto & c: all_charges) {
        if (c != 0) {
            result += all_InvMasses.at(i).M();
        }
        i++;
    }
    return result;
}
""")

# seperate stable particles wrt their vertex position
ROOT.gInterpreter.Declare("""
ROOT::VecOps::RVec<int> get_scalar_FS_particles(ROOT::VecOps::RVec<int> list_of_scalar_decay_indices,
                                                ROOT::VecOps::RVec<edm4hep::MCParticleData> in,
                                                ROOT::VecOps::RVec<edm4hep::MCParticleData> scalar_MC) {
    ROOT::VecOps::RVec<int> result;
    ROOT::VecOps::RVec<edm4hep::MCParticleData> scalar_decays_MC = get_MCparticles(list_of_scalar_decay_indices, in);
    ROOT::VecOps::RVec<edm4hep::Vector3d> scalar_decays_vertices = FCCAnalyses::MCParticle::get_vertex(scalar_decays_MC);
    ROOT::VecOps::RVec<edm4hep::Vector3d> scalar_vertices = FCCAnalyses::MCParticle::get_vertex(scalar_MC);

    for(size_t i = 0; i < list_of_scalar_decay_indices.size(); ++i) {
        if (scalar_decays_vertices[i].x == scalar_vertices[i].x && scalar_decays_vertices[i].y == scalar_vertices[i].y && scalar_decays_vertices[i].z == scalar_vertices[i].z) {
            result.push_back(list_of_scalar_decay_indices[i]);
        }
    }

    return result;
}
""")

# get the endpoint like they do but modified to get it for only one MC particle
ROOT.gInterpreter.Declare("""
ROOT::VecOps::RVec<edm4hep::Vector3d> get_endPoint(ROOT::VecOps::RVec<edm4hep::MCParticleData> in, ROOT::VecOps::RVec<int> ind, ROOT::VecOps::RVec<edm4hep::MCParticleData> MC_particle)  {
        // ( carefull : if a Bs has oscillated into a Bsbar, this returns the production vertex of the Bsbar )
  ROOT::VecOps::RVec<edm4hep::Vector3d> result;
  for (auto & p: MC_particle) {
    edm4hep::Vector3d vertex(1e12, 1e12, 1e12);  // a default value for stable particles
    int db = p.daughters_begin ;
    int de = p.daughters_end;
    if (db != de) { // particle unstable
        int d1 = ind[db] ;   // first daughter
        if ( d1 >= 0 && d1 < in.size() ) {
            vertex = in.at(d1).vertex ;
        }
    }
    result.push_back(vertex);
  }
  return result;
}
""")



##-------- END USER DEFINED CODE -----------------

#Mandatory: RDFanalysis class where the use defines the operations on the TTree
class RDFanalysis():

    #__________________________________________________________
    #Mandatory: analysers funtion to define the analysers to process, please make sure you return the last dataframe, in this example it is df2
    def analysers(df):
        df2 = (
            df

            .Alias("Particle1", "Particle#1.index")
            .Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
            .Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")

            # MC event primary vertex
            .Define("MC_PrimaryVertex",  "FCCAnalyses::MCParticle::get_EventPrimaryVertex(21)( Particle )" )

            # number of tracks
            .Define("n_tracks","ReconstructedParticle2Track::getTK_n(EFlowTrack_1)")

            # Vertex fitting

            # First, reconstruct a vertex from all tracks 
            # Input parameters are 1 = primary vertex, EFlowTrack_1 contains all tracks, bool beamspotconstraint = true, bsc sigma x/y/z
            .Define("VertexObject_allTracks",  "VertexFitterSimple::VertexFitter_Tk ( 1, EFlowTrack_1, true, 4.5, 20e-3, 300)")

            # Select the tracks that are reconstructed  as primaries
            .Define("RecoedPrimaryTracks",  "VertexFitterSimple::get_PrimaryTracks( VertexObject_allTracks, EFlowTrack_1, true, 4.5, 20e-3, 300, 0., 0., 0., 0)")
            .Define("n_RecoedPrimaryTracks",  "ReconstructedParticle2Track::getTK_n( RecoedPrimaryTracks )")

            # the final primary vertex :
            .Define("PrimaryVertexObject",   "VertexFitterSimple::VertexFitter_Tk ( 1, RecoedPrimaryTracks, true, 4.5, 20e-3, 300) ")
            .Define("PrimaryVertex",   "VertexingUtils::get_VertexData( PrimaryVertexObject )")

            .Define("SecondaryTracks", "VertexFitterSimple::get_NonPrimaryTracks( EFlowTrack_1, RecoedPrimaryTracks)")

            # Displaced vertex reconstruction
            # no selection, using LCFI+ on all the tracks
            .Define("DV_no_sel", "VertexFinderLCFIPlus::get_SV_event(SecondaryTracks, EFlowTrack_1, PrimaryVertexObject, true, 9., 40., 5.)")
            # number of DVs no selection applied
            # number of DVs
            .Define('n_nosel_DVs', 'VertexingUtils::get_n_SV(DV_no_sel)')
            # invariant mass at the DVs (assuming the tracks to be pions)
            .Define('invMass_nosel_DVs', 'VertexingUtils::get_invM(DV_no_sel)')
            
            # select tracks with pT > 1 GeV
            .Define('sel_tracks_pt', 'VertexingUtils::sel_pt_tracks(1)(EFlowTrack_1)')
            # select tracks with |d0 |> 2 mm
            .Define('sel_tracks', 'VertexingUtils::sel_d0_tracks(2)(EFlowTrack_1)')        #sel_tracks_p
            # find the DVs
            .Define("DV_evt_seltracks", "VertexFinderLCFIPlus::get_SV_event(sel_tracks, EFlowTrack_1, PrimaryVertexObject, true, 9., 40., 5.)")
            # number of DVs
            .Define('n_seltracks_DVs', 'VertexingUtils::get_n_SV(DV_evt_seltracks)')
            # invariant mass at the DVs (assuming the tracks to be pions)
            .Define('invMass_seltracks_DVs', 'VertexingUtils::get_invM(DV_evt_seltracks)')

            # merge vertices with position within 10*error-of-position, get the tracks from the merged vertices and refit
            .Define('merged_DVs', 'VertexingUtils::mergeVertices(DV_evt_seltracks, EFlowTrack_1)')
            # number of merged DVs
            .Define("merged_DVs_n", "VertexingUtils::get_n_SV(merged_DVs)")
            # invariant mass at the merged DVs
            .Define('invMass_merged_DVs', 'VertexingUtils::get_invM(merged_DVs)')

            ## --- Truthmatching, reconstruction efficiency, fake rate, etc -----------------
            # returns a vector with the MC indices, like [index reco] = index MC
            .Define('RP_MC_index', "ReconstructedParticle2MC::getRP2MC_index(MCRecoAssociations0, MCRecoAssociations1, ReconstructedParticles)")

            # get the indices of the Higgs and the scalars, sorted in order h hs hs
            .Define('H2HSHS_indices', 'MCParticle::get_indices_ExclusiveDecay(25, {35, 35}, false, false)(Particle, Particle1)')

            # get the indices of the b quarks from related scalar, 1 are from the "first" scalar and 2 from the second one
            .Define('HS1_to_bb_indices', 'MCParticle::get_indices_MotherByIndex(H2HSHS_indices[1], {5, -5}, false, false, true, Particle, Particle1)')
            .Define('HS2_to_bb_indices', 'MCParticle::get_indices_MotherByIndex(H2HSHS_indices[2], {5, -5}, false, false, true, Particle, Particle1)')

            # method to retrieve the decay chain from each scalar to stable particles
            .Define('scalar1_decays_indices', 'get_list_of_FSparticles_indices(HS1_to_bb_indices[0], Particle, Particle1)')
            .Define('scalar2_decays_indices', 'get_list_of_FSparticles_indices(HS2_to_bb_indices[0], Particle, Particle1)')

            .Define('scalar1_decays_MC', 'get_MCparticles(scalar1_decays_indices, Particle)')
            .Define('scalar2_decays_MC', 'get_MCparticles(scalar2_decays_indices, Particle)')

            .Define('scalar1_decays_MC_charged_InvMass', 'sum_charged_InvMass(scalar1_decays_MC)')
            .Define('scalar2_decays_MC_charged_InvMass', 'sum_charged_InvMass(scalar2_decays_MC)')

            .Define('scalar1_decays_MC_vertices', 'MCParticle::get_vertex(scalar1_decays_MC)')
            .Define('scalar2_decays_MC_vertices', 'MCParticle::get_vertex(scalar2_decays_MC)')

            # .Define('scalar1_decays_MC_endPoints', 'MCParticle::get_endPoint(scalar1_decays_MC, Particle1)')
            # .Define('scalar2_decays_MC_endPoints', 'MCParticle::get_endPoint(scalar2_decays_MC, Particle1)')

            # .Define('scalar1_decays_MC_vertex_x', 'MCParticle::get_vertex_x(scalar1_decays_MC)')
            # .Define('scalar1_decays_MC_vertex_y', 'MCParticle::get_vertex_y(scalar1_decays_MC)')
            # .Define('scalar1_decays_MC_vertex_z', 'MCParticle::get_vertex_z(scalar1_decays_MC)')

            # .Define('scalar2_decays_MC_vertex_x', 'MCParticle::get_vertex_x(scalar2_decays_MC)')
            # .Define('scalar2_decays_MC_vertex_y', 'MCParticle::get_vertex_y(scalar2_decays_MC)')
            # .Define('scalar2_decays_MC_vertex_z', 'MCParticle::get_vertex_z(scalar2_decays_MC)')

            .Define('s1', 'myUtils::selMC_leg(0) (HS1_to_bb_indices, Particle)')
            .Define('s2', 'myUtils::selMC_leg(0) (HS2_to_bb_indices, Particle)')

            # get the b quarks, picking the first daughter from the list of decay particles from hs
            .Define('b1_s1', 'myUtils::selMC_leg(1) (HS1_to_bb_indices, Particle)')
            .Define('b1_s2', 'myUtils::selMC_leg(1) (HS2_to_bb_indices, Particle)')

            .Define('b2_s1', 'myUtils::selMC_leg(2) (HS1_to_bb_indices, Particle)')
            .Define('b2_s2', 'myUtils::selMC_leg(2) (HS2_to_bb_indices, Particle)')

            # get the production vertex of the b's, to retrieve info like position and radius
            # returns list of vertices, should check that the list only contains one vertex...
            .Define('bquarks_s1_vertices', 'MCParticle::get_vertex(b1_s1)')
            .Define('bquarks_s2_vertices', 'MCParticle::get_vertex(b2_s2)')

            # get the distance between the DVs reconstructed with track selection and the MC b production vertices (for both scalar 1 and 2)
            .Define('distance_b1_seltracks', 'VertexingUtils::get_d3d_SV_obj(DV_evt_seltracks, bquarks_s1_vertices.at(0))')
            .Define('distance_b2_seltracks', 'VertexingUtils::get_d3d_SV_obj(DV_evt_seltracks, bquarks_s2_vertices.at(0))')

            # ---- get the decay lengths of the scalars at generating level, code copied from MC anlysis -------
            # could be simplified/integrated with the code above, using get_d3d_... with the scalar vertices and the origin (0,0,0)

            # select generated scalar HS
            .Define('GenHS_PID',  'MCParticle::sel_pdgID(35, false)(Particle)')

            # get the production vertex for the 2 HS in x y z
            .Define('HS_vertex_x', 'MCParticle::get_vertex_x(GenHS_PID)')
            .Define('HS_vertex_y', 'MCParticle::get_vertex_y(GenHS_PID)')
            .Define('HS_vertex_z', 'MCParticle::get_vertex_z(GenHS_PID)')

            # # get the vertex position for the first group of b quarks
            .Define('s1b1_vertex_x', 'MCParticle::get_vertex_x(b1_s1)')
            .Define('s1b1_vertex_y', 'MCParticle::get_vertex_y(b1_s1)')
            .Define('s1b1_vertex_z', 'MCParticle::get_vertex_z(b1_s1)')

            .Define('s1b2_vertex_x', 'MCParticle::get_vertex_x(b2_s1)')
            .Define('s1b2_vertex_y', 'MCParticle::get_vertex_y(b2_s1)')
            .Define('s1b2_vertex_z', 'MCParticle::get_vertex_z(b2_s1)')

            # # get the vertex position for the second group of b quarks
            .Define('s2b1_vertex_x', 'MCParticle::get_vertex_x(b1_s2)')
            .Define('s2b1_vertex_y', 'MCParticle::get_vertex_y(b1_s2)')
            .Define('s2b1_vertex_z', 'MCParticle::get_vertex_z(b1_s2)')

            .Define('s2b2_vertex_x', 'MCParticle::get_vertex_x(b2_s2)')
            .Define('s2b2_vertex_y', 'MCParticle::get_vertex_y(b2_s2)')
            .Define('s2b2_vertex_z', 'MCParticle::get_vertex_z(b2_s2)')

            # get the decay length of HS 1
            .Define('decayLengthHS1', 'return sqrt((s1b1_vertex_x - HS_vertex_x.at(0))*(s1b1_vertex_x - HS_vertex_x.at(0)) + (s1b1_vertex_y - HS_vertex_y.at(0))*(s1b1_vertex_y - HS_vertex_y.at(0)) + (s1b1_vertex_z - HS_vertex_z.at(0))*(s1b1_vertex_z - HS_vertex_z.at(0)))')

            # get the decay length of HS 2
            .Define('decayLengthHS2', 'return sqrt((s2b1_vertex_x - HS_vertex_x.at(1))*(s2b1_vertex_x - HS_vertex_x.at(1)) + (s2b1_vertex_y - HS_vertex_y.at(1))*(s2b1_vertex_y - HS_vertex_y.at(1)) + (s2b1_vertex_z - HS_vertex_z.at(1))*(s2b1_vertex_z - HS_vertex_z.at(1)))')

            # get decay length of both scalars in a vector for each event
            .Define('decayLengthsHS', 'myUtils::get_both_scalars(decayLengthHS1, decayLengthHS2)')

            # Acceptance true particles, decays with decaylength < 2000 mm
            .Define('acceptance_true_scalars', 'get_acceptance_trueparticles(decayLengthsHS)')
            .Define('acceptance_true_scalars_new', 'get_acceptance_trueparticles_new(decayLengthsHS)')

            .Define('n_true_in_ID', 'get_n_trueparticles_vol(decayLengthsHS)')


            .Define('scalar1_position_decays_indices', 'get_scalar_FS_particles(scalar1_decays_indices, Particle, b1_s1)')
            .Define('scalar2_position_decays_indices', 'get_scalar_FS_particles(scalar2_decays_indices, Particle, b2_s2)')

            # # Truthmatch to the FS particles produced at the scalar decay position only
            # .Define('truthmatchedDVs_s1_position', 'truthmatchDVs_vertex(DV_evt_seltracks, ReconstructedParticles, scalar1_position_decays_indices, RP_MC_index, EFlowTrack_1)')
            # .Define('truthmatchedDVs_s2_position', 'truthmatchDVs_vertex(DV_evt_seltracks, ReconstructedParticles, scalar2_position_decays_indices, RP_MC_index, EFlowTrack_1)')
            # .Define('truthmatchedDVs_position', 'get_all_DVs(truthmatchedDVs_s1_position, truthmatchedDVs_s2_position)')

            # # number of truthmatched DVs
            # .Define('n_truthmatchedDVs_position', 'truthmatchedDVs_position.size()')
            
            # # Get the distance between the truthmatched DVs and the true scalar decays
            # .Define('distance_s1_truthmatchedDVs_position', 'VertexingUtils::get_d3d_SV_obj(truthmatchedDVs_s1_position, bquarks_s1_vertices.at(0))')
            # .Define('distance_s2_truthmatchedDVs_position', 'VertexingUtils::get_d3d_SV_obj(truthmatchedDVs_s2_position, bquarks_s2_vertices.at(0))')
            
            # # Get the purity of all DVs
            # .Define('truthmatch_purity_position', 'truthmatchDV_purity_all_weighted(DV_evt_seltracks, ReconstructedParticles, scalar1_position_decays_indices, scalar2_position_decays_indices, RP_MC_index, EFlowTrack_1)')
            
            # # Get the fake rate
            # .Define('fake_rate_position', 'get_fake_rate(truthmatchedDVs_s1_position, truthmatchedDVs_s2_position, DV_evt_seltracks)')

            # # Get the reconstruction efficiency for each scalar
            # .Define('rec_eff_position_s1', 'get_reconstruction_efficiency(truthmatchedDVs_s1_position, decayLengthHS1.at(0))')
            # .Define('rec_eff_position_s2', 'get_reconstruction_efficiency(truthmatchedDVs_s2_position, decayLengthHS2.at(0))')

            # # Get the total reconstruction efficiency, i.e the reconstruction efficiency for both scalars
            # .Define('rec_eff_position_total', 'get_vector(rec_eff_position_s1, rec_eff_position_s2)')



            # Truthmatch to the scalars seperatly
            .Define('truthmatchedDVs_nosel_s1', 'truthmatchDVs_vertex(DV_no_sel, ReconstructedParticles, scalar1_position_decays_indices, RP_MC_index, EFlowTrack_1)')
            .Define('truthmatchedDVs_nosel_s2', 'truthmatchDVs_vertex(DV_no_sel, ReconstructedParticles, scalar2_position_decays_indices, RP_MC_index, EFlowTrack_1)')
            .Define('truthmatchedDVs_nosel', 'get_all_DVs(truthmatchedDVs_nosel_s1, truthmatchedDVs_nosel_s2)')

            # number of truthmatched DVs
            .Define('n_truthmatchedDVs_nosel', 'truthmatchedDVs_nosel.size()')
            # invariant mass of truthmatched DVs
            .Define('invMass_truthmatched_nosel', 'VertexingUtils::get_invM(truthmatchedDVs_nosel)')
            
            # Get the distance between the truthmatched DVs and the true scalar decays
            .Define('distance_s1_truthmatchedDVs_nosel', 'VertexingUtils::get_d3d_SV_obj(truthmatchedDVs_nosel_s1, bquarks_s1_vertices.at(0))')
            .Define('distance_s2_truthmatchedDVs_nosel', 'VertexingUtils::get_d3d_SV_obj(truthmatchedDVs_nosel_s2, bquarks_s2_vertices.at(0))')

            # Get the absolute number of matched tracks of all DVs
            .Define('n_matchedDVtracks_nosel', 'n_matchedDVtracks_all(scalar1_position_decays_indices, scalar2_position_decays_indices, DV_no_sel, ReconstructedParticles, RP_MC_index)')
            
            # Get the purity of all DVs
            .Define('truthmatch_purity_nosel', 'truthmatchDV_purity_all_weighted(DV_no_sel, ReconstructedParticles, scalar1_position_decays_indices, scalar2_position_decays_indices, RP_MC_index, EFlowTrack_1)')
            
            # Get the fake rate
            .Define('fake_rate_nosel', 'get_fake_rate(truthmatchedDVs_nosel_s1, truthmatchedDVs_nosel_s2, DV_no_sel)')

            # Get the reconstruction efficiency for each scalar
            .Define('rec_eff_nosel_s1', 'get_reconstruction_efficiency(truthmatchedDVs_nosel_s1, decayLengthHS1.at(0))')
            .Define('rec_eff_nosel_s2', 'get_reconstruction_efficiency(truthmatchedDVs_nosel_s2, decayLengthHS2.at(0))')

            # Get the total reconstruction efficiency, i.e the reconstruction efficiency for both scalars
            .Define('rec_eff_nosel_total', 'get_vector(rec_eff_nosel_s1, rec_eff_nosel_s2)')

            
            # Truthmatch to the scalars seperatly
            .Define('truthmatchedDVs_seltracks_s1', 'truthmatchDVs_vertex(DV_evt_seltracks, ReconstructedParticles, scalar1_position_decays_indices, RP_MC_index, EFlowTrack_1)')
            .Define('truthmatchedDVs_seltracks_s2', 'truthmatchDVs_vertex(DV_evt_seltracks, ReconstructedParticles, scalar2_position_decays_indices, RP_MC_index, EFlowTrack_1)')
            .Define('truthmatchedDVs_seltracks', 'get_all_DVs(truthmatchedDVs_seltracks_s1, truthmatchedDVs_seltracks_s2)')

            # number of truthmatched DVs
            .Define('n_truthmatchedDVs_seltracks', 'truthmatchedDVs_seltracks.size()')
            # Invariant mass of truthmatched DVs
            .Define('invMass_truthmatched_seltracks', 'VertexingUtils::get_invM(truthmatchedDVs_seltracks)')
            
            # Get the distance between the truthmatched DVs and the true scalar decays
            .Define('distance_s1_truthmatchedDVs_seltracks', 'VertexingUtils::get_d3d_SV_obj(truthmatchedDVs_seltracks_s1, bquarks_s1_vertices.at(0))')
            .Define('distance_s2_truthmatchedDVs_seltracks', 'VertexingUtils::get_d3d_SV_obj(truthmatchedDVs_seltracks_s2, bquarks_s2_vertices.at(0))')

            # Get the absolute number of matched tracks of all DVs
            .Define('n_matchedDVtracks_seltracks', 'n_matchedDVtracks_all(scalar1_position_decays_indices, scalar2_position_decays_indices, DV_evt_seltracks, ReconstructedParticles, RP_MC_index)')
            
            # Get the purity of all DVs
            .Define('truthmatch_purity_seltracks', 'truthmatchDV_purity_all_weighted(DV_evt_seltracks, ReconstructedParticles, scalar1_position_decays_indices, scalar2_position_decays_indices, RP_MC_index, EFlowTrack_1)')
            
            # Get the fake rate
            .Define('fake_rate_seltracks', 'get_fake_rate(truthmatchedDVs_seltracks_s1, truthmatchedDVs_seltracks_s2, DV_evt_seltracks)')

            # Get the reconstruction efficiency for each scalar
            .Define('rec_eff_seltracks_s1', 'get_reconstruction_efficiency(truthmatchedDVs_seltracks_s1, decayLengthHS1.at(0))')
            .Define('rec_eff_seltracks_s2', 'get_reconstruction_efficiency(truthmatchedDVs_seltracks_s2, decayLengthHS2.at(0))')

            # Get the total reconstruction efficiency, i.e the reconstruction efficiency for both scalars
            .Define('rec_eff_seltracks_total', 'get_vector(rec_eff_seltracks_s1, rec_eff_seltracks_s2)')



            # Truthmatch to the scalars seperatly
            .Define('truthmatchedDVs_merged_s1', 'truthmatchDVs_vertex(merged_DVs, ReconstructedParticles, scalar1_position_decays_indices, RP_MC_index, EFlowTrack_1)')
            .Define('truthmatchedDVs_merged_s2', 'truthmatchDVs_vertex(merged_DVs, ReconstructedParticles, scalar2_position_decays_indices, RP_MC_index, EFlowTrack_1)')
            .Define('truthmatchedDVs_merged', 'get_all_DVs(truthmatchedDVs_merged_s1, truthmatchedDVs_merged_s2)')

            # number of truthmatched DVs
            .Define('n_truthmatchedDVs_merged', 'truthmatchedDVs_merged.size()')
            # Invariant mass of truthmatched DVs
            .Define('invMass_truthmatched_merged', 'VertexingUtils::get_invM(truthmatchedDVs_merged)')
            
            # Get the distance between the truthmatched DVs and the true scalar decays
            .Define('distance_s1_truthmatchedDVs_merged', 'VertexingUtils::get_d3d_SV_obj(truthmatchedDVs_merged_s1, bquarks_s1_vertices.at(0))')
            .Define('distance_s2_truthmatchedDVs_merged', 'VertexingUtils::get_d3d_SV_obj(truthmatchedDVs_merged_s2, bquarks_s2_vertices.at(0))')

            # Get the absolute number of matched tracks of all DVs
            .Define('n_matchedDVtracks_merged', 'n_matchedDVtracks_all(scalar1_position_decays_indices, scalar2_position_decays_indices, merged_DVs, ReconstructedParticles, RP_MC_index)')

            
            # Get the purity of all DVs
            .Define('truthmatch_purity_merged', 'truthmatchDV_purity_all_weighted(merged_DVs, ReconstructedParticles, scalar1_position_decays_indices, scalar2_position_decays_indices, RP_MC_index, EFlowTrack_1)')
            
            # Get the fake rate
            .Define('fake_rate_merged', 'get_fake_rate(truthmatchedDVs_merged_s1, truthmatchedDVs_merged_s2, merged_DVs)')

            # Get the reconstruction efficiency for each scalar
            .Define('rec_eff_merged_s1', 'get_reconstruction_efficiency(truthmatchedDVs_merged_s1, decayLengthHS1.at(0))')
            .Define('rec_eff_merged_s2', 'get_reconstruction_efficiency(truthmatchedDVs_merged_s2, decayLengthHS2.at(0))')

            # Get the total reconstruction efficiency, i.e the reconstruction efficiency for both scalars
            .Define('rec_eff_merged_total', 'get_vector(rec_eff_merged_s1, rec_eff_merged_s2)')
        )
        return df2


    #__________________________________________________________
    #Mandatory: output function, please make sure you return the branchlist as a python list
    def output():
        branchList = [
            # "n_tracks",
            # "n_RecoedPrimaryTracks",

            # 'n_nosel_DVs',
            # 'invMass_nosel_DVs',
            'n_seltracks_DVs',
            'invMass_seltracks_DVs',
            # "merged_DVs_n",
            # 'invMass_merged_DVs',

            # 'decayLengthHS1',
            # 'decayLengthHS2',
            'decayLengthsHS',

            # 'acceptance_true_scalars',
            # 'acceptance_true_scalars_new',

            # 'scalar1_decays_MC_charged_InvMass',
            # 'scalar2_decays_MC_charged_InvMass',
            # 'n_true_in_ID',

            # 's1b1_vertex_x',
            # 's1b1_vertex_y',
            # 's1b1_vertex_z',
            # 's1b2_vertex_x',
            # 's1b2_vertex_y',
            # 's1b2_vertex_z',

            # 's2b1_vertex_x',
            # 's2b1_vertex_y',
            # 's2b1_vertex_z',
            # 's2b2_vertex_x',
            # 's2b2_vertex_y',
            # 's2b2_vertex_z',

            'scalar1_decays_MC',
            'scalar2_decays_MC',

            # 'scalar1_decays_MC_vertex_x',
            # 'scalar1_decays_MC_vertex_y',
            # 'scalar1_decays_MC_vertex_z',
            # 'scalar2_decays_MC_vertex_x',
            # 'scalar2_decays_MC_vertex_y',
            # 'scalar2_decays_MC_vertex_z',

            # 'n_truthmatchedDVs_position',

            # 'distance_s1_truthmatchedDVs_position',
            # 'distance_s2_truthmatchedDVs_position',

            # 'truthmatch_purity_position',
            # 'fake_rate_position',
            # 'rec_eff_position_s1',
            # 'rec_eff_position_s2',
            # 'rec_eff_position_total',

            # 'n_truthmatchedDVs_nosel',
            # 'invMass_truthmatched_nosel',

            # 'distance_s1_truthmatchedDVs_nosel',
            # 'distance_s2_truthmatchedDVs_nosel',

            # 'n_matchedDVtracks_nosel',
            # 'truthmatch_purity_nosel',
            # 'fake_rate_nosel',
            # 'rec_eff_nosel_s1',
            # 'rec_eff_nosel_s2',
            # 'rec_eff_nosel_total',

            
            'n_truthmatchedDVs_seltracks',
            'invMass_truthmatched_seltracks',

            'distance_s1_truthmatchedDVs_seltracks',
            'distance_s2_truthmatchedDVs_seltracks',

            'n_matchedDVtracks_seltracks',
            'truthmatch_purity_seltracks',
            'fake_rate_seltracks',
            'rec_eff_seltracks_s1',
            'rec_eff_seltracks_s2',
            'rec_eff_seltracks_total',

            # 'n_truthmatchedDVs_merged',
            # 'invMass_truthmatched_merged',

            # 'distance_s1_truthmatchedDVs_merged',
            # 'distance_s2_truthmatchedDVs_merged',

            # 'n_matchedDVtracks_merged',
            # 'truthmatch_purity_merged',
            # 'fake_rate_merged',
            # 'rec_eff_merged_s1',
            # 'rec_eff_merged_s2',
            # 'rec_eff_merged_total',
        ]
        return branchList