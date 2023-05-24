import ROOT

testFile = "/eos/experiment/fcc/ee/analyses/case-studies/bsm/LLPs/H_SS_4b/output_MadgraphPythiaDelphes/exoticHiggs_scalar_ms20GeV_sine-6.root"

#Mandatory: List of processes
processList = {

        #privately-produced signals
        # 'exoticHiggs_scalar_ms20GeV_sine-5':{},
        'exoticHiggs_scalar_ms20GeV_sine-6':{},
        # 'exoticHiggs_scalar_ms20GeV_sine-7':{},
        # 'exoticHiggs_scalar_ms60GeV_sine-5':{},
        # 'exoticHiggs_scalar_ms60GeV_sine-6':{},
        # 'exoticHiggs_scalar_ms60GeV_sine-7':{},

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
inputDir = "/eos/experiment/fcc/ee/analyses/case-studies/bsm/LLPs/H_SS_4b/output_MadgraphPythiaDelphes/"


#Optional: output directory, default is local dir
#outputDir = "/eos/experiment/fcc/ee/analyses/case-studies/bsm/LLPs/H_SS_4b/Reco_output_stage1/"
#outputDirEos = "/eos/experiment/fcc/ee/analyses/case-studies/bsm/LLPs/H_SS_4b/Reco_output_stage1/"
outputDir = "Reco_output_stage1_230524/"

#Optional: ncpus, default is 4
nCPUS       = 8

#Optional running on HTCondor, default is False
runBatch    = False
#runBatch    = True

#Optional batch queue name when running on HTCondor, default is workday
#batchQueue = "longlunch"

#Optional computing account when running on HTCondor, default is group_u_FCC.local_gen
#compGroup = "group_u_FCC.local_gen"

#USER DEFINED CODE
# For costum displaced vertex selection, apply selection on the DVs with invariant mass higher than 1 GeV and distance from PV to DV less than 2000 mm, but longer than 4 mm
# and count the number of DVs that fulfill this selection
ROOT.gInterpreter.Declare("""
int filter_n_DVs(ROOT::VecOps::RVec<double> distanceDV, ROOT::VecOps::RVec<double> invMassDV) {
    int result = 0;
    for (size_t i = 0; i < invMassDV.size(); ++i) {
        // filter DVs with invariant mass > 1 GeV and distance from PV between 4 mm and 2000 mm
        if (invMassDV.at(i) > 1 && distanceDV.at(i) > 4 && distanceDV.at(i) < 2000)
            result += 1;
    }
    return result;
}

// method to get the corresponding MC indices from the track indices of the DVs
ROOT::VecOps::RVec<int> get_VertexTrack2MC(ROOT::VecOps::RVec<int> vertexRecoInd, ROOT::VecOps::RVec<int> RPMCindex) {
    ROOT::VecOps::RVec<int> result;
    result.resize(vertexRecoInd.size(),-1.);
    for (size_t i = 0; i < vertexRecoInd.size(); ++i) {
        result[i] = RPMCindex[vertexRecoInd[i]];
    }
    return result;
}

// method to truthmatch the DVs by comparing the DVs mcind to the mcind of the b quarks.
// Should be rewritten such that if one DV track is matched to one of the b quark pairs, it cannot belong to the other.
// Could also rewrite in more efficient list searching...
// Returns a score of matched tracks in DVs, i.e 0 = no matched DV tracks, 1 = one matched, 2 = two matched. 
// Can set a score/threshold for truthmatching
// Important that the input indices are the MC indices, i.e the indices of the particles/tracks in the MC collection

int truthmatchDV_score(ROOT::VecOps::RVec<int> MCtruthind, ROOT::VecOps::RVec<int> DVMCind) {
    int matchedDV_tracks = 0;

    for (size_t i = 0; i < DVMCind.size(); ++i) {
        for (size_t j = 0; j < MCtruthind.size(); ++j) {
            if (DVMCind[i] == MCtruthind[j]) {
                matchedDV_tracks += 1; 
            }
        }
    }
    return matchedDV_tracks;
}

ROOT::VecOps::RVec<int> truthmatchDV_score_all(ROOT::VecOps::RVec<int> MCtruthind1, ROOT::VecOps::RVec<int> MCtruthind2,
                                                ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> DVs,
                                                const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& reco,
                                                ROOT::VecOps::RVec<int> RPMCindex) {
    ROOT::VecOps::RVec<int> result;
    result.resize(DVs.size(),-1.);

    for (size_t i = 0; i < DVs.size(); ++i) {
        int truthmatch_score = 0;
        ROOT::VecOps::RVec<int> DVsRecoIndices;
        ROOT::VecOps::RVec<int> DVsMCIndices;
        DVsRecoIndices = FCCAnalyses::VertexingUtils::get_VertexRecoParticlesInd(DVs[i], reco);
        DVsMCIndices = get_VertexTrack2MC(DVsRecoIndices, RPMCindex);
        truthmatch_score = truthmatchDV_score(MCtruthind1, DVsMCIndices);
        if (truthmatch_score == 0) {
            truthmatch_score = truthmatchDV_score(MCtruthind2, DVsMCIndices);
        }
        result[i] = truthmatch_score;
    }
    return result;
}

ROOT::VecOps::RVec<double> truthmatchDV_purity_all(ROOT::VecOps::RVec<int> truthmatch_score_all,
                                                ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> DVs) {
    ROOT::VecOps::RVec<double> result;
    result.resize(DVs.size(),-1.);

    for (size_t i = 0; i < DVs.size(); ++i) {
        double purity = 0.;
        double n_total_tracks_DV = FCCAnalyses::VertexingUtils::get_VertexNtrk(DVs[i]);
        double n_truthmatched_tracks_DV = truthmatch_score_all[i];
        purity = n_truthmatched_tracks_DV/n_total_tracks_DV;
        result[i] = purity;
        }
    return result;
}

double truthmatchDV_purity(int truthmatch_score, FCCAnalyses::VertexingUtils::FCCAnalysesVertex DV) {
    double result;
    double n_total_tracks_DV = FCCAnalyses::VertexingUtils::get_VertexNtrk(DV);
    double n_truthmatched_tracks_DV = truthmatch_score;
    result = n_truthmatched_tracks_DV/n_total_tracks_DV;

    return result;
}

float get_acceptance_trueparticles(ROOT::VecOps::RVec<float> decaylengths_scalar) {
    float result;
    float n_scalar_in_ID = 0.;
    for (size_t i = 0; i < decaylengths_scalar.size(); ++i) {
        if (decaylengths_scalar.at(0) < 2000) 
            n_scalar_in_ID += 1.;
    }
    result = n_scalar_in_ID/decaylengths_scalar.size();
    return result;
}

// Important that the input indices are the MC indices, i.e the indices of the particles/tracks in the MC collection
ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> truthmatchDVs_event(ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> DVs,
                                                                                        const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& reco,
                                                                                        ROOT::VecOps::RVec<int> MCtruthind1, ROOT::VecOps::RVec<int> MCtruthind2,
                                                                                        ROOT::VecOps::RVec<int> RPMCindex) {

    ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> matchedDVs;

    for (size_t i = 0; i < DVs.size(); ++i) {
        ROOT::VecOps::RVec<int> DVsRecoIndices;
        ROOT::VecOps::RVec<int> DVsMCIndices;
        int truthmatch_score_1 = 0;
        int truthmatch_score_2 = 0;
        DVsRecoIndices = FCCAnalyses::VertexingUtils::get_VertexRecoParticlesInd(DVs[i], reco);
        // std::cout << "Index in DV list: " << i << " and DV track indices " << DVsRecoIndices;
        DVsMCIndices = get_VertexTrack2MC(DVsRecoIndices, RPMCindex);
        truthmatch_score_1 = truthmatchDV_score(MCtruthind1, DVsMCIndices);
        truthmatch_score_2 = truthmatchDV_score(MCtruthind2, DVsMCIndices);

        double n_truthmatched_tracks_DV = truthmatch_score_1 + truthmatch_score_2;

        // change threshold here. Either put an absolute one where the threshold is set requiring a minimum number of truth matched tracks
        // or have a "purity" threshold, i.e the fraction of truthmatched tracks of all tracks in the DV should be more than e.g 0.5
        double n_total_tracks_DV = FCCAnalyses::VertexingUtils::get_VertexNtrk(DVs[i]);
        double purity = n_truthmatched_tracks_DV/n_total_tracks_DV;
        if (purity > 0.5)   //absolute threshold requring at least 1 track to be truth matched
            matchedDVs.push_back(DVs[i]);
    }
    return matchedDVs;
}

// Important that the input indices are the MC indices, i.e the indices of the particles/tracks in the MC collection
ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> truthmatchDVs_vertex(ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> DVs,
                                                                                        const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& reco,
                                                                                        ROOT::VecOps::RVec<int> MCtruthind,
                                                                                        ROOT::VecOps::RVec<int> RPMCindex) {

    ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> matchedDVs;

    for (size_t i = 0; i < DVs.size(); ++i) {
        ROOT::VecOps::RVec<int> DVsRecoIndices;
        ROOT::VecOps::RVec<int> DVsMCIndices;
        DVsRecoIndices = FCCAnalyses::VertexingUtils::get_VertexRecoParticlesInd(DVs[i], reco);
        DVsMCIndices = get_VertexTrack2MC(DVsRecoIndices, RPMCindex);

        double n_matchedDV_tracks = 0.;
            for (size_t i = 0; i < DVsMCIndices.size(); ++i) {
                for (size_t j = 0; j < MCtruthind.size(); ++j) {
                    if (DVsMCIndices[i] == MCtruthind[j])
                        n_matchedDV_tracks += 1.;
                }    
            }

        // change threshold here. Either put an absolute one where the threshold is set requiring a minimum number of truth matched tracks
        // or have a "purity" threshold, i.e the fraction of truthmatched tracks of all tracks in the DV should be more than e.g 0.5
        double n_total_tracks_DV = FCCAnalyses::VertexingUtils::get_VertexNtrk(DVs[i]);
        double purity = n_matchedDV_tracks/n_total_tracks_DV;
        if (purity > 0.8)   //absolute threshold requring at least 1 track to be truth matched
            matchedDVs.push_back(DVs[i]);
    }
    return matchedDVs;
}

double get_min_distance(ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> DVs, edm4hep::Vector3d MC_vertices) {

    ROOT::VecOps::RVec<double> distance = FCCAnalyses::VertexingUtils::get_d3d_SV_obj(DVs, MC_vertices);
    double distance_min = ROOT::VecOps::Min(distance);

    return distance_min;
}

//ROOT::VecOps::RVec<float> 

std::vector<double> get_vector(double in1, double in2) {
    std::vector<double> result;
    result.push_back(in1);
    result.push_back(in2);
    return result;
}

std::size_t get_min_distance_index(ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> DVs, edm4hep::Vector3d MC_vertices) {

    ROOT::VecOps::RVec<double> distance = FCCAnalyses::VertexingUtils::get_d3d_SV_obj(DVs, MC_vertices);
    std::size_t distance_min_index = ROOT::VecOps::ArgMin(distance);

    return distance_min_index;
}


ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> truthmatchDVs_distance_event(ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> DVs, edm4hep::Vector3d scalar1_vertex, edm4hep::Vector3d scalar2_vertex) {
    ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> matchedDVs;
    int n_DVs = DVs.size();

    if (n_DVs == 0.) {
        return matchedDVs;
    }

    double min_distance_s1 = get_min_distance(DVs, scalar1_vertex);
    double min_distance_s2 = get_min_distance(DVs, scalar2_vertex);

    if (min_distance_s1 <= 2) {
        std::size_t index1 = get_min_distance_index(DVs, scalar1_vertex);
        matchedDVs.push_back(DVs.at(index1));
    }

    if (min_distance_s2 <= 2) {
        std::size_t index2 = get_min_distance_index(DVs, scalar2_vertex);
        matchedDVs.push_back(DVs.at(index2));
    }

    return matchedDVs;
}

ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> truthmatchDVs_distance_per_truthdecay(ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> DVs, edm4hep::Vector3d scalar_vertex) {
    ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> matchedDVs;
    int n_DVs = DVs.size();

    if (n_DVs == 0.) {
        return matchedDVs;
    }

    double min_distance = get_min_distance(DVs, scalar_vertex);

    if (min_distance <= 2) {
        std::size_t index = get_min_distance_index(DVs, scalar_vertex);
        matchedDVs.push_back(DVs.at(index));
    }

    return matchedDVs;
}

float get_reconstruction_efficency_per_truthdecay(ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> truthmatchedDVs, ROOT::VecOps::RVec<float> decaylength_scalar) {
    float n_truthmatchedDVs = truthmatchedDVs.size();

    if (n_truthmatchedDVs == 0.) {
        return 0.;
    }

    if (decaylength_scalar.at(0) > 2000) {
        return 0.;
    }
    return n_truthmatchedDVs;
}

ROOT::VecOps::RVec<float> get_reconstruction_efficency_radius(ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> truthmatchedDVs,
                                                                ROOT::VecOps::RVec<float> decaylengths_scalar,
                                                                FCCAnalyses::VertexingUtils::FCCAnalysesVertex PV) {
    float n_truthmatchedDVs = truthmatchedDVs.size();
    ROOT::VecOps::RVec<double> distances_DVs = FCCAnalyses::VertexingUtils::get_d3d_SV(truthmatchedDVs, PV);
    ROOT::VecOps::RVec<float> result;
    int volume_area = 2000/100; //radius of the ID = 2000 mm, binning = 100
    result.resize(volume_area);

    for (int i = 0; i < result.size(); ++i) {
        float j = i*volume_area;
        float k = (i+1)*volume_area;
        for (size_t m = 0; m < decaylengths_scalar.size(); ++m) {
            if (decaylengths_scalar.at(m) >= j || decaylengths_scalar.at(m) < k) {
                for (size_t n = 0; n < n_truthmatchedDVs; ++n) {
                    if (distances_DVs.at(n) >= j || distances_DVs.at(n) < k) result[i] = 1.;
                    //else result[i] = 0.;
                }
            }
        }
    }
    return result;
}

float get_reconstruction_efficency(ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> truthmatchedDVs, ROOT::VecOps::RVec<float> decaylengths_scalar) {
    float n_truthmatchedDVs = truthmatchedDVs.size();

    if (n_truthmatchedDVs == 0.) {
        return 0.;
    }

    float truthdecays = 0.;
    for (size_t i = 0; i < decaylengths_scalar.size(); ++i) {
        if (decaylengths_scalar.at(i) < 2000) {
            truthdecays += 1.;
        }
    }
    if (truthdecays == 0.) {
        return 0.;
    }
    float result = n_truthmatchedDVs/truthdecays;
    return result;
}

float get_fake_rate(ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> truthmatchedDVs, ROOT::VecOps::RVec<FCCAnalyses::VertexingUtils::FCCAnalysesVertex> DVs) {
    float n_DVs = DVs.size();

    if (n_DVs == 0.) {
        return 0.;
    }
    
    float nonmatchedDVs = n_DVs - truthmatchedDVs.size();
    float fake_rate = nonmatchedDVs/n_DVs;
    return fake_rate;
}

ROOT::VecOps::RVec<int> erase_first(ROOT::VecOps::RVec<int> vector) {
    ROOT::VecOps::RVec<int> result;
    for (int i = 1; i < vector.size(); ++i)
        result.push_back(vector.at(i));
    return result;
}

ROOT::VecOps::RVec<edm4hep::MCParticleData> get_MCparticles(ROOT::VecOps::RVec<int> list_of_indices,  ROOT::VecOps::RVec<edm4hep::MCParticleData> in) {
    ROOT::VecOps::RVec<edm4hep::MCParticleData> result;
    for(size_t i = 0; i < list_of_indices.size(); ++i) {
        result.push_back(FCCAnalyses::MCParticle::sel_byIndex(list_of_indices[i], in));
    }
    return result;
}

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
#END USER DEFINED CODE

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

            # .Define("particle_indices", "Particle1")
            # .Define("MCRecoAssoc0", "MCRecoAssociations0")
            # .Define("MCRecoAssoc1", "MCRecoAssociations1")
            # .Define("daughter_begin_indices", "Particle.daughters_begin")
            # .Define("daughter_end_indices", "Particle.daughters_end")

            # .Define("tracks_begin", "ReconstructedParticles.tracks_begin")
            # .Define("tracks_end", "ReconstructedParticles.tracks_end")

            # .Define("MC_FS", 'MCParticle::sel_genStatus(1)(Particle)')
            # .Define('MC_FS_pdg', 'MCParticle::get_pdg(MC_FS)')

            # MC event primary vertex
            #.Define("MC_PrimaryVertex",  "FCCAnalyses::MCParticle::get_EventPrimaryVertex(21)( Particle )" )

            # number of tracks
            #.Define("n_tracks","ReconstructedParticle2Track::getTK_n(EFlowTrack_1)")

            # Vertex fitting

            # Select the tracks that are reconstructed  as primaries
            .Define("RecoedPrimaryTracks",  "VertexFitterSimple::get_PrimaryTracks(EFlowTrack_1, true, 4.5, 20e-3, 300, 0., 0., 0.)")
            .Define("n_RecoedPrimaryTracks",  "ReconstructedParticle2Track::getTK_n( RecoedPrimaryTracks )")

            # the primary vertex
            # Input parameters are 1 = primary vertex, EFlowTrack_1 contains all tracks, bool beamspotconstraint = true, bsc sigma x/y/z
            # ATM RETURNING "ERROR": VertexFit::RegInv: null determinant for N = 2, 
            .Define("PrimaryVertexObject",   "VertexFitterSimple::VertexFitter_Tk (1, RecoedPrimaryTracks, true, 4.5, 20e-3, 300) ")
            .Define("PrimaryVertex",   "VertexingUtils::get_VertexData(PrimaryVertexObject)")



            # # Displaced vertex reconstruction
            
            # # select tracks with pT > 1 GeV
            # .Define('sel_tracks_pt', 'VertexingUtils::sel_pt_tracks(1)(EFlowTrack_1)')
            # # select tracks with |d0 |> 2 mm
            # .Define('sel_tracks', 'VertexingUtils::sel_d0_tracks(2)(sel_tracks_pt)')
            # # find the DVs
            # # ATM RETURNING: CheckDefPos: found <= 0 eigenvalue E(4) = -0.309843. CheckDefPos: input matrix NOT posite definite. Printing normalized matrix.
            # .Define("DV_evt_seltracks", "VertexFinderLCFIPlus::get_SV_event(sel_tracks, EFlowTrack_1, PrimaryVertexObject, true, 9., 40., 5.)")
            # # number of DVs
            # .Define('n_seltracks_DVs', 'VertexingUtils::get_n_SV(DV_evt_seltracks)')
            # # number of tracks from the DVs
            # .Define('n_trks_seltracks_DVs', 'VertexingUtils::get_VertexNtrk(DV_evt_seltracks)')
            # # invariant mass at the DVs (assuming the tracks to be pions)
            # .Define('invMass_seltracks_DVs', 'VertexingUtils::get_invM(DV_evt_seltracks)')

            # # get the chi2 distributions of the DVs from selected tracks
            # .Define("DV_evt_seltracks_chi2",    "VertexingUtils::get_chi2_SV(DV_evt_seltracks)")
            # .Define("DV_evt_seltracks_normchi2","VertexingUtils::get_norm_chi2_SV(DV_evt_seltracks)") # DV chi2 (normalised)

            # # get the decay radius of all the DVs from selected tracks
            # .Define("Reco_seltracks_DVs_Lxy","VertexingUtils::get_dxy_SV(DV_evt_seltracks, PrimaryVertexObject)")
            # .Define("Reco_seltracks_DVs_Lxyz","VertexingUtils::get_d3d_SV(DV_evt_seltracks, PrimaryVertexObject)")

            
            # # merge vertices with position within 10*error-of-position, get the tracks from the merged vertices and refit
            # .Define('merged_DVs', 'VertexingUtils::mergeVertices(DV_evt_seltracks, EFlowTrack_1)')
            # # number of merged DVs
            # .Define("merged_DVs_n", "VertexingUtils::get_n_SV(merged_DVs)")
            # # number of tracks from the merged DVs
            # .Define('n_trks_merged_DVs', 'VertexingUtils::get_VertexNtrk(merged_DVs)')
            # # invariant mass at the merged DVs
            # .Define('invMass_merged_DVs', 'VertexingUtils::get_invM(merged_DVs)')

            # # get the chi2 distributions of the merged DVs
            # .Define("merged_DVs_chi2",    "VertexingUtils::get_chi2_SV(merged_DVs)")
            # .Define("merged_DVs_normchi2","VertexingUtils::get_norm_chi2_SV(merged_DVs)") # DV chi2 (normalised)

            # # get the decay radius of all the merged DVs
            # .Define("Reco_DVs_merged_Lxy","VertexingUtils::get_dxy_SV(merged_DVs, PrimaryVertexObject)")
            # .Define("Reco_DVs_merged_Lxyz","VertexingUtils::get_d3d_SV(merged_DVs, PrimaryVertexObject)")

             #### Under construction

            # find the indices of the tracks belonging to each DV, atm the method only works for one FCCAnalyses vertex, 
            # should rewrite the method to work for a vector of vertices or write my own method at the top of this script
            # .Define('recoind_seltracks_DVs','VertexingUtils::get_VertexRecoParticlesInd(DV_evt_seltracks[0], ReconstructedParticles)')

            # returns a vector with the MC indices, like [index reco] = index MC
            # .Define('RP_MC_index', "ReconstructedParticle2MC::getRP2MC_index(MCRecoAssociations0, MCRecoAssociations1, ReconstructedParticles)")

            # # add code/function which returns the MC particles/indices to each vertex, i.e write a function which loops through the RP_MC_index with the relevant track indices
            # .Define('Vertex_MC_index', 'get_VertexTrack2MC(recoind_seltracks_DVs, RP_MC_index)')
            
            # get the indices of the Higgs and the scalars, sorted in order h hs hs
            # .Define('H2HSHS_indices', 'MCParticle::get_indices_ExclusiveDecay(25, {35, 35}, false, false)(Particle, Particle1)')

            # # get the indices of the b quarks from related scalar, 1 are from the "first" scalar and 2 from the second one
            # .Define('HS1_to_bb_indices', 'MCParticle::get_indices_MotherByIndex(H2HSHS_indices[1], {5, -5}, false, false, true, Particle, Particle1)')
            # .Define('HS2_to_bb_indices', 'MCParticle::get_indices_MotherByIndex(H2HSHS_indices[2], {5, -5}, false, false, true, Particle, Particle1)')

            # # remove the the scalar index from each list, keep only the b indices
            # .Define('bquarks1_indices', 'erase_first(HS1_to_bb_indices)')
            # .Define('bquarks2_indices', 'erase_first(HS2_to_bb_indices)')

            # # method to retrieve the decay chain from each scalar to stable particles
            # .Define('scalar1_decays_indices', 'get_list_of_FSparticles_indices(HS1_to_bb_indices[0], Particle, Particle1)')
            # .Define('scalar2_decays_indices', 'get_list_of_FSparticles_indices(HS2_to_bb_indices[0], Particle, Particle1)')

            # # method to retrieve the decay chain with the first branching of daughters from the b quarks of scalar 1
            # .Define('s1b1_decays_indices', 'get_list_of_FSparticles_indices(HS1_to_bb_indices[1], Particle, Particle1)')
            # .Define('s1b2_decays_indices', 'get_list_of_FSparticles_indices(HS1_to_bb_indices[2], Particle, Particle1)')

            # # method to retrieve the decay chain with the first branching of daughters from the b quarks of scalar 2
            # .Define('s2b1_decays_indices', 'get_list_of_FSparticles_indices(HS2_to_bb_indices[1], Particle, Particle1)')
            # .Define('s2b2_decays_indices', 'get_list_of_FSparticles_indices(HS2_to_bb_indices[2], Particle, Particle1)')


            # # get the b quarks, right now only picking the first daughter (hopefully b quark) from the list of decay particles from hs, should update this section maybe
            # .Define('bquarks1', 'myUtils::selMC_leg(0) (bquarks1_indices, Particle)')
            # .Define('bquarks2', 'myUtils::selMC_leg(0) (bquarks2_indices, Particle)')

            # # get the production vertex of the b's, to retrieve info like position and radius
            # # returns list of vertices, should check that the list only contains one vertex...
            # .Define('bquarks1_vertices', 'MCParticle::get_vertex(bquarks1)')
            # .Define('bquarks2_vertices', 'MCParticle::get_vertex(bquarks2)')

            # # get the distance between the DVs reconstructed with track selection and the MC b production vertices (for both scalar 1 and 2)
            # .Define('distance_b1_seltracks', 'VertexingUtils::get_d3d_SV_obj(DV_evt_seltracks, bquarks1_vertices.at(0))')
            # .Define('distance_b2_seltracks', 'VertexingUtils::get_d3d_SV_obj(DV_evt_seltracks, bquarks2_vertices.at(0))')

            # # get the distance between the DVs reconstructed including merging and the MC b production vertices (for both scalar 1 and 2)
            # .Define('distance_b1_merged', 'VertexingUtils::get_d3d_SV_obj(merged_DVs, bquarks1_vertices.at(0))')
            # .Define('distance_b2_merged', 'VertexingUtils::get_d3d_SV_obj(merged_DVs, bquarks2_vertices.at(0))')

            # # methods to see the minimal distance between MC vertex and DVs, i.e the closest DV
            # .Define('min_distance_b1_seltracks', 'get_min_distance(DV_evt_seltracks, bquarks1_vertices.at(0))')
            # .Define('min_distance_b2_seltracks', 'get_min_distance(DV_evt_seltracks, bquarks2_vertices.at(0))')
            # .Define('min_distance_seltracks', 'get_vector(min_distance_b1_seltracks, min_distance_b2_seltracks)')

            # .Define('min_distance_b1_merged', 'get_min_distance(merged_DVs, bquarks1_vertices.at(0))')
            # .Define('min_distance_b2_merged', 'get_min_distance(merged_DVs, bquarks2_vertices.at(0))')
            # .Define('min_distance_merged', 'get_vector(min_distance_b1_merged, min_distance_b2_merged)')


            # # select generated scalar HS
            # .Define('GenHS_PID',  'MCParticle::sel_pdgID(35, false)(Particle)')

            # # get the production vertex for the 2 HS in x y z
            # .Define('HS_vertex_x', 'MCParticle::get_vertex_x(GenHS_PID)')
            # .Define('HS_vertex_y', 'MCParticle::get_vertex_y(GenHS_PID)')
            # .Define('HS_vertex_z', 'MCParticle::get_vertex_z(GenHS_PID)')

            # # # get the vertex position for the first group of b quarks
            # .Define('b1_vertex_x', 'MCParticle::get_vertex_x(bquarks1)')
            # .Define('b1_vertex_y', 'MCParticle::get_vertex_y(bquarks1)')
            # .Define('b1_vertex_z', 'MCParticle::get_vertex_z(bquarks1)')

            # # # get the vertex position for the second group of b quarks
            # .Define('b2_vertex_x', 'MCParticle::get_vertex_x(bquarks2)')
            # .Define('b2_vertex_y', 'MCParticle::get_vertex_y(bquarks2)')
            # .Define('b2_vertex_z', 'MCParticle::get_vertex_z(bquarks2)')

            # # get the decay length of HS 1
            # .Define('decayLengthHS1', 'return sqrt((b1_vertex_x - HS_vertex_x.at(0))*(b1_vertex_x - HS_vertex_x.at(0)) + (b1_vertex_y - HS_vertex_y.at(0))*(b1_vertex_y - HS_vertex_y.at(0)) + (b1_vertex_z - HS_vertex_z.at(0))*(b1_vertex_z - HS_vertex_z.at(0)))')

            # # get the decay length of HS 2
            # .Define('decayLengthHS2', 'return sqrt((b2_vertex_x - HS_vertex_x.at(1))*(b2_vertex_x - HS_vertex_x.at(1)) + (b2_vertex_y - HS_vertex_y.at(1))*(b2_vertex_y - HS_vertex_y.at(1)) + (b2_vertex_z - HS_vertex_z.at(1))*(b2_vertex_z - HS_vertex_z.at(1)))')

            # # get decay length of both scalars in a vector for each event
            # .Define('decayLengthsHS', 'myUtils::get_both_scalars(decayLengthHS1, decayLengthHS2)')

            # ## Truth matching with requiring minimal distance method
            # ## get reconstriction efficiency and fake rate
            # .Define('truthmatchedDVs_mindist_seltracks', 'truthmatchDVs_distance_event(DV_evt_seltracks, bquarks1_vertices.at(0), bquarks2_vertices.at(0))')
            # .Define('distance_truthmatchedDVs_mindist_seltracks', 'VertexingUtils::get_d3d_SV(truthmatchedDVs_mindist_seltracks, PrimaryVertexObject)')
            # .Define('rec_eff_mindist_seltracks', 'get_reconstruction_efficency(truthmatchedDVs_mindist_seltracks, decayLengthsHS)')
            # .Define('fake_rate_mindist_seltracks', 'get_fake_rate(truthmatchedDVs_mindist_seltracks, DV_evt_seltracks)')
            
            # .Define('truthmatchedDVs_mindist_merged', 'truthmatchDVs_distance_event(merged_DVs, bquarks1_vertices.at(0), bquarks2_vertices.at(0))')
            # .Define('distance_truthmatchedDVs_mindist_merged', 'VertexingUtils::get_d3d_SV(truthmatchedDVs_mindist_merged, PrimaryVertexObject)')
            # .Define('rec_eff_mindist_merged', 'get_reconstruction_efficency(truthmatchedDVs_mindist_merged, decayLengthsHS)')
            # .Define('fake_rate_mindist_merged', 'get_fake_rate(truthmatchedDVs_mindist_merged, merged_DVs)')

            # # Acceptance true particles, decays with decaylength < 2000 mm
            # .Define('acceptance_true_scalars', 'get_acceptance_trueparticles(decayLengthsHS)')

            # truthmatch to the scalars seperatly to get the distance between the truthmatched DVs and the true decays
            #.Define('trutmatchedDVs_seltracks_s1', 'truthmatchDVs_vertex(DV_evt_seltracks, ReconstructedParticles, scalar1_decays_indices, RP_MC_index)')
            #.Define('distance_s1_truthmatchedDVs_seltracks_stable', 'VertexingUtils::get_d3d_SV_obj(trutmatchedDVs_seltracks_s1, bquarks1_vertices.at(0))')
            # .Define('truthmatch_score_seltracks_s1', 'truthmatchDV_score_all(scalar1_decays_indices, scalar2_decays_indices, DV_evt_seltracks, ReconstructedParticles, RP_MC_index)')
            # .Define('truthmatch_purity_seltracks_s1', 'truthmatchDV_purity_all(truthmatch_score_seltracks_s1, DV_evt_seltracks)')
            # .Define('rec_eff_seltracks_s1', 'get_reconstruction_efficency(truthmatchedDVs_seltracks_s1, decayLengthsHS)')
            # .Define('fake_rate_seltracks_s1', 'get_fake_rate(truthmatchedDVs_seltracks_s1, DV_evt_seltracks)')

            #.Define('trutmatchedDVs_seltracks_s2', 'truthmatchDVs_vertex(DV_evt_seltracks, ReconstructedParticles, scalar2_decays_indices, RP_MC_index)')
            #.Define('distance_s2_truthmatchedDVs_seltracks_stable', 'VertexingUtils::get_d3d_SV_obj(trutmatchedDVs_seltracks_s2, bquarks2_vertices.at(0))')
            # .Define('truthmatch_score_seltracks_s2', 'truthmatchDV_score_all(scalar1_decays_indices, scalar2_decays_indices, DV_evt_seltracks, ReconstructedParticles, RP_MC_index)')
            # .Define('truthmatch_purity_seltracks_s2', 'truthmatchDV_purity_all(truthmatch_score_seltracks_s2, DV_evt_seltracks)')
            # .Define('rec_eff_seltracks_s2', 'get_reconstruction_efficency(truthmatchedDVs_seltracks_s2, decayLengthsHS)')
            # .Define('fake_rate_seltracks_s2', 'get_fake_rate(truthmatchedDVs_seltracks_s2, DV_evt_seltracks)')

            # ## Truth matching with comparing tracks to MC method, comparing to true particles from the full decay chain to stable particles
            # ## get reconstriction efficiency and fake rate
            # .Define('truthmatchedDVs_seltracks_stable', 'truthmatchDVs_event(DV_evt_seltracks, ReconstructedParticles, scalar1_decays_indices_stable, scalar2_decays_indices_stable, RP_MC_index)')
            # .Define('truthmatch_score_seltracks_stable', 'truthmatchDV_score_all(scalar1_decays_indices_stable, scalar2_decays_indices_stable, DV_evt_seltracks, ReconstructedParticles, RP_MC_index)')
            # .Define('truthmatch_purity_seltracks_stable', 'truthmatchDV_purity_all(truthmatch_score_seltracks_stable, DV_evt_seltracks)')
            # .Define('rec_eff_seltracks_stable', 'get_reconstruction_efficency(truthmatchedDVs_seltracks_stable, decayLengthsHS)')
            # .Define('fake_rate_seltracks_stable', 'get_fake_rate(truthmatchedDVs_seltracks_stable, DV_evt_seltracks)')
            
            # .Define('truthmatchedDVs_merged_stable', 'truthmatchDVs_event(merged_DVs, ReconstructedParticles, scalar1_decays_indices_stable, scalar2_decays_indices_stable, RP_MC_index)')
            # .Define('truthmatch_score_merged_stable', 'truthmatchDV_score_all(scalar1_decays_indices_stable, scalar2_decays_indices_stable, merged_DVs, ReconstructedParticles, RP_MC_index)')
            # .Define('truthmatch_purity_merged_stable', 'truthmatchDV_purity_all(truthmatch_score_merged_stable, merged_DVs)')
            # .Define('rec_eff_merged_stable', 'get_reconstruction_efficency(truthmatchedDVs_merged_stable, decayLengthsHS)')
            # .Define('fake_rate_merged_stable', 'get_fake_rate(truthmatchedDVs_merged_stable, merged_DVs)')

            # ## Truth matching with comparing tracks to MC method, comparing to the first branching of decay particles from the scalar decays
            # ## get reconstriction efficiency and fake rate
            # .Define('truthmatchedDVs_seltracks_s', 'truthmatchDVs_event(DV_evt_seltracks, ReconstructedParticles, scalar1_decays_indices, scalar2_decays_indices, RP_MC_index)')
            # .Define('truthmatch_score_seltracks_s', 'truthmatchDV_score_all(scalar1_decays_indices, scalar2_decays_indices, DV_evt_seltracks, ReconstructedParticles, RP_MC_index)')
            # .Define('truthmatch_purity_seltracks_s', 'truthmatchDV_purity_all(truthmatch_score_seltracks_s, DV_evt_seltracks)')
            # .Define('rec_eff_seltracks_s', 'get_reconstruction_efficency(truthmatchedDVs_seltracks_s, decayLengthsHS)')
            # .Define('fake_rate_seltracks_s', 'get_fake_rate(truthmatchedDVs_seltracks_s, DV_evt_seltracks)')
            
            # .Define('truthmatchedDVs_merged_s', 'truthmatchDVs_event(merged_DVs, ReconstructedParticles, scalar1_decays_indices, scalar2_decays_indices, RP_MC_index)')
            # .Define('truthmatch_score_merged_s', 'truthmatchDV_score_all(scalar1_decays_indices, scalar2_decays_indices, merged_DVs, ReconstructedParticles, RP_MC_index)')
            # .Define('truthmatch_purity_merged_s', 'truthmatchDV_purity_all(truthmatch_score_merged_s, merged_DVs)')
            # .Define('rec_eff_merged_s', 'get_reconstruction_efficency(truthmatchedDVs_merged_s, decayLengthsHS)')
            # .Define('fake_rate_merged_s', 'get_fake_rate(truthmatchedDVs_merged_s, merged_DVs)')

            # ## Truth matching with comparing tracks to MC method, comparing to the first branching of decay particles from the b's from scalar 1 decay
            # ## get reconstriction efficiency and fake rate
            # .Define('truthmatchedDVs_seltracks_s1b', 'truthmatchDVs_event(DV_evt_seltracks, ReconstructedParticles, s1b1_decays_indices, s1b2_decays_indices, RP_MC_index)')
            # .Define('truthmatch_score_seltracks_s1b', 'truthmatchDV_score_all(s1b1_decays_indices, s1b2_decays_indices, DV_evt_seltracks, ReconstructedParticles, RP_MC_index)')
            # .Define('truthmatch_purity_seltracks_s1b', 'truthmatchDV_purity_all(truthmatch_score_seltracks_s1b, DV_evt_seltracks)')
            # .Define('rec_eff_seltracks_s1b', 'get_reconstruction_efficency(truthmatchedDVs_seltracks_s1b, decayLengthsHS)')
            # .Define('fake_rate_seltracks_s1b', 'get_fake_rate(truthmatchedDVs_seltracks_s1b, DV_evt_seltracks)')
            
            # .Define('truthmatchedDVs_merged_s1b', 'truthmatchDVs_event(merged_DVs, ReconstructedParticles, s1b1_decays_indices, s1b2_decays_indices, RP_MC_index)')
            # .Define('truthmatch_score_merged_s1b', 'truthmatchDV_score_all(s1b1_decays_indices, s1b2_decays_indices, merged_DVs, ReconstructedParticles, RP_MC_index)')
            # .Define('truthmatch_purity_merged_s1b', 'truthmatchDV_purity_all(truthmatch_score_merged_s1b, merged_DVs)')
            # .Define('rec_eff_merged_s1b', 'get_reconstruction_efficency(truthmatchedDVs_merged_s1b, decayLengthsHS)')
            # .Define('fake_rate_merged_s1b', 'get_fake_rate(truthmatchedDVs_merged_s1b, merged_DVs)')

            # ## Truth matching with comparing tracks to MC method, comparing to the first branching of decay particles from the b's from scalar 2 decay
            # ## get reconstriction efficiency and fake rate
            # .Define('truthmatchedDVs_seltracks_s2b', 'truthmatchDVs_event(DV_evt_seltracks, ReconstructedParticles, s2b1_decays_indices, s2b2_decays_indices, RP_MC_index)')
            # .Define('truthmatch_score_seltracks_s2b', 'truthmatchDV_score_all(s2b1_decays_indices, s2b2_decays_indices, DV_evt_seltracks, ReconstructedParticles, RP_MC_index)')
            # .Define('truthmatch_purity_seltracks_s2b', 'truthmatchDV_purity_all(truthmatch_score_seltracks_s2b, DV_evt_seltracks)')
            # .Define('rec_eff_seltracks_s2b', 'get_reconstruction_efficency(truthmatchedDVs_seltracks_s2b, decayLengthsHS)')
            # .Define('fake_rate_seltracks_s2b', 'get_fake_rate(truthmatchedDVs_seltracks_s2b, DV_evt_seltracks)')
            
            # .Define('truthmatchedDVs_merged_s2b', 'truthmatchDVs_event(merged_DVs, ReconstructedParticles, s2b1_decays_indices, s2b2_decays_indices, RP_MC_index)')
            # .Define('truthmatch_score_merged_s2b', 'truthmatchDV_score_all(s2b1_decays_indices, s2b2_decays_indices, merged_DVs, ReconstructedParticles, RP_MC_index)')
            # .Define('truthmatch_purity_merged_s2b', 'truthmatchDV_purity_all(truthmatch_score_merged_s2b, merged_DVs)')
            # .Define('rec_eff_merged_s2b', 'get_reconstruction_efficency(truthmatchedDVs_merged_s2b, decayLengthsHS)')
            # .Define('fake_rate_merged_s2b', 'get_fake_rate(truthmatchedDVs_merged_s2b, merged_DVs)')

            # # Number of DVs with distance and invariant mass selection applied
            # .Define("filter_n_DVs_seltracks", "filter_n_DVs(Reco_seltracks_DVs_Lxyz, invMass_seltracks_DVs)")
            # .Define("filter_n_DVs_merged", "filter_n_DVs(Reco_DVs_merged_Lxyz, invMass_merged_DVs)")
            
            # ## Reconstructed electrons and muons

            # # Electrons
            # .Alias('Electron0', 'Electron#0.index')
            # .Define('RecoElectrons',  'ReconstructedParticle::get(Electron0, ReconstructedParticles)') 
            # .Define('n_RecoElectrons',  'ReconstructedParticle::get_n(RecoElectrons)') #count how many electrons are in the event in total

            # # some kinematics of the reconstructed electrons and positrons
            # .Define("RecoElectron_e", "ReconstructedParticle::get_e(RecoElectrons)")
            # .Define("RecoElectron_p", "ReconstructedParticle::get_p(RecoElectrons)")
            # .Define("RecoElectron_pt", "ReconstructedParticle::get_pt(RecoElectrons)")
            # .Define("RecoElectron_px", "ReconstructedParticle::get_px(RecoElectrons)")
            # .Define("RecoElectron_py", "ReconstructedParticle::get_py(RecoElectrons)")
            # .Define("RecoElectron_pz", "ReconstructedParticle::get_pz(RecoElectrons)")
            # .Define("RecoElectron_charge",  "ReconstructedParticle::get_charge(RecoElectrons)")

            # # finding the invariant mass of the reconstructed electron and positron pair
            # .Define("Reco_ee_energy", "if ((n_RecoElectrons>1) && (RecoElectron_charge.at(0) != RecoElectron_charge.at(1))) return (RecoElectron_e.at(0) + RecoElectron_e.at(1)); else return float(-1.);")
            # .Define("Reco_ee_px", "if ((n_RecoElectrons>1) && (RecoElectron_charge.at(0) != RecoElectron_charge.at(1))) return (RecoElectron_px.at(0) + RecoElectron_px.at(1)); else return float(-1.);")
            # .Define("Reco_ee_py", "if ((n_RecoElectrons>1) && (RecoElectron_charge.at(0) != RecoElectron_charge.at(1))) return (RecoElectron_py.at(0) + RecoElectron_py.at(1)); else return float(-1.);")
            # .Define("Reco_ee_pz", "if ((n_RecoElectrons>1) && (RecoElectron_charge.at(0) != RecoElectron_charge.at(1))) return (RecoElectron_pz.at(0) + RecoElectron_pz.at(1)); else return float(-1.);")
            # .Define("Reco_ee_invMass", "if ((n_RecoElectrons>1) && (RecoElectron_charge.at(0) != RecoElectron_charge.at(1))) return sqrt(Reco_ee_energy*Reco_ee_energy - Reco_ee_px*Reco_ee_px - Reco_ee_py*Reco_ee_py - Reco_ee_pz*Reco_ee_pz ); else return float(-1.);")


            # # Muons
            # .Alias('Muon0', 'Muon#0.index')
            # .Define('RecoMuons',  'ReconstructedParticle::get(Muon0, ReconstructedParticles)')
            # .Define('n_RecoMuons',  'ReconstructedParticle::get_n(RecoMuons)') #count how many muons are in the event in total

            # # some kinematics of the reconstructed muons
            # .Define("RecoMuon_e",      "ReconstructedParticle::get_e(RecoMuons)")
            # .Define("RecoMuon_p",      "ReconstructedParticle::get_p(RecoMuons)")
            # .Define("RecoMuon_pt",      "ReconstructedParticle::get_pt(RecoMuons)")
            # .Define("RecoMuon_px",      "ReconstructedParticle::get_px(RecoMuons)")
            # .Define("RecoMuon_py",      "ReconstructedParticle::get_py(RecoMuons)")
            # .Define("RecoMuon_pz",      "ReconstructedParticle::get_pz(RecoMuons)")
            # .Define("RecoMuon_charge",  "ReconstructedParticle::get_charge(RecoMuons)")

            # # finding the invariant mass of the reconstructed muon pair
            # .Define("Reco_mumu_energy", "if ((n_RecoMuons>1) && (RecoMuon_charge.at(0) != RecoMuon_charge.at(1))) return (RecoMuon_e.at(0) + RecoMuon_e.at(1)); else return float(-1.);")
            # .Define("Reco_mumu_px", "if ((n_RecoMuons>1) && (RecoMuon_charge.at(0) != RecoMuon_charge.at(1))) return (RecoMuon_px.at(0) + RecoMuon_px.at(1)); else return float(-1.);")
            # .Define("Reco_mumu_py", "if ((n_RecoMuons>1) && (RecoMuon_charge.at(0) != RecoMuon_charge.at(1))) return (RecoMuon_py.at(0) + RecoMuon_py.at(1)); else return float(-1.);")
            # .Define("Reco_mumu_pz", "if ((n_RecoMuons>1) && (RecoMuon_charge.at(0) != RecoMuon_charge.at(1))) return (RecoMuon_pz.at(0) + RecoMuon_pz.at(1)); else return float(-1.);")
            # .Define("Reco_mumu_invMass", "if ((n_RecoMuons>1) && (RecoMuon_charge.at(0) != RecoMuon_charge.at(1))) return sqrt(Reco_mumu_energy*Reco_mumu_energy - Reco_mumu_px*Reco_mumu_px - Reco_mumu_py*Reco_mumu_py - Reco_mumu_pz*Reco_mumu_pz ); else return float(-1.);")

        )
        return df2


    #__________________________________________________________
    #Mandatory: output function, please make sure you return the branchlist as a python list
    def output():
        branchList = [
            # "particle_indices",
            # "MCRecoAssoc0",
            # "MCRecoAssoc1",
            # "daughter_begin_indices",
            # "daughter_end_indices",
            # "tracks_begin",
            # "tracks_end",
            
            # "n_tracks",
            # "n_RecoedPrimaryTracks",

            "PrimaryVertex",

            ## DVs reconstructed from selected tracks
            # 'n_seltracks_DVs',
            # 'n_trks_seltracks_DVs',
            # 'invMass_seltracks_DVs',
            # "DV_evt_seltracks_chi2",
            # "DV_evt_seltracks_normchi2",
            # "Reco_seltracks_DVs_Lxy",
            # "Reco_seltracks_DVs_Lxyz",

            ## DVs reconstructed including merging
            # "merged_DVs_n",
            # 'n_trks_merged_DVs',
            # 'invMass_merged_DVs',
            # "merged_DVs_chi2",
            # "merged_DVs_normchi2",
            # "Reco_DVs_merged_Lxy",
            # "Reco_DVs_merged_Lxyz",

            ## Some truthmatching indices checking
            # 'recoind_seltracks_DVs',
            # 'RP_MC_index',
            # 'Vertex_MC_index',

            # ## MC indices of true particles
            # 'H2HSHS_indices',
            # 'HS1_to_bb_indices',
            # 'HS2_to_bb_indices',
            # 'bquarks1_indices',
            # 'bquarks2_indices',
            # 'scalar1_decays_indices',
            # 'scalar2_decays_indices',
            # 's1b1_decays_indices',
            # 's1b2_decays_indices',
            # 's2b1_decays_indices',
            # 's2b2_decays_indices',

            
            # ## Distances between true decay and DVs
            # 'distance_b1_seltracks',
            # 'distance_b2_seltracks',
            # 'distance_b1_merged',
            # 'distance_b2_merged',
            # 'min_distance_b1_seltracks',
            # 'min_distance_b2_seltracks',
            # 'min_distance_b1_merged',
            # 'min_distance_b2_merged',
            # 'min_distance_seltracks',
            # 'min_distance_merged',

            # ## Decay lengths of true particles
            # 'decayLengthHS1',
            # 'decayLengthHS2',
            # 'decayLengthsHS',

            # # Acceptance
            # 'acceptance_true_scalars',

            # ## Reconstruction efficiency and fake rate, etc
            # 'distance_truthmatchedDVs_mindist_seltracks',
            # 'rec_eff_mindist_seltracks',
            # 'fake_rate_mindist_seltracks',
            # 'distance_truthmatchedDVs_mindist_merged',
            # 'rec_eff_mindist_merged',
            # 'fake_rate_mindist_merged',

            # 'truthmatch_score_seltracks_stable',
            # 'truthmatch_purity_seltracks_stable',
            # 'rec_eff_seltracks_stable',
            # 'fake_rate_seltracks_stable',
            # # 'truthmatch_score_merged_stable',
            # # 'truthmatch_purity_merged_stable',
            # # 'rec_eff_merged_stable',
            # # 'fake_rate_merged_stable',


            #'truthmatch_score_seltracks_s1',
            #'truthmatch_purity_seltracks_s1',
            #'truthmatch_score_seltracks_s2',
            #'truthmatch_purity_seltracks_s2',

            # 'truthmatch_score_seltracks_s',
            # 'truthmatch_purity_seltracks_s',
            # 'rec_eff_seltracks_s',
            # 'fake_rate_seltracks_s',
            # # 'truthmatch_score_merged_s',
            # # 'truthmatch_purity_merged_s',
            # # 'rec_eff_merged_s',
            # # 'fake_rate_merged_s',

            # 'truthmatch_score_seltracks_s1b',
            # 'truthmatch_purity_seltracks_s1b',
            # 'rec_eff_seltracks_s1b',
            # 'fake_rate_seltracks_s1b',
            # # 'truthmatch_score_merged_s1b',
            # # 'truthmatch_purity_merged_s1b',
            # # 'rec_eff_merged_s1b',
            # # 'fake_rate_merged_s1b',

            # 'truthmatch_score_seltracks_s2b',
            # 'truthmatch_purity_seltracks_s2b',
            # 'rec_eff_seltracks_s2b',
            # 'fake_rate_seltracks_s2b',
            # # 'truthmatch_score_merged_s2b',
            # # 'truthmatch_purity_merged_s2b',
            # # 'rec_eff_merged_s2b',
            # # 'fake_rate_merged_s2b',

            ## Vertex selection
            # "filter_n_DVs_seltracks",
            # "filter_n_DVs_merged",

            ## Reconstructed electrons
            # 'n_RecoElectrons',
            # "RecoElectron_e",
            # "RecoElectron_p",
            # "RecoElectron_pt",
            # "RecoElectron_px",
            # "RecoElectron_py",
            # "RecoElectron_pz",
            # "RecoElectron_charge",
            # "Reco_ee_invMass",

            ## Reconstructed muons
            # 'n_RecoMuons',
            # "RecoMuon_e",
            # "RecoMuon_p",
            # "RecoMuon_pt",
            # "RecoMuon_px",
            # "RecoMuon_py",
            # "RecoMuon_pz",
            # "RecoMuon_charge",
            # "Reco_mumu_invMass",
        ]
        return branchList