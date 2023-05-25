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
        # 'p8_ee_ZH_ecm240':{},
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
outputDir = "Vertex_output_stage1/"

#Optional: ncpus, default is 4
nCPUS       = 8

#Optional running on HTCondor, default is False
runBatch    = False
#runBatch    = True

#Optional batch queue name when running on HTCondor, default is workday
#batchQueue = "longlunch"

#Optional computing account when running on HTCondor, default is group_u_FCC.local_gen
#compGroup = "group_u_FCC.local_gen"

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

            # number of tracks
            .Define("n_tracks","ReconstructedParticle2Track::getTK_n(EFlowTrack_1)")

            # Vertex fitting

            # Select the tracks that are reconstructed  as primaries
            .Define("RecoedPrimaryTracks",  "VertexFitterSimple::get_PrimaryTracks(EFlowTrack_1, true, 4.5, 20e-3, 300, 0., 0., 0.)")
            .Define("n_RecoedPrimaryTracks",  "ReconstructedParticle2Track::getTK_n( RecoedPrimaryTracks )")

            # Fit the primary vertex
            # Input parameters are 1 = primary vertex, EFlowTrack_1 contains all tracks, bool beamspotconstraint = true, bsc sigma x/y/z
            .Define("PrimaryVertexObject",   "VertexFitterSimple::VertexFitter_Tk (1, RecoedPrimaryTracks, true, 4.5, 20e-3, 300) ")
            .Define("PrimaryVertex",   "VertexingUtils::get_VertexData(PrimaryVertexObject)")

        )
        return df2


    #__________________________________________________________
    #Mandatory: output function, please make sure you return the branchlist as a python list
    def output():
        branchList = [
            "n_tracks",
            "n_RecoedPrimaryTracks",
            "PrimaryVertex",
        ]
        return branchList
 