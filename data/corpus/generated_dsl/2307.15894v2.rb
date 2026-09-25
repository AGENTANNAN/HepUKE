# =============================================================================
# e+e- -> Sigma+ Sigma- cross-section scan
# sqrt(s) = 2.3960, 2.6454, 2.9000 GeV  (BESIII R-scan, dataset 713)
# =============================================================================

### ---------------------------------------------------------------------------
### Dataset description
### ---------------------------------------------------------------------------
# Real data at the three energy points.  The 2.6454 GeV point is formed by
# combining the 2.6444 (Rscan_2644) and 2.6464 (Rscan_2646) datasets.
data_2396 = DatasetManager.real_data.find("713_Rscan_2396")
data_2644 = DatasetManager.real_data.find("713_Rscan_2644")
data_2646 = DatasetManager.real_data.find("713_Rscan_2646")
data_2900 = DatasetManager.real_data.find("713_Rscan_2900")

# Corresponding inclusive MC samples
incMC_2396 = DatasetManager.inclusive_mc.find("713_Rscan_2396")
incMC_2644 = DatasetManager.inclusive_mc.find("713_Rscan_2644")
incMC_2646 = DatasetManager.inclusive_mc.find("713_Rscan_2646")
incMC_2900 = DatasetManager.inclusive_mc.find("713_Rscan_2900")

# ConExc decay card (continuum / R-scan Born cross section):
# vpho -> Sigma+ anti-Sigma-,  Sigma+ -> p+ pi0,  anti-Sigma- -> anti-p- pi0,  pi0 -> gamma gamma.
# NOTE: no explicit `Particle vpho` line - for a multi-energy scan the DSL injects
# `Particle vpho <ECMS> 0.0` per energy point automatically.
decay_card_sigma_conexc = <<~DECAYCARD
    Decay vpho
    1.0000  Sigma+  anti-Sigma-   ConExc;
    Enddecay

    Decay Sigma+
    1.0000  p+  pi0   PHSP;
    Enddecay

    Decay anti-Sigma-
    1.0000  anti-p-  pi0   PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma   PHSP;
    Enddecay

    End
DECAYCARD

# 6M-event exclusive ConExc signal MC at each energy point
# (same decay card, several distinct energy points -> create_exclusive_mc_for)
exMCs_signal = DatasetManager.create_exclusive_mc_for([data_2396, data_2644, data_2646, data_2900]) do |config|
  config.sample_name   = "exmc_sigma_sigmabar_conexc"
  config.events        = 6_000_000
  config.decay_card    = decay_card_sigma_conexc
  config.cross_section = :default
end
exMCs_signal.each { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }

### ---------------------------------------------------------------------------
### Single-tag selection at 2.3960 GeV
### (reconstruct only anti-Sigma- -> anti-p pi0)
### ---------------------------------------------------------------------------
alg_st = Algorithm.new("SigmaSigmaST2396")
alg_st.set_header(["SigmaSigmaST2396Alg/SigmaSigmaST2396.h"])
      .set_constant({"ECMS" => [:double, 2.396]})
      .set_alias({"std::vector<double>" => "Vdouble"})
      .note(:pid_correction_method, "at 2.396 GeV the proton PID uses dE/dx only (TOF/MUC information not used); protons identified against kaons and pions")
      .note(:delta_e_mbc_window, "anti-Sigma- tag selected with dE in [-0.013, 0.005] GeV and beam-constrained mass Mbc in [1.185, 1.191] GeV/c^2; both tag variables are stored and windowed in the ROOT analysis")

st_selection = Selection.new
    .select_track {          # charged tracks
      cos_theta 0.93         # |cos(theta)| < 0.93
      Vz        10.0         # |Vz| < 10 cm
      Vr        1.0          # Vr < 1 cm
      nChrn     ">=1"        # at least one negative track
    }
    .select_photon {         # photons
      tdc_emc_start   0
      tdc_emc_end     14
      angle_to_track  20.0   # > 20 degrees from any charged track
      energyThreshold_b 0.025  # 25 MeV (barrel)
      energyThreshold_e 0.050  # 50 MeV (endcap)
      nGam ">=2"             # at least two photons
    }
    .pid(method: :probability) {   # PID
      prob_cut 0.001
      identify :proton, against: [:kaon, :pion]   # p+ and anti-p- (dE/dx only at this energy)
      nprm ">=1"             # at least one antiproton
    }
    .select_isolated_photon {      # isolated photons
      angle_to_prm_track 20.0      # > 20 degrees from the antiproton track
      nGam ">=2"                   # at least two isolated photons
    }
    .kalman_kinematic_fit([:gamma, :gamma]) {   # gamma-gamma fit to the pi0 mass
      invariant_mass_of(:gamma, :gamma).within(0.126, 0.139)   # pi0 mass window
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 25
      npi0 ">=1"                 # at least one pi0
    }
    .kinematic_fit([:prm, :pi0]) {   # combine anti-p and pi0 (anti-Sigma- candidate)
      nominal
      chi2_cut 200
    }

alg_st.with_decay_card(decay_card_sigma_conexc).apply(st_selection)

### ---------------------------------------------------------------------------
### Double-tag selection at 2.6454 GeV and 2.9000 GeV
### (reconstruct Sigma+ -> p pi0 and anti-Sigma- -> anti-p pi0, one pi0 missing)
### ---------------------------------------------------------------------------
dt_selection = Selection.new
    .select_track {          # charged tracks
      cos_theta 0.93
      Vz        10.0
      Vr        1.0
      nChrp     ">=1"        # at least one positive track
      nChrn     ">=1"        # at least one negative track
    }
    .select_photon {         # same photon criteria as the single tag
      tdc_emc_start   0
      tdc_emc_end     14
      angle_to_track  20.0
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      nGam ">=2"
    }
    .pid(method: :probability) {   # same PID criteria, both charges required
      prob_cut 0.001
      identify :proton, against: [:kaon, :pion]
      nprp ">=1"             # at least one proton
      nprm ">=1"             # at least one antiproton
    }
    .select_isolated_photon {      # isolated photons (> 20 deg from proton OR antiproton track)
      angle_to_prp_track 20.0
      angle_to_prm_track 20.0
      nGam ">=2"
    }
    .kalman_kinematic_fit([:gamma, :gamma]) {   # gamma-gamma fit to the pi0 mass
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 25
      npi0 ">=1"
    }
    .kinematic_fit([:prp, :prm, :pi0]) {   # 2C fit to p anti-p pi0 with one pi0 missing
      nominal
      miss_track_of(:pi0)                  # second pi0 undetected
      constrain_four_momentum
      invariant_mass_of(:prp, :pi0).within(1.175, 1.200)   # M_Sigma tag window
      invariant_mass_of(:prm, :pi0).within(1.175, 1.200)   # M_Sigma tag window
      chi2_cut 15
    }

# 2.6454 GeV algorithm (runs over the combined 2.6444 + 2.6464 data)
alg_dt_2645 = Algorithm.new("SigmaSigmaDT2645")
alg_dt_2645.set_header(["SigmaSigmaDT2645Alg/SigmaSigmaDT2645.h"])
          .set_constant({"ECMS" => [:double, 2.6454]})
          .set_alias({"std::vector<double>" => "Vdouble"})
          .note(:tag_fit, "2C kinematic fit to p anti-p pi0 with one pi0 missing, chi2 < 15")
          .note(:best_candidate_selection, "among the Sigma candidates the one with minimum |M(p/anti-p pi0) - M_Sigma| is chosen; M_Sigma tag window [1.175, 1.200] GeV/c^2")
          .note(:recoil_mass_observable, "the recoil mass against the reconstructed p/anti-p pi0 system is the primary observable for the cross-section extraction")

# 2.9000 GeV algorithm
alg_dt_2900 = Algorithm.new("SigmaSigmaDT2900")
alg_dt_2900.set_header(["SigmaSigmaDT2900Alg/SigmaSigmaDT2900.h"])
          .set_constant({"ECMS" => [:double, 2.900]})
          .set_alias({"std::vector<double>" => "Vdouble"})
          .note(:tag_fit, "2C kinematic fit to p anti-p pi0 with one pi0 missing, chi2 < 15")
          .note(:best_candidate_selection, "among the Sigma candidates the one with minimum |M(p/anti-p pi0) - M_Sigma| is chosen; M_Sigma tag window [1.175, 1.200] GeV/c^2")
          .note(:recoil_mass_observable, "the recoil mass against the reconstructed p/anti-p pi0 system is the primary observable for the cross-section extraction")

alg_dt_2645.with_decay_card(decay_card_sigma_conexc).apply(dt_selection)
alg_dt_2900.with_decay_card(decay_card_sigma_conexc).apply(dt_selection)

### ---------------------------------------------------------------------------
### Execute on real data, inclusive MC and exclusive ConExc signal MC
### (exMCs_signal order follows [data_2396, data_2644, data_2646, data_2900])
### ---------------------------------------------------------------------------
root_files_st   = alg_st.execute_on([data_2396, incMC_2396, exMCs_signal[0]])

root_files_dt_2645 = alg_dt_2645.execute_on([
  data_2644, data_2646,
  incMC_2644, incMC_2646,
  exMCs_signal[1], exMCs_signal[2]
])

root_files_dt_2900 = alg_dt_2900.execute_on([data_2900, incMC_2900, exMCs_signal[3]])