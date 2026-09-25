# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# Seven energy points of the 4.600-4.700 GeV scan (range 4.5995 - 4.6988 GeV)
scan_samples = %w[703_4600 706_4610 706_4620 706_4640 706_4660 706_4680 706_4700]
data_points  = scan_samples.map { |s| DatasetManager.real_data.find(s) }
incmc_points = scan_samples.map { |s| DatasetManager.inclusive_mc.find(s) }

# Decay card for the signal mode: e+e- -> Lambda_c+ anti-Lambda_c-, Lambda_c+ -> p pi0, pi0 -> gamma gamma
decay_card_pi0 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay

  Decay Lambda_c+
  1.0000 p+ pi0 PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Decay card for the reference mode: Lambda_c+ -> p eta, eta -> gamma gamma
decay_card_eta = <<~DECAYCARD
  Decay psi(4260)
  1.0000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay

  Decay Lambda_c+
  1.0000 p+ eta PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# 50k-event exclusive MC for each of the two decay modes, generated at every energy point
exMC_pi0 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "Lc_p_pi0"
  config.events        = 50000
  config.decay_card    = decay_card_pi0
  config.cross_section = :default
end

exMC_eta = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "Lc_p_eta"
  config.events        = 50000
  config.decay_card    = decay_card_eta
  config.cross_section = :default
end

### Event selection (BOSS) ###
# Common selection chain shared by both modes (they differ only in the pi0/eta mass constraint)
common_selection = Selection.new
  .select_track {
    cos_theta 0.93        # |cos(theta)| < 0.93
    Vz        10.0        # |Vz| < 10 cm
    Vr        1.0         # Vr < 1 cm
    nChrp    ">=1"        # at least one positively charged track
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025   # EMC barrel: E > 25 MeV
    energyThreshold_e 0.050   # EMC endcap: E > 50 MeV
    angle_to_track    10.0    # angle to the nearest charged track > 10 deg
    nGam             ">=2"    # at least two photons
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]   # p+ / anti-p- separated from K and pi
    nprp ">=1"                                  # at least one proton
  }
  .remove([:prp <= :chrgp])                     # proton-vs-charged-pion momentum-ordering removal
  .assign({:chrgp => :pip, :chrgn => :pim})     # remaining tracks assigned as pi+ / pi-
  .select_isolated_photon {
    angle_to_prp_track 20.0   # isolated from the primary proton track
    nGam ">=2"
  }

# Signal mode: 1C fit of the photon pair to the pi0 mass, then 4C fit of (p, pi0)
sel_pi0 = common_selection.dup
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200
    npi0 ">=1"
  }
  .kinematic_fit([:prp, :pi0]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

# Reference mode: identical chain but 1C fit of the photon pair to the eta mass
sel_eta = common_selection.dup
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 200
    neta ">=1"
  }
  .kinematic_fit([:prp, :eta]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

# --- Signal-mode algorithm (Lambda_c+ -> p pi0) ---
alg_pi0 = Algorithm.new("LcToPPi0")
alg_pi0.set_header(["LcToPPi0Alg/LcToPPi0.h"])
       .set_constant({"ECMS" => [:double, 4.649]})   # representative scan energy (see :ecms_scan note)
       .note(:ecms_scan, "CMS energy varies across the seven scan points (4.5995-4.6988 GeV); the fixed ECMS constant must be re-set per energy point for each job")
       .note(:helix_correction, "helix-parameter correction applied to all charged tracks before the kinematic fit")
       .note(:proton_pion_removal, "proton candidates removed from the positive-track list by a momentum-ordering criterion against charged pions before the remaining tracks are assigned as pi+/pi-")
       .note(:background_veto, "Particle Transformer DNN (20-model ensemble, iterative weighting, trained jointly for p-eta and p-pi0) required to give score > 0.95, evaluated on all charged tracks and on isolated showers not associated with the Lambda_c+ candidate; not expressible in the BOSS DSL")
alg_pi0.with_decay_card(decay_card_pi0).apply(sel_pi0)
root_files_pi0 = alg_pi0.execute_on(data_points + incmc_points + exMC_pi0)

# --- Reference-mode algorithm (Lambda_c+ -> p eta) ---
alg_eta = Algorithm.new("LcToPEta")
alg_eta.set_header(["LcToPEtaAlg/LcToPEta.h"])
       .set_constant({"ECMS" => [:double, 4.649]})   # representative scan energy (see :ecms_scan note)
       .note(:ecms_scan, "CMS energy varies across the seven scan points (4.5995-4.6988 GeV); the fixed ECMS constant must be re-set per energy point for each job")
       .note(:helix_correction, "helix-parameter correction applied to all charged tracks before the kinematic fit")
       .note(:proton_pion_removal, "proton candidates removed from the positive-track list by a momentum-ordering criterion against charged pions before the remaining tracks are assigned as pi+/pi-")
       .note(:background_veto, "Particle Transformer DNN (20-model ensemble, iterative weighting, trained jointly for p-eta and p-pi0) required to give score > 0.95, evaluated on all charged tracks and on isolated showers not associated with the Lambda_c+ candidate; not expressible in the BOSS DSL")
alg_eta.with_decay_card(decay_card_eta).apply(sel_eta)
root_files_eta = alg_eta.execute_on(data_points + incmc_points + exMC_eta)