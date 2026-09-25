### Dataset description ###
# J/psi data at ECM = 3.097 GeV (~1.0087e10 events) and matching inclusive MC
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Full signal decay card (EvtGen):
#   J/psi -> Xi- anti-Xi+ ;  Xi- -> Sigma+ e- e-  (LNV, |dS|=|dL|=2)
#   anti-Xi+ -> anti-Lambda0 pi+ ; anti-Lambda0 -> anti-p- pi+ ; Sigma+ -> p+ pi0 ; pi0 -> gamma gamma
decay_card_signal = <<~DECAYCARD
  Decay J/psi
  1.0000 Xi- anti-Xi+ PHSP;
  Enddecay

  Decay Xi-
  1.0000 Sigma+ e- e- PHSP;
  Enddecay

  Decay anti-Xi+
  1.0000 anti-Lambda0 pi+ PHSP;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+ HypWK;
  Enddecay

  Decay Sigma+
  1.0000 p+ pi0 HypWK;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# 1,000,000 exclusive MC events for the full signal decay card
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_Xim_SigmaPee_LNV"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "XimLNV"
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})

# Single tag: anti-Xi+ -> anti-Lambda pi+, anti-Lambda -> anti-p pi+ (reconstructed from scratch),
# with the signal side Xi- -> Sigma+ e- e-, Sigma+ -> p pi0, pi0 -> gamma gamma.
event_selection = Selection.new
event_selection
  .select_track {
    cos_theta 0.93     # |cos(theta)| < 0.93
    Vz        100.0    # |Vz| < 100 cm
    Vr        10.0     # Vr < 10 in the transverse plane
    nChrp     ">=3"    # at least three positive tracks
    nChrn     ">=2"    # at least two negative tracks
  }
  .select_photon {
    tdc_emc_start     0        # EMC TDC window 0-14
    tdc_emc_end       14
    angle_to_track    10.0     # > 10 deg from any charged track
    energyThreshold_b 0.025    # E > 25 MeV (barrel)
    energyThreshold_e 0.050    # E > 50 MeV (endcap)
    nGam              ">=2"    # at least two photons
  }
  .pid(method: :probability) {
    prob_cut 0.001
    # electrons cannot be separated by identify(); use the high-momentum lepton branch
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    identify :proton, against: [:kaon, :pion]   # p+ and anti-p vs K and pi
    nprp ">=1"     # at least one proton
    nprm ">=1"     # at least one anti-proton
    nlm  ">=1"     # at least one electron (negative lepton)
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn, :lp <= :chrgp, :lm <= :chrgn])  # drop identified p, anti-p, e from track lists
  .select_isolated_photon {
    angle_to_prp_track 20.0    # > 20 deg from the primary (proton) track
    nGam ">=2"
  }
  .assign({:chrgp => :pip, :chrgn => :pim})   # remaining tracks are pions
  # anti-Lambda -> anti-p pi+ (secondary vertex, mass difference minimised)
  .secondary_vertex_fit([:prm, :pip]) {
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # anti-Xi+ -> anti-Lambda pi+ (secondary vertex, mass difference minimised)
  .secondary_vertex_fit([:Lambda_bar, :pip]) {
    build_virtual_particle(:Xi_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # pi0 -> gamma gamma (Kalman fit, gamma-gamma mass constrained to nominal pi0)
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 30
    npi0 ">=1"
  }
  # Nominal 4C fit on anti-Xi+, p, pi0 and e-, with the second electron treated as missing
  .kinematic_fit([:Xi_bar, :prp, :pi0, :lm, :lm]) {
    nominal
    constrain_four_momentum
    miss_track_of :lm    # the second electron is not reconstructed
    chi2_cut 20
  }

my_Algorithm
  .note(:signal_model, "the LNV decay Xi- -> Sigma+ e- e- is modelled with the Barbero et al. LNV width (B=0); the decay card uses PHSP as a placeholder until the dedicated LNV EvtGen model is available")
  .with_decay_card(decay_card_signal).apply(event_selection)

root_files = my_Algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])