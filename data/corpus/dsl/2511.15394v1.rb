# Search for Xi- -> Sigma+ e- e- (LNV, |Delta S|=|Delta L|=2)
# J/psi data, ~1.0087e10 events, ECM = 3.097 GeV
# UL < 2.0e-5 at 90% CL

jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

decay_card_sig = <<~DECAYCARD
    Decay J/psi
    1.0000 Xi- anti-Xi+ PHSP;
    Enddecay

    Decay Xi-
    1.0000 Sigma+ e- e- PHSP;
    Enddecay

    Decay anti-Xi+
    1.0000 anti-Lambda0 pi+ PHSP;
    Enddecay

    Decay Sigma+
    1.0000 p+ pi0 PHSP;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_XiToSigmaEE_3097"
  config.related_dataset = jpsi_data
  config.events = 1_000_000
  config.decay_card = decay_card_sig
  config.cross_section = :default
end

alg = Algorithm.new("XiToSigmaEE")
alg.set_header(["XiToSigmaEEAlg/XiToSigmaEE.h"])
   .set_constant({ "ECMS" => [:double, 3.097] })
   .note(:single_tag, "ST: reconstruct anti-Xi+ -> anti-Lambda pi+, anti-Lambda -> anti-p pi+; Xi- yield from recoil mass fit")
   .note(:blind_analysis, "10% of data used to validate procedure before unblinding full sample")
   .note(:signal_model, "LNV width distribution from Barbero et al. Phys.Rev.D 87, 036010 (2013) with B=0")

event_selection = Selection.new

event_selection.select_track {
  cos_theta 0.93
  Vz 100.0
  Vr 10.0
  nChrp ">=3"
  nChrn ">=2"
}
.select_photon {
  tdc_emc_start 0
  tdc_emc_end 14
  angle_to_track 10.0
  angle_to_prm_track 20.0
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  nGam ">=2"
}
.pid(method: :probability) {
  prob_cut 0.001
  identify :proton, against: [:kaon, :pion]
  identify :electron, against: [:pion, :kaon, :proton]
  nprp ">=1"
  nprm ">=1"
  nem ">=1"
}
.remove([:prp <= :chrgp, :prm <= :chrgn])
.remove([:em <= :chrgn])

# Reconstruct anti-Lambda -> anti-p pi+
event_selection.secondary_vertex_fit([:prm, :pip]) {
  build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
}

# Reconstruct anti-Xi+ -> anti-Lambda pi+
event_selection.secondary_vertex_fit([:Lambda_bar, :pip]) {
  build_virtual_particle(:anti_Xi_plus).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
}

# Reconstruct pi0 -> gamma gamma
event_selection.kalman_kinematic_fit([:gamma, :gamma]) {
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 30
  npi0 ">=1"
}

# Main kinematic fit with missing electron (one e- not reconstructed)
event_selection.kinematic_fit([:anti_Xi_plus, :prp, :pi0, :em]) {
  nominal
  constrain_four_momentum
  miss_track_of(:em)
  chi2_cut 20
}

alg.with_decay_card(decay_card_sig).apply(event_selection)
alg.execute_on([jpsi_data, jpsi_incMC, exMC_signal])