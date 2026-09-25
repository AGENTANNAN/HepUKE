# Measurement of Xi0 -> gamma Sigma0 weak radiative decay at J/psi
# Double Tag technique: J/psi -> Xi0 anti-Xi0
# Tag: anti-Xi0 -> anti-Lambda pi0, anti-Lambda -> anti-p pi+, pi0 -> gamma gamma
# Signal: Xi0 -> gamma Sigma0, Sigma0 -> gamma Lambda, Lambda -> p pi-
# J/psi data, ECM = 3.097 GeV
# BF = (3.69 +/- 0.21 +/- 0.12) x 10^-3
# alpha_gamma = -0.807 +/- 0.095 +/- 0.011

jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

decay_card_sig = <<~DECAYCARD
    Decay J/psi
    1.0000 Xi0 anti-Xi0 PHSP;
    Enddecay

    Decay Xi0
    1.0000 gamma Sigma0 PHSP;
    Enddecay

    Decay anti-Xi0
    1.0000 anti-Lambda0 pi0 PHSP;
    Enddecay

    Decay Sigma0
    1.0000 gamma Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- PHSP;
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
  config.sample_name = "sig_Xi0ToGammaSigma0_3097"
  config.related_dataset = jpsi_data
  config.events = 1_000_000
  config.decay_card = decay_card_sig
  config.cross_section = :default
end

alg = Algorithm.new("Xi0ToGammaSigma0")
alg.set_header(["Xi0ToGammaSigma0Alg/Xi0ToGammaSigma0.h"])
   .set_constant({ "ECMS" => [:double, 3.097] })
   .note(:double_tag, "DT: ST via anti-Xi0 -> anti-Lambda pi0 tag; signal Xi0 -> gamma Sigma0, Sigma0 -> gamma Lambda")
   .note(:polarization, "alpha_gamma parameter extracted at ROOT level from angular distribution of decay products")
   .note(:background, "dominant background from Xi0 -> Lambda pi0 with misidentified photons evaluated via MC and sidebands")

event_selection = Selection.new

event_selection.select_track {
  cos_theta 0.93
  Vz 100.0
  Vr 10.0
  nChrgp ">=2"
  nChrgn ">=2"
}
.select_photon {
  tdc_emc_start 0
  tdc_emc_end 14
  angle_to_track 10.0
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  nGam ">=4"
}
.pid(method: :probability) {
  prob_cut 0.001
  identify :proton, against: [:kaon, :pion]
  identify :pion, against: [:kaon, :proton]
  nprp ">=1"
  nprm ">=1"
  npip ">=1"
  npim ">=1"
}
.remove([:prp <= :chrgp, :prm <= :chrgn])
.remove([:pip <= :chrgp, :pim <= :chrgn])

# Reconstruct Lambda -> p pi- (signal side)
event_selection.secondary_vertex_fit([:prp, :pim]) {
  build_virtual_particle(:Lambda).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
}

# Reconstruct anti-Lambda -> anti-p pi+ (tag side)
event_selection.secondary_vertex_fit([:prm, :pip]) {
  build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
}

# Reconstruct pi0 -> gamma gamma (tag side: anti-Xi0 -> anti-Lambda pi0)
event_selection.kalman_kinematic_fit([:gamma, :gamma]) {
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 30
  npi0 ">=1"
}

# Tag: reconstruct anti-Xi0 -> anti-Lambda pi0
event_selection.secondary_vertex_fit([:Lambda_bar, :pi0]) {
  build_virtual_particle(:anti_Xi0).by_minimizing_mass_difference
  remove_used_particle_from_candidate_list
}

# Final kinematic fit: full event with all constraints
# Tag side: anti-Xi0 (reconstructed)
# Signal side: Xi0 -> gamma Sigma0, Sigma0 -> gamma Lambda
# Remaining particles: Lambda (signal-side), and two photons (radiative + from Sigma0)
event_selection.kinematic_fit([:anti_Xi0, :Lambda, :gamma, :gamma]) {
  nominal
  constrain_four_momentum
  invariant_mass_of(:gamma, :Lambda).constrain_to_nominal_mass_of(:Sigma0)
  chi2_cut 200
}

alg.with_decay_card(decay_card_sig).apply(event_selection)
alg.execute_on([jpsi_data, jpsi_incMC, exMC_signal])