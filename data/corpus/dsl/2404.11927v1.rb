# Paper 2404.11927v1: Search for e+e- -> K+K- psi(3770) at 4.84-4.95 GeV
# Reconstructs K+K- + D meson; recoiling mass method

# --- Datasets (BOSS 707, 4.84-4.95 GeV) ---
data_4840 = DatasetManager.real_data.find("707_4840")
data_4914 = DatasetManager.real_data.find("707_4914")
data_4946 = DatasetManager.real_data.find("707_4946")
incMC_4840 = DatasetManager.inclusive_mc.find("707_4840")
incMC_4914 = DatasetManager.inclusive_mc.find("707_4914")
incMC_4946 = DatasetManager.inclusive_mc.find("707_4946")

data_points = [data_4840, data_4914, data_4946]

# --- Decay card: e+e- -> K+K- D Dbar (psi(3770) -> D Dbar) ---
# Signal: K+K- + D meson reconstruction
decay_card = <<~DECAYCARD
    Decay psi(4260)
    1.000 K+ K- D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 K- pi+ PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# --- Exclusive MC ---
exMC = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "sig_KKpsi3770"
  config.events = 200000
  config.decay_card = decay_card
  config.cross_section = :default
end

# ==========================================================================
# Algorithm: e+e- -> K+K- psi(3770), psi(3770)-> D Dbar
# Reconstruct D meson + K+K- for recoil mass analysis
# ==========================================================================
alg = Algorithm.new("KKpsi3770Search")
alg.set_header(["KKpsi3770SearchAlg/KKpsi3770Search.h"])
   .set_constant({ "ECMS" => [:double, 4.840] })
   .note(:multi_energy, "Data taken at three energy points: 4.84, 4.91, 4.95 GeV; run separately per dataset")
   .note(:recoil_mass_method, "Signal yield extracted from 2D fit to RM(K+K-D) vs RM(K+K-); D presence inferred from recoiling mass")
   .note(:phi_veto, "Events with |M(K+K-)-M_phi| < 0.02 GeV/c^2 vetoed to suppress phi backgrounds")
   .note(:d_tag_modes, "Nine D decay modes used: D0->K+pi-, D0->K+pi-pi0, D0->K+pi+pi-pi-; D-->K+pi-pi-, D-->K+pi-pi-pi0, D-->KSpi-, D-->KSpipi0, D-->KSpi-pi-pi+, D-->K+K-pi-")
   .note(:best_d_candidate, "If multiple D candidates, the one with mass closest to nominal D mass is selected")
   .note(:d_1c_kinfit, "1C kinematic fit constraining Dbar mass to nominal value, chi2 < 13; optimized by Punzi FOM")
   .note(:lowest_momentum_kaon, "Bachelor kaons are the lowest-momentum kaon pair not used in D reconstruction")

event_selection = Selection.new

event_selection.select_track {
  cos_theta 0.93
  Vz 10.0
  Vr 1.0
  nChrp ">=2"
  nChrn ">=2"
}
.select_photon {
  tdc_emc_start 0
  tdc_emc_end 14
  angle_to_track 10.0
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  nGam ">=2"
}
.pid(method: :probability) {
  prob_cut 0.001
  identify :kaon, against: [:pion]
  identify :pion, against: [:kaon]
}
.kalman_kinematic_fit([:gamma, :gamma]) {
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  chi2_cut 25
  npi0 ">=0"
}
.kinematic_fit([:kp, :km, :pip, :pim, :gamma, :gamma]) {
  nominal
  constrain_four_momentum
  chi2_cut 200
}

alg.with_decay_card(decay_card).apply(event_selection)
alg.execute_on(data_points + [incMC_4840, incMC_4914, incMC_4946] + exMC)