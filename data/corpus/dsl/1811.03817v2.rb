# Paper: 1811.03817v2 — Dalitz plot analysis of omega → pi+ pi- pi0
# BESIII, J/psi → omega eta, 1.3×10⁹ J/psi events
# Ordinary analysis: 6C kinematic fit (4C + π0 mass + η mass)

### Dataset preparation ###
jpsi_data = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

decay_card = <<~DECAYCARD
  Decay J/psi
  1.000 omega eta PHSP;
  Enddecay

  Decay omega
  1.000 pi+ pi- pi0 OMEGA_DALITZ;
  Enddecay

  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Jpsi_omega_eta"
  config.related_dataset = jpsi_data
  config.events = 24_000_000
  config.decay_card = decay_card
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg = Algorithm.new("JpsiOmegaEta")
alg.set_header(["JpsiOmegaEtaAlg/JpsiOmegaEta.h"])
  .set_constant({ "ECMS" => [:double, 3.097] })
  .note(:background_veto, "Lambda' veto: |M(π+π-π0 γ_low) - m_etap| > 0.04 GeV/c² applied at ROOT level")
  .note(:tag_mode_unavailable, "omega signal region |M(π+π-π0) - m_omega| < 0.04 GeV/c² applied at ROOT level")

event_selection = Selection.new
event_selection.select_track {
  cos_theta   0.93
  Vz          10.0
  Vr          1.0
  nChrp       ">=1"
  nChrn       ">=1"
  nNet        "==0"
}
.select_photon {
  tdc_emc_start     0
  tdc_emc_end       14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track    10.0
  nGam              ">=4"
}
.pid(method: :probability) {
  prob_cut   0.001
  identify :pion, against: [:kaon, :proton]
  npip       ">=1"
  npim       ">=1"
}
# 6C kinematic fit: 4C + π0 mass constraint + η mass constraint
.kinematic_fit([:pip, :pim, :gamma, :gamma, :gamma, :gamma]) {
  nominal
  constrain_four_momentum
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  chi2_cut 200
}

alg.with_decay_card(decay_card).apply(event_selection)
alg.execute_on([jpsi_data, jpsi_incMC, exMC_signal])