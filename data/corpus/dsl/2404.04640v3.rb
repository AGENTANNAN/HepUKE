# Dataset description
jpsi_data = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Decay card for signal: J/psi -> gamma a, a -> gamma gamma
decay_card = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma a PHSP;
  Enddecay
  Decay a
  1.0000 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# Exclusive MC for signal process
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "jpsi_gamma_a_gammagamma"
  config.related_dataset = jpsi_data
  config.events = 1000000
  config.decay_card = decay_card
  config.cross_section = :default
end

# Algorithm
alg_name = "JpsiToGammaALP"
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
my_Algorithm.set_constant({"ECMS" => [:double, 3.097]})

# Event selection: purely neutral final state
event_selection = Selection.new
event_selection.select_track {
  nChrp "==0"   # Zero positively charged tracks
  nChrn "==0"   # Zero negatively charged tracks
}
.select_photon {
  tdc_emc_start 0
  tdc_emc_end 14
  energyThreshold_b 0.025       # Barrel: |cosTheta| < 0.80, min 25 MeV
  energyThreshold_e 1000.0      # Disable endcap (barrel-only reconstruction)
  angle_to_track 10.0
  nGam ">=3"                    # At least 3 photon candidates
}
# Nominal 4C kinematic fit: 3 photons constrained to CMS energy
.kinematic_fit([:gamma, :gamma, :gamma]) {
  nominal
  constrain_four_momentum
  chi2_cut 30                   # Require chi2_4C < 30
}
# Competing hypothesis: 2-photon fit (to veto e+e- -> gamma gamma)
.kinematic_fit([:gamma, :gamma]) {
  constrain_four_momentum
}
# Competing hypothesis: 4-photon fit (to veto J/psi -> gamma pi0 pi0 etc.)
.kinematic_fit([:gamma, :gamma, :gamma, :gamma]) {
  constrain_four_momentum
}
# Competing hypothesis: 5-photon fit
.kinematic_fit([:gamma, :gamma, :gamma, :gamma, :gamma]) {
  constrain_four_momentum
}

my_Algorithm
  .note(:photon_time_diff, "EMC time difference between any two photons: -500 < DeltaT < 500 ns to suppress electronics noise")
  .note(:chi2_comparison, "require chi2_4C_3gamma < chi2_4C_2gamma AND chi2_4C_3gamma < chi2_4C_ngamma (n=4,5) to suppress QED and J/psi -> gamma pi0 pi0 backgrounds")
  .note(:deltaE_cut, "energy difference between 1st(2nd) and 3rd photon: DeltaE_13 < 1.46 GeV, DeltaE_23 < 1.41 GeV; photons sorted by decreasing energy")
  .note(:deltaphi_cut, "absolute azimuthal angle difference between 3rd and 1st photon: Deltaphi_31 > 1 radian")
  .note(:background_veto, "reject events where any di-photon invariant mass falls in pi0(0.11-0.16), eta(0.52-0.56), etap(0.92-0.99), or etac(2.92-3.04) GeV/c^2 windows (J/psi -> gamma P suppression)")
  .with_decay_card(decay_card)
  .apply(event_selection)

root_files = my_Algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])