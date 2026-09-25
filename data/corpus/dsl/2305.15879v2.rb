# Paper: 2305.15879v2
# Title: Amplitude analysis and branching fraction measurement of D+ -> K_S0 pi+ pi0 pi0
# Energy: 3.773 GeV (psi(3770), single energy point)
# Double-tag (DT) method:
#   ST: D- -> hadronic tag modes
#   DT: D+ -> K_S0 pi+ pi0 pi0 (signal)
# Amplitude analysis with isobar model
# K_S0 -> pi+ pi-, pi0 -> gamma gamma

### Dataset preparation ###
data_712_3773 = DatasetManager.load_real_data.find("712_3773")
incMC_712_3773 = DatasetManager.load_inclusive_mc.find("712_3773")

all_data = [data_712_3773]
all_incMC = [incMC_712_3773]

exMC_signal = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "Dp_to_KShort_piplus_pi0_pi0"
  config.events         = 500000
  config.decay_card     = nil  # TagAnalysis uses decay card via Algorithm, not exMC
end

### TagAnalysis: Double-tag with D- ST hadronic tags ###
alg = TagAnalysis.new("Dp_to_KShort_piplus_pi0_pi0")
alg.set_header(["DpKShortPiPi0Pi0Alg/DpKShortPiPi0Pi0.h"])

# Six ST tag modes for D- hadronic decays
alg.tag_side(:Dm) {
  modes(:DmtoKPiPi, :DmtoKPiPiPi0, :DmtoKSPi, :DmtoKSPiPi0,
        :DmtoKKPi, :DmtoKPiPiPi)
  charm -1
}

# Signal side: D+ -> K_S0 pi+ pi0 pi0
# K_S0 -> pi+pi-, pi0 -> gamma gamma
alg.signal_side {
  charged(pip: 1)
  photons 4  # two pi0 -> 4 photons
  require_charge 1
}

# Fit: constrain ST and DT D masses, constrain K_S0 and pi0 masses
alg.fit {
  constrain_four_momentum
  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)

  # Actually, for the amplitude analysis section (Sec 4.1), a 4C kinematic
  # fit constraining D+, K_S0, and pi0 masses is used.
  # For BF measurement (Sec 5), ST: DeltaE cut + M_BC fit; DT: DeltaE sig cut
  chi2_cut 200
}

# Additional notes specific to this analysis
alg.note(:amplitude_analysis_selection,
  "Amplitude analysis (Sec 4.1): ST mode D- -> K+pi-pi- only. DeltaE: tag [-25,20] MeV, sig [-40,20] MeV. K_S0 veto: M(pi0pi0) not in [0.46,0.52] GeV/c2. D0D0bar wrong-combination veto: M_BC^W not in [1.862,1.870] GeV. 4C kinematic fit constraining D+, K_S0, pi0 masses. 1458 events in signal region, purity 96.86%. Isobar model fit with 6 intermediate processes. Applied in ROOT.")

alg.note(:bf_measurement,
  "BF measurement (Sec 5): 6 ST tag modes. ST yield from M_BC fit. DT yield from 2D M_BC^sig vs M_BC^tag fit. K_S0 sideband subtraction. BF(D+ -> K_S0 pi+pi0pi0) = (2.888+/-0.058+/-0.069)%. DT efficiency from amplitude-model MC. Applied in ROOT.")

# 2.93 fb^-1 at psi(3770)
alg.note(:luminosity, "2.93 fb^-1 at psi(3770).")

all_datasets = all_data + all_incMC + exMC_signal
root_files = alg.execute_on([all_datasets])