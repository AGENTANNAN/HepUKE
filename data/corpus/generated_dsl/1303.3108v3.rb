### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi real data at 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # Corresponding inclusive J/psi MC

# Signal decay card: J/psi -> omega p pbar, omega -> gamma pi0, pi0 -> gamma gamma (phase space)
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.000 omega p+ anti-p- PHSP;
    Enddecay

    Decay omega
    1.000 gamma pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Irreducible background decay card: J/psi -> gamma pi0 p pbar (phase space);
# its p pbar shape serves as the non-resonant template
decay_card_bkg = <<~DECAYCARD
    Decay J/psi
    1.000 gamma pi0 p+ anti-p- PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 500k-event phase-space exclusive MC for the signal chain
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_omega_ppbar_phsp"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

# 500k-event phase-space exclusive MC for the irreducible background
exMC_bkg = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_gamma_pi0_ppbar_phsp"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_bkg
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "OmegaPPbar"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.097]})   # CMS energy = 3.097 GeV (J/psi peak)
   .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {          # Charged track selection
    cos_theta 0.8                       # |cos(theta)| < 0.8
    Vz        10.0                      # |Vz| < 10 cm
    Vr        1.0                       # Vr < 1 cm in transverse plane
    nChrp     "==1"                     # Exactly 1 positive track
    nChrn     "==1"                     # Exactly 1 negative track
    nNet      "==0"                     # Net charge zero
  }
  .pid(method: :probability) {          # Probability PID
    prob_cut 0.001                      # Confidence-level cut 0.001
    identify :proton, against: [:pion, :kaon]  # p and pbar vs pi/K
    nprp "==1"                          # Exactly 1 proton
    nprm "==1"                          # Exactly 1 antiproton
  }
  .select_photon {                      # Photon selection
    tdc_emc_start 0                     # EMC timing window 0
    tdc_emc_end 14                      # EMC timing window 14
    angle_to_track 10.0                 # Min angle to nearest charged track (deg)
    energyThreshold_b 0.025             # E > 25 MeV in barrel (|cos(theta)|<0.8)
    energyThreshold_e 0.050             # E > 50 MeV in endcap (0.86<|cos(theta)|<0.92)
    nGam ">=3"                          # At least 3 photons
  }
  .select_isolated_photon {             # Each photon > 30 deg from p and pbar
    angle_to_prp_track 30.0             # Min angle to nearest proton track (deg)
    angle_to_prm_track 30.0             # Min angle to nearest antiproton track (deg)
    nGam ">=3"                          # At least 3 photons remain
  }
  .remove(:prp) { condition "three_momentum_of(:prp) < 0.3" }  # Drop p with p < 0.3 GeV/c
  .remove(:prm) { condition "three_momentum_of(:prm) < 0.3" }  # Drop pbar with p < 0.3 GeV/c
  .kinematic_fit([:gamma, :gamma, :gamma, :prp, :prm]) {       # Nominal 4C fit to gamma gamma gamma p pbar
    nominal                             # Mark as the nominal fit
    constrain_four_momentum             # 4C energy-momentum conservation; min chi2 combination kept
    chi2_cut 200                        # chi2 < 200 as encoded
  }

# Attach the signal decay card and render the selection
alg.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on real data, inclusive MC, signal exclusive MC, and background exclusive MC
root_files = alg.execute_on([jpsi_data, jpsi_incMC, exMC_signal, exMC_bkg])