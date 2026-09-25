### Dataset description ###
# J/psi(3097) real data and the corresponding inclusive MC (225.2M J/psi events)
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# ---------------- Decay cards (EvtGen format) ----------------
# Signal: J/psi -> p pbar
decay_card_ppbar = <<~DECAYCARD
    Decay J/psi
    1.0000 p+ anti-p- PHSP;
    Enddecay
    End
DECAYCARD

# Signal: J/psi -> n nbar
decay_card_nnbar = <<~DECAYCARD
    Decay J/psi
    1.0000 n0 anti-n0 PHSP;
    Enddecay
    End
DECAYCARD

# Background: J/psi -> gamma eta_c, eta_c -> p pbar
decay_card_bkg_gam_etac_ppbar = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta_c PHSP;
    Enddecay

    Decay eta_c
    1.0000 p+ anti-p- PHSP;
    Enddecay
    End
DECAYCARD

# Background: J/psi -> gamma eta_c, eta_c -> n nbar
decay_card_bkg_gam_etac_nnbar = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta_c PHSP;
    Enddecay

    Decay eta_c
    1.0000 n0 anti-n0 PHSP;
    Enddecay
    End
DECAYCARD

# Background: J/psi -> pi0 p pbar (pi0 -> gamma gamma)
decay_card_bkg_pi0_ppbar = <<~DECAYCARD
    Decay J/psi
    1.0000 pi0 p+ anti-p- PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

# Background: J/psi -> pi0 n nbar (pi0 -> gamma gamma)
decay_card_bkg_pi0_nnbar = <<~DECAYCARD
    Decay J/psi
    1.0000 pi0 n0 anti-n0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

# ---------------- Exclusive MC samples ----------------
# Signal MC: J/psi -> p pbar (1,000,000 events)
exMC_ppbar = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_ppbar"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_ppbar
  config.cross_section   = :default
end

# Signal MC: J/psi -> n nbar (1,000,000 events)
exMC_nnbar = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_nnbar"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_nnbar
  config.cross_section   = :default
end

# Background MC: J/psi -> gamma eta_c, eta_c -> p pbar (200,000 events)
exMC_bkg_gam_etac_ppbar = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_gametac_ppbar"
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_bkg_gam_etac_ppbar
  config.cross_section   = :default
end

# Background MC: J/psi -> gamma eta_c, eta_c -> n nbar (200,000 events)
exMC_bkg_gam_etac_nnbar = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_gametac_nnbar"
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_bkg_gam_etac_nnbar
  config.cross_section   = :default
end

# Background MC: J/psi -> pi0 p pbar (200,000 events)
exMC_bkg_pi0_ppbar = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_pi0_ppbar"
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_bkg_pi0_ppbar
  config.cross_section   = :default
end

# Background MC: J/psi -> pi0 n nbar (200,000 events)
exMC_bkg_pi0_nnbar = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_pi0_nnbar"
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_bkg_pi0_nnbar
  config.cross_section   = :default
end

########################################################################
### Channel 1: J/psi -> p pbar                                        ###
########################################################################
alg_name_ppbar = "JpsiToPPbar"
alg_ppbar = Algorithm.new(alg_name_ppbar)
alg_ppbar.set_header(["#{alg_name_ppbar}Alg/#{alg_name_ppbar}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})   # J/psi(3097) centre-of-mass energy in GeV

selection_ppbar = Selection.new
  .select_track {
    cos_theta 0.80        # |cos(theta)| < 0.80
    Vz        10.0        # |Vz| < 10 cm
    Vr        1.0         # Vr < 1 cm
    nChrp     "==1"       # exactly one positive track
    nChrn     "==1"       # exactly one negative track
    nNet      "==0"       # net charge zero
  }
  # PID applied to the positive track only: proton hypothesis over pi and K
  .pid(method: :probability) {
    prob_cut 0.001
    identify :prp, against: [:pion, :kaon]
  }
  # The negative track is taken as the antiproton without any PID
  .assign({:chrgn => :prm})
  # Vertex fit to the two tracks, chi2 < 200, flagged as the nominal fit
  .kinematic_fit([:prp, :prm]) {
    nominal
    vertex_fit([0, 1])    # constrain prp (index 0) and prm (index 1) to a common vertex
    chi2_cut 200
  }
  # NOTE: the opening-angle cut (>178 deg) and the momentum requirement
  # (|p - 1.232| < 0.030 GeV/c) are applied at ROOT level, after the fit.

alg_ppbar.with_decay_card(decay_card_ppbar).apply(selection_ppbar)
# Same data / inclusive-MC samples, plus the p pbar signal and both backgrounds
root_files_ppbar = alg_ppbar.execute_on([
  jpsi_data,
  jpsi_incMC,
  exMC_ppbar,
  exMC_bkg_gam_etac_ppbar,
  exMC_bkg_pi0_ppbar
])

########################################################################
### Channel 2: J/psi -> n nbar                                        ###
########################################################################
alg_name_nnbar = "JpsiToNNbar"
alg_nnbar = Algorithm.new(alg_name_nnbar)
alg_nnbar.set_header(["#{alg_name_nnbar}Alg/#{alg_name_nnbar}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})

selection_nnbar = Selection.new
  .select_track {
    nChrp "==0"           # zero good charged tracks
    nChrn "==0"
    nNet  "==0"
  }
  .select_photon {
    tdc_emc_start     0    # TDC window 0-14
    tdc_emc_end       14
    angle_to_track    10.0 # photon-track angle > 10 degrees
    energyThreshold_b 0.025 # 25 MeV barrel threshold
    energyThreshold_e 0.050 # 50 MeV endcap threshold
    nGam ">=2"             # at least two showers (nbar and n candidates)
  }
# The nbar / n identification from the EMC shower properties and the E_extra = 0
# requirement, as well as the deliberate absence of any kinematic fit or
# recoil-mass constraint, cannot be expressed with the current DSL, so they are
# preserved as notes instead of being silently dropped.
alg_nnbar
  .note(:neutron_shower_selection,
        "nbar candidate = most energetic EMC shower with 0.6 < E < 2.0 GeV, " \
        "|cos(theta)| < 0.8, second moment > 20 cm^2 and > 40 hits inside a 50 degree " \
        "cone; n candidate = opposite-side shower with 0.06 < E < 0.6 GeV; " \
        "E_extra = 0 required (no additional showers beyond the nbar/n pair). " \
        "Shower second moment and hit multiplicity are not expressible in the DSL.")
  .note(:no_kinematic_fit,
        "No kinematic fit and no recoil-mass constraint are applied in this channel; " \
        "the J/psi -> n nbar signal is extracted at ROOT level from the " \
        "nbar-n opening-angle distribution near 180 degrees.")

alg_nnbar.with_decay_card(decay_card_nnbar).apply(selection_nnbar)
# Same data / inclusive-MC samples, plus the n nbar signal and both backgrounds
root_files_nnbar = alg_nnbar.execute_on([
  jpsi_data,
  jpsi_incMC,
  exMC_nnbar,
  exMC_bkg_gam_etac_nnbar,
  exMC_bkg_pi0_nnbar
])