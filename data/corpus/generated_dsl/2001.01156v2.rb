# =============================================================================
# e+e- -> gamma X(3872) on the psi(4260) energy scan (4.180, 4.230, 4.260 GeV)
#   Mode A (signal) : X(3872) -> gamma J/psi        , J/psi -> l+l-
#   Mode B (signal) : X(3872) -> D*0 anti-D0        , D*0 -> gamma D0 , D0 -> K- pi+
#   Mode C (norm.)  : X(3872) -> pi+ pi- J/psi      , J/psi -> l+l-
# =============================================================================

### Datasets ###
data_4180 = DatasetManager.real_data.find("703_4180")
data_4230 = DatasetManager.real_data.find("703_4230")
data_4260 = DatasetManager.real_data.find("703_4260")
scan_points = [data_4180, data_4230, data_4260]

incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")
incMCs = [incMC_4180, incMC_4230, incMC_4260]

### Decay cards (EvtGen) ###
decay_card_gammaJpsi = <<~DECAYCARD
    Decay psi(4260)
    1.000 gamma X(3872) PHSP;
    Enddecay
    Decay X(3872)
    1.000 gamma J/psi PHSP;
    Enddecay
    Decay J/psi
    1.000 e+ e- PHOTOS VLL;
    Enddecay
    End
DECAYCARD

decay_card_DstarD0 = <<~DECAYCARD
    Decay psi(4260)
    1.000 gamma X(3872) PHSP;
    Enddecay
    Decay X(3872)
    1.000 D*0 anti-D0 PHSP;
    Enddecay
    Decay D*0
    1.000 gamma D0 PHSP;
    Enddecay
    Decay D0
    1.000 K- pi+ PHSP;
    Enddecay
    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay
    End
DECAYCARD

decay_card_pipiJpsi = <<~DECAYCARD
    Decay psi(4260)
    1.000 gamma X(3872) PHSP;
    Enddecay
    Decay X(3872)
    1.000 pi+ pi- J/psi PHSP;
    Enddecay
    Decay J/psi
    1.000 e+ e- PHOTOS VLL;
    Enddecay
    End
DECAYCARD

### Exclusive MC: 200k events per mode, one sample per scan point ###
exMC_gammaJpsi = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "exmc_gammaJpsi"
  config.events        = 200_000
  config.decay_card    = decay_card_gammaJpsi
  config.cross_section = :default
end

exMC_DstarD0 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "exmc_DstarD0"
  config.events        = 200_000
  config.decay_card    = decay_card_DstarD0
  config.cross_section = :default
end

exMC_pipiJpsi = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "exmc_pipiJpsi"
  config.events        = 200_000
  config.decay_card    = decay_card_pipiJpsi
  config.cross_section = :default
end

# =============================================================================
# Mode A : e+e- -> gamma X(3872) -> gamma gamma J/psi, J/psi -> l+l-
# =============================================================================
alg_gammaJpsi = Algorithm.new("GammaJpsi")
alg_gammaJpsi.set_header(["GammaJpsiAlg/GammaJpsi.h"])
             .set_constant({"ECMS" => [:double, 4.260]})
             .set_alias({"std::vector<double>" => "Vdouble"})
             .note(:beam_energy_scan,
                   "ECMS differs across the three scan points (4.180 / 4.230 / 4.260 GeV); \
the ECMS constant is set to 4.260 GeV and must be re-set per energy point")
             .note(:background_veto,
                   "pi0/eta veto on M(gamma_L gamma_H): the two photons are ordered by energy \
(higher = gamma_H, lower = gamma_L); events with M(gamma gamma) inside the pi0 or eta mass window are rejected")
             .note(:bhabha_suppression,
                   "Bhabha suppression |cos theta_gamma| in [-0.7, 0.7] applied for the J/psi -> e+e- channel \
to reject e+e- -> e+e- gamma")

sel_gammaJpsi = Selection.new
  .select_track {                         # charged track selection
      cos_theta 0.93                      # |cos(theta)| < 0.93
      Vz        10.0                      # |Vz| < 10 cm
      Vr        1.0                       # Vr < 1 cm
      nNet      "==0"                     # net charge zero
      nChrp     ">=1"                     # at least one positive track
      nChrn     ">=1"                     # at least one negative track
  }
  .select_photon {                        # photon selection
      tdc_emc_start     0
      tdc_emc_end       14
      angle_to_track    10.0              # angle to nearest charged track > 10 deg
      energyThreshold_b 0.025             # barrel threshold 25 MeV
      energyThreshold_e 0.050             # endcap threshold 50 MeV
      nGam              ">=2"             # at least two photons
  }
  .pid(method: :probability) {            # PID (probability method)
      prob_cut 0.001
      identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                     treat_as_electron_if_energy_above: 0.6
      nlp ">=1"                           # at least one lepton+
      nlm ">=1"                           # at least one lepton-
  }
  .remove([:lp <= :chrgp, :lm <= :chrgn])       # remove identified leptons
  .assign({:chrgp => :pip, :chrgn => :pim})     # remaining charged tracks are pions
  .kinematic_fit([:lp, :lm, :gamma, :gamma]) {  # nominal 4C fit to l+l- gamma gamma
      nominal
      constrain_four_momentum
      chi2_cut 40                               # chi2 < 40
      invariant_mass_of(:lp, :lm).within(3.077, 3.117)          # |M(l+l-) - m_J/psi| < 0.02
      invariant_mass_of(:gamma, :lp, :lm).out_of(3.491, 3.531)  # chi_c veto |M(gamma_L J/psi) - m_chi_c| > 0.02
  }

alg_gammaJpsi.with_decay_card(decay_card_gammaJpsi).apply(sel_gammaJpsi)
root_files_gammaJpsi = alg_gammaJpsi.execute_on(scan_points + incMCs + exMC_gammaJpsi)

# =============================================================================
# Mode B : X(3872) -> D*0 anti-D0 -> gamma K+ K- pi+ pi-
# =============================================================================
alg_DstarD0 = Algorithm.new("DstarD0")
alg_DstarD0.set_header(["DstarD0Alg/DstarD0.h"])
           .set_constant({"ECMS" => [:double, 4.260]})
           .set_alias({"std::vector<double>" => "Vdouble"})
           .note(:beam_energy_scan,
                 "ECMS differs across the three scan points (4.180 / 4.230 / 4.260 GeV); \
the ECMS constant is set to 4.260 GeV and must be re-set per energy point")
           .note(:background_veto,
                 "pi0/eta veto applied on M(gamma gamma); D-meson mass windows applied on the reconstructed \
D0 / anti-D0; D*0 mass window +/-0.004 GeV for the (noted) pi0 D0 sub-mode; the smallest-chi2 combination \
is chosen among the multi-photon ambiguities")

sel_DstarD0 = Selection.new
  .select_track {                         # charged track selection
      cos_theta 0.93
      Vz        10.0
      Vr        1.0
      nNet      "==0"
      nChrp     ">=2"                     # at least two positive tracks
      nChrn     ">=2"                     # at least two negative tracks
  }
  .select_photon {                        # photon selection
      tdc_emc_start     0
      tdc_emc_end       14
      angle_to_track    10.0
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      nGam              ">=2"             # at least two photons
  }
  .pid(method: :probability) {            # PID (probability method), pi/K separation
      prob_cut 0.001
      identify :kaon, against: [:pion, :proton]
      identify :pion, against: [:kaon, :proton]
      nkp  ">=1"
      nkm  ">=1"
      npip ">=1"
      npim ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {  # intermediate gamma gamma fit to pi0 mass
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 25                            # chi2 < 25 (no pi0-count requirement)
  }
  .kinematic_fit([:kp, :km, :pip, :pim, :gamma, :gamma]) {  # nominal 4C fit to K+K-pi+pi-gamma gamma
      nominal
      constrain_four_momentum
      chi2_cut 60                                              # chi2 < 60
      invariant_mass_of(:gamma, :km, :pip).within(2.001, 2.013)  # D*0 window +/-0.006 GeV (gamma D0 mode)
  }

alg_DstarD0.with_decay_card(decay_card_DstarD0).apply(sel_DstarD0)
root_files_DstarD0 = alg_DstarD0.execute_on(scan_points + incMCs + exMC_DstarD0)

# =============================================================================
# Mode C (normalization) : X(3872) -> pi+ pi- J/psi -> gamma pi+ pi- l+ l-
# =============================================================================
alg_pipiJpsi = Algorithm.new("PipiJpsi")
alg_pipiJpsi.set_header(["PipiJpsiAlg/PipiJpsi.h"])
            .set_constant({"ECMS" => [:double, 4.260]})
            .set_alias({"std::vector<double>" => "Vdouble"})
            .note(:beam_energy_scan,
                  "ECMS differs across the three scan points (4.180 / 4.230 / 4.260 GeV); \
the ECMS constant is set to 4.260 GeV and must be re-set per energy point")
            .note(:background_veto,
                  "pi0/eta veto applied on M(gamma gamma) in the normalization mode")

sel_pipiJpsi = Selection.new
  .select_track {                         # charged track selection
      cos_theta 0.93
      Vz        10.0
      Vr        1.0
      nNet      "==0"
      nChrp     ">=2"                     # at least two positive tracks
      nChrn     ">=2"                     # at least two negative tracks
  }
  .select_photon {                        # photon selection
      tdc_emc_start     0
      tdc_emc_end       14
      angle_to_track    10.0
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      nGam              ">=1"             # at least one photon
  }
  .pid(method: :probability) {            # PID (probability method): leptons + pions
      prob_cut 0.001
      identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                     treat_as_electron_if_energy_above: 0.6
      identify :pion, against: [:kaon, :proton]
      nlp  ">=1"
      nlm  ">=1"
      npip ">=1"
      npim ">=1"
  }
  .kinematic_fit([:pip, :pim, :lp, :lm, :gamma]) {  # nominal 4C fit to pi+pi-l+l-gamma
      nominal
      constrain_four_momentum
      chi2_cut 200                                   # chi2 < 200
      invariant_mass_of(:lp, :lm).within(3.077, 3.117)  # |M(l+l-) - m_J/psi| < 0.02
  }

alg_pipiJpsi.with_decay_card(decay_card_pipiJpsi).apply(sel_pipiJpsi)
root_files_pipiJpsi = alg_pipiJpsi.execute_on(scan_points + incMCs + exMC_pipiJpsi)