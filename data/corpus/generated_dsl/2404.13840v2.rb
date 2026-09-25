# =====================================================================
#  e+e- -> omega X(3872)  and  e+e- -> gamma X(3872)
#  at 4.66, 4.68, 4.70, 4.74, 4.75, 4.78, 4.84 GeV
#  BOSS part: dataset preparation + event selection (up to final 4C fit)
# =====================================================================

### ------------------------------------------------------------------
### Dataset description (real data + inclusive MC at each energy point)
### ------------------------------------------------------------------
data_4660 = DatasetManager.real_data.find("706_4660")   # 4.661 GeV
data_4680 = DatasetManager.real_data.find("706_4680")   # 4.682 GeV
data_4700 = DatasetManager.real_data.find("706_4700")   # 4.699 GeV
data_4740 = DatasetManager.real_data.find("707_4740")   # 4.740 GeV
data_4750 = DatasetManager.real_data.find("707_4750")   # 4.750 GeV
data_4780 = DatasetManager.real_data.find("707_4780")   # 4.781 GeV
data_4840 = DatasetManager.real_data.find("707_4840")   # 4.843 GeV
data_points = [data_4660, data_4680, data_4700, data_4740, data_4750, data_4780, data_4840]

incMC_4660 = DatasetManager.inclusive_mc.find("706_4660")
incMC_4680 = DatasetManager.inclusive_mc.find("706_4680")
incMC_4700 = DatasetManager.inclusive_mc.find("706_4700")
incMC_4740 = DatasetManager.inclusive_mc.find("707_4740")
incMC_4750 = DatasetManager.inclusive_mc.find("707_4750")
incMC_4780 = DatasetManager.inclusive_mc.find("707_4780")
incMC_4840 = DatasetManager.inclusive_mc.find("707_4840")
incMC_points = [incMC_4660, incMC_4680, incMC_4700, incMC_4740, incMC_4750, incMC_4780, incMC_4840]

### ------------------------------------------------------------------
### Decay cards (EvtGen format).  Top mother = psi(4260) (KKMC default).
### J/psi -> l+l- is generated for both lepton flavours (l = e, mu).
### ------------------------------------------------------------------

# --- Mode 1: e+e- -> omega X(3872), X(3872) -> gamma J/psi -------------
decay_card_mode1 = <<~DECAYCARD
    Decay psi(4260)
    1.0000 omega X(3872) PHSP;
    Enddecay

    Decay X(3872)
    1.0000 gamma J/psi PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    Decay J/psi
    0.500 e+ e- PHOTOS VLL;
    0.500 mu+ mu- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# --- Mode 2: e+e- -> omega chi_c1, chi_c1 -> gamma J/psi ---------------
decay_card_mode2 = <<~DECAYCARD
    Decay psi(4260)
    1.0000 omega chi_c1 PHSP;
    Enddecay

    Decay chi_c1
    1.0000 gamma J/psi PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    Decay J/psi
    0.500 e+ e- PHOTOS VLL;
    0.500 mu+ mu- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# --- Mode 3: e+e- -> omega X(3872), X(3872) -> pi+ pi- J/psi -----------
decay_card_mode3 = <<~DECAYCARD
    Decay psi(4260)
    1.0000 omega X(3872) PHSP;
    Enddecay

    Decay X(3872)
    1.0000 pi+ pi- J/psi PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    Decay J/psi
    0.500 e+ e- PHOTOS VLL;
    0.500 mu+ mu- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# --- Mode 4: e+e- -> gamma X(3872), X(3872) -> pi+ pi- J/psi -----------
decay_card_mode4 = <<~DECAYCARD
    Decay psi(4260)
    1.0000 gamma X(3872) PHSP;
    Enddecay

    Decay X(3872)
    1.0000 pi+ pi- J/psi PHSP;
    Enddecay

    Decay J/psi
    0.500 e+ e- PHOTOS VLL;
    0.500 mu+ mu- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

### ------------------------------------------------------------------
### Exclusive MC  (200k events per mode, distributed over the 7 energy
### points by luminosity weight via the multi-round interface)
### ------------------------------------------------------------------
exMC_mode1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_wX3872_gamJpsi"
  config.related_dataset = data_points          # Array -> one MC per energy point
  config.events          = 200000
  config.decay_card      = decay_card_mode1
  config.cross_section   = :default
end

exMC_mode2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_wchic1_gamJpsi"
  config.related_dataset = data_points
  config.events          = 200000
  config.decay_card      = decay_card_mode2
  config.cross_section   = :default
end

exMC_mode3 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_wX3872_pipiJpsi"
  config.related_dataset = data_points
  config.events          = 200000
  config.decay_card      = decay_card_mode3
  config.cross_section   = :default
end

exMC_mode4 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_gX3872_pipiJpsi"
  config.related_dataset = data_points
  config.events          = 200000
  config.decay_card      = decay_card_mode4
  config.cross_section   = :default
end

### ==================================================================
### Event selection (BOSS)
### ==================================================================

# =====================================================================
# Channel A :  omega gamma J/psi   (Modes 1 & 2 share identical final
#              state pi+pi-gamma gamma gamma l+l- and identical cuts)
# =====================================================================
alg_wgamJpsi = Algorithm.new("WGammaJpsi")
alg_wgamJpsi.set_header(["WGammaJpsiAlg/WGammaJpsi.h"])
            .set_constant({"ECMS" => [:double, 4.68]})   # nominal CMS energy (per-point energy read at run time)
            .set_alias({"std::vector<double>" => "Vdouble"})

sel_wgamJpsi = Selection.new
  .select_track {                 # charged-track quality cuts
    cos_theta 0.93                # |cos(theta)| < 0.93
    Vz        10.0                # |Vz| < 10 cm
    Vr        1.0                 # Vr < 1 cm
    nChrp     ">=2"               # >= 2 positive tracks
    nChrn     ">=2"               # >= 2 negative tracks
  }
  .select_photon {                # photon selection
    tdc_emc_start     0           # EMC timing window 0 - 700 ns
    tdc_emc_end       14
    energyThreshold_b 0.025       # E > 25 MeV (barrel)
    energyThreshold_e 0.050       # E > 50 MeV (endcap)
    angle_to_track    10.0        # >= 10 deg from any charged track
    nGam              ">=2"       # at least two photons for omega gamma J/psi
  }
  .pid(method: :probability) {    # probability-based PID
    prob_cut 0.001                # prob > 0.001
    # p > 1 GeV tracks are treated as leptons (p < 1 GeV remain pions); the
    # lepton is an electron if EMC eraw is above threshold, otherwise a muon.
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 1.0
    identify :pion, against: [:kaon]   # pi/K separation (pi+/pi-)
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {        # pi0 : gamma gamma mass-constrained fit
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0     ">=1"
  }
  .kinematic_fit([:pip, :pim, :pi0, :gamma, :lp, :lm]) {   # nominal 4C fit
    nominal
    constrain_four_momentum
    chi2_cut 200
    invariant_mass_of(:pip, :pim, :pi0).within(0.74, 0.82)      # omega mass window
    invariant_mass_of(:lp, :lm).within(3.06, 3.14)              # J/psi mass window
    invariant_mass_of(:pip, :pim).out_of(3.656, 3.716)          # psi(3686) veto on M(pi+pi-)
    invariant_mass_of(:pip, :pim, :lp, :lm).out_of(3.676, 3.696) # psi(3686) veto on M(pi+pi-J/psi)
  }

alg_wgamJpsi
  .note(:muon_muc_depth, "require at least one muon candidate with MUC hit depth > 30 cm; " \
        "MUC hit depth is not expressible in the current DSL and is applied as an extra " \
        "PID refinement on the high-momentum lepton list (index_lp/index_lm)")
  .note(:lepton_emc_windows, "electron (muon) lepton hypothesis requires EMC energy > 1.0 GeV " \
        "(< 0.4 GeV); the DSL identify_high_momentum_leptons exposes a single EMC threshold, " \
        "so the exact two-sided e/mu window is applied during the PID correction step")
  .note(:ecms_per_point, "analysis runs over 7 energy points (4.66-4.84 GeV); the ECMS constant " \
        "is set to a nominal value and the true per-run beam energy is used for the 4C constraint")

alg_wgamJpsi.with_decay_card(decay_card_mode1).apply(sel_wgamJpsi)
root_files_wgamJpsi = alg_wgamJpsi.execute_on(data_points + incMC_points + exMC_mode1 + exMC_mode2)

# =====================================================================
# Channel B :  omega pi+ pi- J/psi
#              final state pi+ pi- pi+ pi- gamma gamma l+ l-
# =====================================================================
alg_wpipiJpsi = Algorithm.new("WPipiJpsi")
alg_wpipiJpsi.set_header(["WPipiJpsiAlg/WPipiJpsi.h"])
             .set_constant({"ECMS" => [:double, 4.68]})
             .set_alias({"std::vector<double>" => "Vdouble"})

sel_wpipiJpsi = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=2"
    nChrn     ">=2"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=1"       # at least one photon for the other modes
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 1.0
    identify :pion, against: [:kaon]
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {        # pi0 mass-constrained fit
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0     ">=1"
  }
  .kinematic_fit([:pip, :pim, :pip, :pim, :pi0, :lp, :lm]) {   # nominal 4C fit
    nominal
    constrain_four_momentum
    chi2_cut 200
    invariant_mass_of(:pip, :pim, :pi0).within(0.75, 0.81)      # omega mass window
    invariant_mass_of(:lp, :lm).within(3.07, 3.13)              # J/psi mass window
  }

alg_wpipiJpsi
  .note(:topology_coverage, "selection covers both the five-track (pi+pi-pi+-l+l-gamma gamma) " \
        "and the six-track (pi+pi-pi+pi-l+l-gamma) topologies; the number of reconstructed " \
        "pions/photons is accommodated by the combinatorial candidate sharing in the 4C fit")

alg_wpipiJpsi.with_decay_card(decay_card_mode3).apply(sel_wpipiJpsi)
root_files_wpipiJpsi = alg_wpipiJpsi.execute_on(data_points + incMC_points + exMC_mode3)

# =====================================================================
# Channel C :  gamma pi+ pi- J/psi
#              final state gamma pi+ pi- l+ l-
# =====================================================================
alg_gpipiJpsi = Algorithm.new("GPipiJpsi")
alg_gpipiJpsi.set_header(["GPipiJpsiAlg/GPipiJpsi.h"])
             .set_constant({"ECMS" => [:double, 4.68]})
             .set_alias({"std::vector<double>" => "Vdouble"})

sel_gpipiJpsi = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=2"
    nChrn     ">=2"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=1"       # at least one photon for the other modes
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 1.0
    identify :pion, against: [:kaon]
  }
  .kinematic_fit([:pip, :pim, :gamma, :lp, :lm]) {   # nominal 4C fit
    nominal
    constrain_four_momentum
    chi2_cut 200
    invariant_mass_of(:lp, :lm).within(3.08, 3.12)   # J/psi mass window
  }

alg_gpipiJpsi.with_decay_card(decay_card_mode4).apply(sel_gpipiJpsi)
root_files_gpipiJpsi = alg_gpipiJpsi.execute_on(data_points + incMC_points + exMC_mode4)