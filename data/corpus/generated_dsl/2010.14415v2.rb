# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# Six c.m. energy points: 4.178, 4.226, 4.258, 4.358, 4.416, 4.600 GeV
data_4178 = DatasetManager.real_data.find("703_4180")   # 4.178 GeV
data_4226 = DatasetManager.real_data.find("703_4230")   # 4.226 GeV (sample label 4230)
data_4258 = DatasetManager.real_data.find("703_4260")   # 4.258 GeV (sample label 4260)
data_4358 = DatasetManager.real_data.find("703_4360")   # 4.358 GeV (sample label 4360)
data_4416 = DatasetManager.real_data.find("703_4420")   # 4.416 GeV (sample label 4420)
data_4600 = DatasetManager.real_data.find("703_4600")   # 4.600 GeV

incMC_4178 = DatasetManager.inclusive_mc.find("703_4180")
incMC_4226 = DatasetManager.inclusive_mc.find("703_4230")
incMC_4258 = DatasetManager.inclusive_mc.find("703_4260")
incMC_4358 = DatasetManager.inclusive_mc.find("703_4360")
incMC_4416 = DatasetManager.inclusive_mc.find("703_4420")
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")

scan_data  = [data_4178, data_4226, data_4258, data_4358, data_4416, data_4600]
scan_incMC = [incMC_4178, incMC_4226, incMC_4258, incMC_4358, incMC_4416, incMC_4600]

# Decay card for the signal mode: e+e- -> eta_c pi+ pi- pi0,
# eta_c -> K_S0 K+ pi-, K_S0 -> pi+ pi-, pi0 -> gamma gamma
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.000  eta_c pi+ pi- pi0    PHSP;
    Enddecay

    Decay eta_c
    1.000  K_S0 K+ pi-    PHSP;
    Enddecay

    Decay K_S0
    1.000  pi+ pi-    PHSP;
    Enddecay

    Decay pi0
    1.000  gamma gamma    PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for the upper-limit channel: e+e- -> eta_c pi+ pi-
# (same eta_c and K_S0 decays)
decay_card_upper_limit = <<~DECAYCARD
    Decay psi(4260)
    1.000  eta_c pi+ pi-    PHSP;
    Enddecay

    Decay eta_c
    1.000  K_S0 K+ pi-    PHSP;
    Enddecay

    Decay K_S0
    1.000  pi+ pi-    PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive MC for the signal mode at each of the six energy points
exMCs_signal = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name   = "exmc_etac_pipipi0"
  config.events        = 100_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

# 100k-event exclusive MC for the upper-limit channel at each of the six energy points
exMCs_upper_limit = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name   = "exmc_etac_pipi"
  config.events        = 100_000
  config.decay_card    = decay_card_upper_limit
  config.cross_section = :default
end

### Event selection (BOSS) ###

# ------------------------------------------------------------------
# Mode I: e+e- -> eta_c pi+ pi- pi0 (eta_c -> K_S0 K+ pi-, K_S0 -> pi+ pi-, pi0 -> gamma gamma)
# ------------------------------------------------------------------
alg_name_modeI = "EtacPiPiPi0"
alg_modeI = Algorithm.new(alg_name_modeI)
alg_modeI.set_header(["#{alg_name_modeI}Alg/#{alg_name_modeI}.h"])
         .set_constant({ "ECMS" => [:double, 4.178] })  # nominal value; the measured per-run beam energy is used across the scan
         .set_alias({ "std::vector<double>" => "Vdouble" })
         .note(:background_veto, "the eta_c pi+ pi- pi0 mode is fitted simultaneously in 16 hadronic eta_c decay channels; D-meson / K* / omega / eta vetoes are optimised per channel at the BOSS level")
         .note(:efficiency_curve, "ISR corrections are applied for the eta_c pi+ pi- pi0 channel in the simultaneous per-channel fit")

sel_modeI = Selection.new
  .select_track {                       # Charged track selection
    cos_theta 0.93                      # |cos(theta)| < 0.93
    Vz        10.0                      # |Vz| < 10 cm
    Vr        1.0                       # Vr < 1 cm
    nChrp     ">=2"                     # at least two positive tracks
    nChrn     ">=1"                     # at least one negative track
  }
  .select_photon {                      # Photon selection
    tdc_emc_start     0                 # TDC window [0, 14]
    tdc_emc_end       14
    energyThreshold_b 0.025             # barrel  E > 25 MeV
    energyThreshold_e 0.050             # endcap  E > 50 MeV
    angle_to_track    10.0              # > 10 deg from any charged track
    nGam              ">=2"             # at least two photons
  }
  .pid(method: :probability) {          # PID with probability method
    prob_cut 0.001                      # prob > 0.001
    identify :kaon, against: [:pion, :proton]   # K+/pi-/p separation
    identify :pion, against: [:kaon, :proton]
    nkp   "==1"                         # exactly one K+
    npim  ">=1"                         # at least one pi-
  }
  .remove([:kp <= :chrgp])              # remove overlapping K+ candidates from the positive charged list
  .secondary_vertex_fit([:pip, :pim]) { # Reconstruct K_S0 from pi+ pi-
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # Reconstruct pi0 from two photons
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0     ">=1"                      # at least one pi0
  }
  .kinematic_fit([:K_S0, :kp, :pip, :pim, :pim, :pi0]) {   # nominal 4C fit to K_S0 K+ pi- pi+ pi- pi0
    nominal
    constrain_four_momentum
    chi2_cut 200                        # loose chi2 cut (tight cut applied later in ROOT)
  }

alg_modeI.with_decay_card(decay_card_signal).apply(sel_modeI)
alg_modeI.execute_on(scan_data + scan_incMC + exMCs_signal)

# ------------------------------------------------------------------
# Mode II: e+e- -> eta_c pi+ pi- (upper-limit channel, same eta_c / K_S0 decays)
# ------------------------------------------------------------------
alg_name_modeII = "EtacPiPi"
alg_modeII = Algorithm.new(alg_name_modeII)
alg_modeII.set_header(["#{alg_name_modeII}Alg/#{alg_name_modeII}.h"])
          .set_constant({ "ECMS" => [:double, 4.178] })  # nominal value; measured per-run beam energy used across the scan
          .set_alias({ "std::vector<double>" => "Vdouble" })
          .note(:z_c_search, "for the Z_c search at 4.23 GeV an additional eta_c (-> K_S0 K+ pi-) invariant-mass window 2.880 < M(eta_c) < 3.080 GeV/c2 is required, and the Z_c mass is scanned over 3625-3805 MeV/c2 and its width over 8-38 MeV in the ROOT fit")

sel_modeII = Selection.new
  .select_track {                       # same charged track selection as Mode I
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=2"
    nChrn     ">=1"
  }
  .pid(method: :probability) {          # same PID selection as Mode I
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
    nkp   "==1"
    npim  ">=1"
  }
  .remove([:kp <= :chrgp])              # remove overlapping K+ candidates from the positive charged list
  .secondary_vertex_fit([:pip, :pim]) { # same K_S0 reconstruction as Mode I
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:K_S0, :kp, :pip, :pim, :pim]) {   # nominal 4C fit to K_S0 K+ pi- pi+ pi-
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_modeII.with_decay_card(decay_card_upper_limit).apply(sel_modeII)
alg_modeII.execute_on(scan_data + scan_incMC + exMCs_upper_limit)