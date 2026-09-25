# =============================================================================
# T_ccbar(4020)^- : three channels at sqrt(s) = 4.395 and 4.416 GeV
#   (1) e+e- -> D*0 D*- pi+   : partial reconstruction (pi- from D*- -> anti-D0 pi-
#                               is not reconstructed); D*0 -> D0 pi0 / D0 gamma,
#                               D0 -> K- pi+, K- pi+ pi0, K- pi+ pi+ pi-
#   (2) e+e- -> pi+ pi- J/psi : J/psi -> e+e- or mu+mu-
#   (3) e+e- -> pi+ pi- h_c   : h_c -> gamma eta_c, eta_c -> 16 hadronic modes
# BOSS part only: dataset preparation + event selection up to and including
# the final kinematic fit.
# =============================================================================

### ------------------------------- Datasets ------------------------------- ###
data_4395  = DatasetManager.real_data.find("703_4390")       # real data, sqrt(s) ~ 4.395 GeV
incMC_4395 = DatasetManager.inclusive_mc.find("703_4390")    # matching inclusive MC
data_4416  = DatasetManager.real_data.find("703_4420")       # real data, sqrt(s) ~ 4.416 GeV
incMC_4416 = DatasetManager.inclusive_mc.find("703_4420")    # matching inclusive MC

energy_points = [data_4395, data_4416]                        # both energy points
all_inputs    = [data_4395, incMC_4395, data_4416, incMC_4416]

### ----------------------------- Decay cards ------------------------------ ###
# ---- Channel 1 : e+e- -> D*0 D*- pi+ (pi- from D*- missed) -----------------
decay_card_ch1 = <<~DECAYCARD
    Decay psi(4260)
    1.0000 D*0 D*- pi+ PHSP;
    Enddecay

    Decay D*0
    0.5000 D0 pi0 PHSP;
    0.5000 D0 gamma PHSP;
    Enddecay

    Decay D*-
    1.0000 anti-D0 pi- PHSP;
    Enddecay

    Decay D0
    0.3400 K- pi+ PHSP;
    0.3300 K- pi+ pi0 PHSP;
    0.3300 K- pi+ pi+ pi- PHSP;
    Enddecay

    Decay anti-D0
    0.3400 K+ pi- PHSP;
    0.3300 K+ pi- pi0 PHSP;
    0.3300 K+ pi- pi+ pi- PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# ---- Channel 2 : e+e- -> pi+ pi- J/psi, J/psi -> e+e- ----------------------
decay_card_ch2_ee = <<~DECAYCARD
    Decay psi(4260)
    1.0000 pi+ pi- J/psi PHSP;
    Enddecay

    Decay J/psi
    1.0000 e+ e- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# ---- Channel 2 : e+e- -> pi+ pi- J/psi, J/psi -> mu+mu- --------------------
decay_card_ch2_mumu = <<~DECAYCARD
    Decay psi(4260)
    1.0000 pi+ pi- J/psi PHSP;
    Enddecay

    Decay J/psi
    1.0000 mu+ mu- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# ---- Channel 3 : e+e- -> pi+ pi- h_c, h_c -> gamma eta_c, eta_c -> hadrons --
# The 16 hadronic eta_c modes are represented by their 4-, 6- and 8-charged-track
# final states (K_S0, pi0 and eta are subsequently reconstructed).
decay_card_ch3 = <<~DECAYCARD
    Decay psi(4260)
    1.0000 pi+ pi- h_c PHSP;
    Enddecay

    Decay h_c
    1.0000 gamma eta_c PHSP;
    Enddecay

    Decay eta_c
    0.0625 K+ K- pi+ pi- PHSP;
    0.0625 K+ K- pi0 PHSP;
    0.0625 pi+ pi- pi+ pi- PHSP;
    0.0625 K_S0 K+ pi- PHSP;
    0.0625 K_S0 K- pi+ PHSP;
    0.0625 K+ K- pi+ pi- pi0 PHSP;
    0.0625 K+ K- pi+ pi- pi+ pi- PHSP;
    0.0625 pi+ pi- pi+ pi- pi+ pi- PHSP;
    0.0625 K+ K- K+ K- PHSP;
    0.0625 K+ K- pi+ pi- pi0 pi0 PHSP;
    0.0625 K_S0 K+ pi- pi0 PHSP;
    0.0625 K_S0 K_S0 PHSP;
    0.0625 pi+ pi- pi0 PHSP;
    0.0625 eta pi+ pi- PHSP;
    0.0625 eta pi+ pi- pi+ pi- PHSP;
    0.0625 eta K+ K- PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

### ------------------- Exclusive MC (500k events / channel) ---------------- ###
exMCs_ch1 = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "exmc_Dstar0DstarmPi"        # one sample per energy point
  config.events        = 500_000
  config.decay_card    = decay_card_ch1
  config.cross_section = :default
end

exMCs_ch2_ee = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "exmc_pipiJpsi_ee"
  config.events        = 500_000
  config.decay_card    = decay_card_ch2_ee
  config.cross_section = :default
end

exMCs_ch2_mumu = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "exmc_pipiJpsi_mumu"
  config.events        = 500_000
  config.decay_card    = decay_card_ch2_mumu
  config.cross_section = :default
end

exMCs_ch3 = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "exmc_pipihc_gamma_etac"
  config.events        = 500_000
  config.decay_card    = decay_card_ch3
  config.cross_section = :default
end

### ========================= Event selection (BOSS) ======================== ###

# -----------------------------------------------------------------------------
# Common charged-track and photon requirements:
#   tracks : |cos(theta)| < 0.93, Vr < 1 cm, |Vz| < 10 cm, net charge 0
#   photons: E > 25 MeV (barrel) / 50 MeV (endcap), TDC 0-14, > 10 deg from any track
# -----------------------------------------------------------------------------

# ============================== Channel 1 ====================================
alg_name_ch1 = "Tcc4020DstarDstarmPi"
alg_ch1 = Algorithm.new(alg_name_ch1)
alg_ch1.set_header(["#{alg_name_ch1}Alg/#{alg_name_ch1}.h"])
       .set_constant({"ECMS" => [:double, 4.395]})

sel_ch1 = Selection.new
sel_ch1.select_track {
          cos_theta 0.93       # |cos(theta)| < 0.93
          Vr        1.0        # Vr < 1 cm
          Vz        10.0       # |Vz| < 10 cm
          nChrp     ">=2"      # >= 1 K+ and >= 1 pi+
          nChrn     ">=2"      # >= 1 K- and >= 1 pi-
          nNet      "==0"      # zero net charge
        }
       .select_photon {
          tdc_emc_start     0
          tdc_emc_end       14
          angle_to_track    10.0
          energyThreshold_b 0.025
          energyThreshold_e 0.050
          nGam              ">=1"   # photons from pi0 (D*0 -> D0 pi0) or from D*0 -> D0 gamma
        }
       .pid(method: :probability) {
          prob_cut 0.001
          identify :kaon, against: [:pion, :proton]   # K+ and K-
          identify :pion, against: [:kaon, :proton]   # pi+ and pi-
          nkp  ">=1"
          nkm  ">=1"
          npip ">=1"
          npim ">=1"
        }
       # 1C Kalman fit : gamma gamma -> pi0 (chi2 < 200), at least one pi0
       .kalman_kinematic_fit([:gamma, :gamma]) {
          invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
          chi2_cut 200
          npi0     ">=1"
        }
       # Main fit : the pi- from D*- -> anti-D0 pi- is missed (partial reconstruction);
       # the D0 -> K- pi+ invariant mass is constrained to the nominal D0 mass.
       .kinematic_fit([:km, :pip, :kp, :pim, :pi0, :pip]) {
          nominal
          miss_track_of :pim                 # the D*- pi- is not reconstructed
          constrain_four_momentum            # 4C energy-momentum constraint
          invariant_mass_of(:km, :pip).constrain_to_nominal_mass_of(:D0)
          chi2_cut 200
        }

# The same algorithm is run at the second energy point: the CMS energy constant
# must be re-set to 4.416 GeV for that job (per-run beam energy handling).
alg_ch1.note(:cms_energy_per_point,
             "the algorithm is executed on both energy points; ECMS must be set to 4.416 GeV " \
             "for the 4.416 GeV dataset (only the 4.395 GeV value is encoded in the constant)")
       # D*+ veto : mass recoiling against (D0 pi+) must exceed 2.03 GeV/c2
       .note(:background_veto,
             "MQ(D0 pi+) > 2.03 GeV/c2 required to reject D*+ -> D0 pi+ background; " \
             "recoil mass built from the reconstructed D0 and the bachelor pi+")
       # Remaining mass windows are applied on kinematic-fit-corrected quantities
       # in the ROOT analysis (D0 1.850-1.880, pi0 0.120-0.145 GeV/c2, D*0 and
       # recoil-mass windows).
       .with_decay_card(decay_card_ch1)
       .apply(sel_ch1)

alg_ch1.execute_on(all_inputs + exMCs_ch1)

# ============================== Channel 2 ====================================
alg_name_ch2 = "Tcc4020pipiJpsi"
alg_ch2 = Algorithm.new(alg_name_ch2)
alg_ch2.set_header(["#{alg_name_ch2}Alg/#{alg_name_ch2}.h"])
       .set_constant({"ECMS" => [:double, 4.395]})

# Both lepton modes (e+e- and mu+mu-) share the same final-state topology and the
# same selection, so a single algorithm is used for the two exclusive MC samples.
sel_ch2 = Selection.new
sel_ch2.select_track {
          cos_theta 0.93
          Vr        1.0
          Vz        10.0
          nChrp     "==2"      # one pi+ and one l+
          nChrn     "==2"      # one pi- and one l-
          nNet      "==0"
        }
       .pid(method: :probability) {
          # tracks with p > 1.06 GeV/c are leptons; among them, EMC energy > 1.1 GeV
          # -> electron, otherwise -> muon (EMC energy < 0.35 GeV is the muon side
          # of the description; the default muon branch is used here).
          identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.06,
                                         treat_as_electron_if_energy_above: 1.1
          identify :pion, against: [:kaon]     # pi+ and pi-
          npip "==1"
          npim "==1"
          nlp  "==1"
          nlm  "==1"
        }
       # 4C fit with the l+l- invariant mass constrained to the J/psi mass
       .kinematic_fit([:pip, :pim, :lp, :lm]) {
          nominal
          constrain_four_momentum
          invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)
          chi2_cut 40
        }

alg_ch2.note(:cms_energy_per_point,
             "the algorithm is executed on both energy points; ECMS must be set to 4.416 GeV " \
             "for the 4.416 GeV dataset")
       # Muon quality : at least five MUC hit layers for muon candidates
       .note(:muon_muc_hit_layers,
             "muon candidates are required to have at least 5 hit layers in the muon counter (MUC); " \
             "applied on top of the probability-method PID from the BOSS algorithm")
       # Angular vetoes against mis-reconstructed J/psi / leptons and the J/psi mass window
       # (3.090-3.105 GeV/c2) are applied on kinematic-fit-corrected quantities in ROOT.
       .note(:background_veto,
             "cos(theta)(pi+pi-) < 0.98 and cos(theta)(pi e) < 0.98 required to reject " \
             "backgrounds with mis-assigned leptons")
       .with_decay_card(decay_card_ch2_ee)
       .apply(sel_ch2)

alg_ch2.execute_on(all_inputs + exMCs_ch2_ee + exMCs_ch2_mumu)

# ============================== Channel 3 ====================================
alg_name_ch3 = "Tcc4020pipiHcEtaC"
alg_ch3 = Algorithm.new(alg_name_ch3)
alg_ch3.set_header(["#{alg_name_ch3}Alg/#{alg_name_ch3}.h"])
       .set_constant({"ECMS" => [:double, 4.395]})

sel_ch3 = Selection.new
sel_ch3.select_track {
          cos_theta 0.93
          Vr        1.0
          Vz        10.0
          nChrp     ">=2"      # >= 1 pi+ from the recoil pair and >= 1 from eta_c
          nChrn     ">=2"      # >= 1 pi- from the recoil pair and >= 1 from eta_c
          nNet      "==0"
          nTot      ">=4"      # 4, 6 or 8 charged tracks depending on the eta_c mode
        }
       .select_photon {
          tdc_emc_start     0
          tdc_emc_end       14
          angle_to_track    10.0
          energyThreshold_b 0.025
          energyThreshold_e 0.050
          nGam              ">=1"    # photon from h_c -> gamma eta_c
        }
       .pid(method: :probability) {
          prob_cut 0.001
          identify :pion,  against: [:kaon, :proton]
          identify :kaon,  against: [:pion, :proton]
          identify :proton, against: [:kaon, :pion]
          npip ">=1"
          npim ">=1"
        }
       # 4C fit of pi+ pi- pi+ pi- gamma
       .kinematic_fit([:pip, :pim, :pip, :pim, :gamma]) {
          nominal
          constrain_four_momentum
          chi2_cut 35
        }

alg_ch3.note(:cms_energy_per_point,
             "the algorithm is executed on both energy points; ECMS must be set to 4.416 GeV " \
             "for the 4.416 GeV dataset")
       # The 16 hadronic eta_c modes are grouped by charged multiplicity (4, 6, 8 tracks);
       # the K_S0, pi0/eta, h_c and eta_c mass windows are applied on kinematic-fit-corrected
       # quantities in the ROOT analysis.
       .note(:eta_c_mode_grouping,
             "the 16 hadronic eta_c decay modes are reconstructed in three charged-multiplicity " \
             "groups (4, 6 and 8 tracks); K_S0 -> pi+pi- and pi0/eta -> gamma gamma candidates are " \
             "reconstructed before the 4C fit")
       .with_decay_card(decay_card_ch3)
       .apply(sel_ch3)

alg_ch3.execute_on(all_inputs + exMCs_ch3)