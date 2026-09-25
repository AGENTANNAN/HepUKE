# =============================================================================
# e+e- -> omega chi_c1,2  (chi_c1,2 -> gamma J/psi, J/psi -> e+e- / mu+mu-)
# e+e- -> omega chi_c0    (chi_c0 -> pi+pi-  or  K+K-),  omega -> pi+pi- pi0
# at sqrt(s) = 4.416, 4.467, 4.527, 4.574 and 4.599 GeV
# BOSS part only: decay cards + exclusive MC + event selection up to the 5C fit
# =============================================================================

### Dataset description ###
data_points = [
  DatasetManager.real_data.find("703_4420"),   # 4.416 GeV
  DatasetManager.real_data.find("703_4470"),   # 4.467 GeV
  DatasetManager.real_data.find("703_4530"),   # 4.527 GeV
  DatasetManager.real_data.find("703_4575"),   # 4.574 GeV
  DatasetManager.real_data.find("703_4600")    # 4.599 GeV
]
incMC_points = [
  DatasetManager.inclusive_mc.find("703_4420"),  # inclusive MC at 4.416 GeV
  DatasetManager.inclusive_mc.find("703_4470"),  # inclusive MC at 4.467 GeV
  DatasetManager.inclusive_mc.find("703_4530"),  # inclusive MC at 4.527 GeV
  DatasetManager.inclusive_mc.find("703_4575"),  # inclusive MC at 4.574 GeV
  DatasetManager.inclusive_mc.find("703_4600")   # inclusive MC at 4.599 GeV
]

# --------------------------- Decay cards (EvtGen) -----------------------------
# --- chi_c1 -> gamma J/psi, J/psi -> e+e- ---
decay_card_c1_ee = <<~DECAYCARD
  Decay psi(4260)
  1.000 omega chi_c1 PHSP;
  Enddecay

  Decay chi_c1
  1.000 gamma J/psi PHSP;
  Enddecay

  Decay J/psi
  1.000 e+ e- PHOTOS VLL;
  Enddecay

  Decay omega
  1.000 pi+ pi- pi0 OMEGA_DALITZ;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# --- chi_c1 -> gamma J/psi, J/psi -> mu+mu- ---
decay_card_c1_mumu = <<~DECAYCARD
  Decay psi(4260)
  1.000 omega chi_c1 PHSP;
  Enddecay

  Decay chi_c1
  1.000 gamma J/psi PHSP;
  Enddecay

  Decay J/psi
  1.000 mu+ mu- PHOTOS VLL;
  Enddecay

  Decay omega
  1.000 pi+ pi- pi0 OMEGA_DALITZ;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# --- chi_c2 -> gamma J/psi, J/psi -> e+e- ---
decay_card_c2_ee = <<~DECAYCARD
  Decay psi(4260)
  1.000 omega chi_c2 PHSP;
  Enddecay

  Decay chi_c2
  1.000 gamma J/psi PHSP;
  Enddecay

  Decay J/psi
  1.000 e+ e- PHOTOS VLL;
  Enddecay

  Decay omega
  1.000 pi+ pi- pi0 OMEGA_DALITZ;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# --- chi_c0 -> pi+pi- ---
decay_card_c0_pipi = <<~DECAYCARD
  Decay psi(4260)
  1.000 omega chi_c0 PHSP;
  Enddecay

  Decay chi_c0
  1.000 pi+ pi- PHSP;
  Enddecay

  Decay omega
  1.000 pi+ pi- pi0 OMEGA_DALITZ;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# --- chi_c0 -> K+K- ---
decay_card_c0_kk = <<~DECAYCARD
  Decay psi(4260)
  1.000 omega chi_c0 PHSP;
  Enddecay

  Decay chi_c0
  1.000 K+ K- PHSP;
  Enddecay

  Decay omega
  1.000 pi+ pi- pi0 OMEGA_DALITZ;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# ------------------------- Exclusive signal MC (100k) -------------------------
# One 100k-event signal MC per mode per energy point (same card/cross section
# across the five energy points -> create_exclusive_mc_for).
exMC_c1_ee = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_omegachic1_jpsiee"
  config.events        = 100_000
  config.decay_card    = decay_card_c1_ee
  config.cross_section = :default
end

exMC_c1_mumu = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_omegachic1_jpsimumu"
  config.events        = 100_000
  config.decay_card    = decay_card_c1_mumu
  config.cross_section = :default
end

exMC_c2_ee = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_omegachic2_jpsiee"
  config.events        = 100_000
  config.decay_card    = decay_card_c2_ee
  config.cross_section = :default
end

exMC_c0_pipi = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_omegachic0_pipi"
  config.events        = 100_000
  config.decay_card    = decay_card_c0_pipi
  config.cross_section = :default
end

exMC_c0_kk = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_omegachic0_kk"
  config.events        = 100_000
  config.decay_card    = decay_card_c0_kk
  config.cross_section = :default
end

# create_exclusive_mc_for returns an Array<ExclusiveMC> -> save element by element
(exMC_c1_ee + exMC_c1_mumu + exMC_c2_ee + exMC_c0_pipi + exMC_c0_kk).each do |m|
  m.save_to_config(format: :yaml, file_path: 'temp_for_test')
end

# ============================ Event selection (BOSS) ==========================

# -----------------------------------------------------------------------------
# Algorithm 1 : the "chi_c1,2 mode"  (chi_c1 -> gamma J/psi and chi_c2 -> gamma J/psi,
#               J/psi -> e+e- or mu+mu-) with omega -> pi+pi- pi0, pi0 -> gamma gamma
# The three lepton sub-channels share an identical final state
# (gamma l+ l- pi+ pi- pi0) and identical selection, so a single Algorithm /
# Selection chain serves them; the header is generated from the chi_c1 card.
# -----------------------------------------------------------------------------
alg_name_lep = "OmegaChiC1C2Lep"
alg_lep = Algorithm.new(alg_name_lep)
alg_lep.set_header(["#{alg_name_lep}Alg/#{alg_name_lep}.h"])
       .set_constant({"ECMS" => [:double, 4.527]})  # scan mid-point; per-run beam energy is set per dataset at execution
       .set_alias({"std::vector<double>" => "Vdouble"})

sel_lep = Selection.new
sel_lep.select_track {                 # four good charged tracks, 2 positive + 2 negative
          cos_theta 0.93               # |cos(theta)| < 0.93
          Vz        10.0               # |Vz| < 10 cm
          Vr        1.0                # Vr < 1 cm
          nChrp     "==2"              # exactly 2 positive tracks
          nChrn     "==2"              # exactly 2 negative tracks
          nNet      "==0"              # net charge zero
       }
       .select_photon {                # photon selection
          tdc_emc_start     0
          tdc_emc_end       14
          angle_to_track    10.0       # at least 10 degrees from any charged track
          energyThreshold_b 0.025      # > 25 MeV in the barrel
          energyThreshold_e 0.050      # > 50 MeV in the endcap
          nGam              ">=3"      # 2 photons from pi0 + 1 prompt photon from chi_c1,2 -> gamma J/psi
       }
       .pid(method: :probability) {    # PID, probability method
          prob_cut 0.001               # probability > 0.001
          identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                         treat_as_electron_if_energy_above: 0.6  # p > 1.0 GeV -> lepton; EMC E > 0.6 GeV -> e, else mu
          identify :pion, against: [:kaon, :proton]   # pi+/pi- vs K, p
          nlp  "==1"                   # exactly one lepton of each charge
          nlm  "==1"
          npip "==1"
          npim "==1"
       }
       # 5C fit: 4-momentum conservation (4C) + pi0 mass constraint on the photon
       # pair with the smallest chi2 (1C)
       .kinematic_fit([:gamma, :gamma, :gamma, :lp, :lm, :pip, :pim]) {
          nominal
          constrain_four_momentum
          invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
          chi2_cut 200                 # BOSS loose cut; published analysis applies chi2_5C < 60 in ROOT
       }

# -----------------------------------------------------------------------------
# Algorithm 2 : chi_c0 -> pi+pi-, omega -> pi+pi- pi0
# -----------------------------------------------------------------------------
alg_name_pipi = "OmegaChiC0PiPi"
alg_pipi = Algorithm.new(alg_name_pipi)
alg_pipi.set_header(["#{alg_name_pipi}Alg/#{alg_name_pipi}.h"])
        .set_constant({"ECMS" => [:double, 4.527]})
        .set_alias({"std::vector<double>" => "Vdouble"})

sel_pipi = Selection.new
sel_pipi.select_track {
           cos_theta 0.93
           Vz        10.0
           Vr        1.0
           nChrp     "==2"
           nChrn     "==2"
           nNet      "==0"
         }
        .select_photon {
           tdc_emc_start     0
           tdc_emc_end       14
           angle_to_track    10.0
           energyThreshold_b 0.025
           energyThreshold_e 0.050
           nGam              ">=2"     # both photons come from pi0 (no prompt photon in chi_c0 modes)
         }
        .pid(method: :probability) {
           prob_cut 0.001
           identify :pion, against: [:kaon, :proton]   # all four tracks are pions
           npip "==2"
           npim "==2"
         }
        .kinematic_fit([:gamma, :gamma, :pip, :pip, :pim, :pim]) {
           nominal
           constrain_four_momentum
           invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
           chi2_cut 200
         }

# -----------------------------------------------------------------------------
# Algorithm 3 : chi_c0 -> K+K-, omega -> pi+pi- pi0
# -----------------------------------------------------------------------------
alg_name_kk = "OmegaChiC0KK"
alg_kk = Algorithm.new(alg_name_kk)
alg_kk.set_header(["#{alg_name_kk}Alg/#{alg_name_kk}.h"])
      .set_constant({"ECMS" => [:double, 4.527]})
      .set_alias({"std::vector<double>" => "Vdouble"})

sel_kk = Selection.new
sel_kk.select_track {
         cos_theta 0.93
         Vz        10.0
         Vr        1.0
         nChrp     "==2"
         nChrn     "==2"
         nNet      "==0"
       }
      .select_photon {
         tdc_emc_start     0
         tdc_emc_end       14
         angle_to_track    10.0
         energyThreshold_b 0.025
         energyThreshold_e 0.050
         nGam              ">=2"
       }
      .pid(method: :probability) {
         prob_cut 0.001
         identify :pion, against: [:kaon, :proton]   # pi+/pi- from omega
         identify :kaon, against: [:pion, :proton]   # K+/K- from chi_c0
         nkp  "==1"
         nkm  "==1"
         npip "==1"
         npim "==1"
       }
      .kinematic_fit([:gamma, :gamma, :kp, :km, :pip, :pim]) {
         nominal
         constrain_four_momentum
         invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
         chi2_cut 200
       }

# --------------------------- Render and execute -------------------------------
alg_lep.with_decay_card(decay_card_c1_ee).apply(sel_lep)
alg_pipi.with_decay_card(decay_card_c0_pipi).apply(sel_pipi)
alg_kk.with_decay_card(decay_card_c0_kk).apply(sel_kk)

# All signal windows (omega, J/psi, chi_c1, chi_c2, chi_c0) and the tight
# chi2_5C < 60 cut are applied later in the ROOT analysis, hence not included here.
root_files_lep  = alg_lep.execute_on(data_points + incMC_points +
                                     exMC_c1_ee + exMC_c1_mumu + exMC_c2_ee)
root_files_pipi = alg_pipi.execute_on(data_points + incMC_points + exMC_c0_pipi)
root_files_kk   = alg_kk.execute_on(data_points + incMC_points + exMC_c0_kk)