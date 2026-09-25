# =====================================================================
# e+e- -> pi+ pi- psi(3686)  @ 16 c.m. energy points (4.008 - 4.600 GeV)
# BOSS part: dataset preparation + event selection (up to and including
# the final kinematic fit).  Best-candidate / MUC / veto / angular and
# the final J/psi mass window are ROOT-level and therefore out of scope.
# =====================================================================

### ------------------------------------------------------------------ ###
### 1. Datasets                                                        ###
### ------------------------------------------------------------------ ###
# 16 c.m. energy points (BOSS 703); sample name = [BOSS]_[Ecms(MeV)]
energy_point_names = %w[
  703_4009 703_4180 703_4190 703_4200
  703_4210 703_4220 703_4230 703_4245
  703_4260 703_4270 703_4280 703_4310
  703_4360 703_4390 703_4420 703_4600
]
data_points = energy_point_names.map { |n| DatasetManager.real_data.find(n) }

# Inclusive MC generated at 4.258 GeV and 4.358 GeV
incMC_4258 = DatasetManager.inclusive_mc.find("703_4260")
incMC_4358 = DatasetManager.inclusive_mc.find("703_4360")

### ------------------------------------------------------------------ ###
### 2. Decay cards                                                     ###
###    (signal produced through ConExc mode 91: e+e- -> psi(2S) pi+ pi-)###
### ------------------------------------------------------------------ ###

# Mode I : psi(3686) -> pi+ pi- J/psi,  J/psi -> e+e- / mu+mu-
decay_card_modeI = <<~DECAYCARD
    Decay vpho
    1.000 ConExc 91;
    Enddecay

    Decay psi(2S)
    1.000 pi+ pi- J/psi PHSP;
    Enddecay

    Decay J/psi
    0.500 e+ e- PHOTOS VLL;
    0.500 mu+ mu- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# Mode II : psi(3686) -> pi0 pi0 J/psi, pi0 J/psi, eta J/psi, gamma gamma J/psi
#           J/psi -> l+ l-
decay_card_modeII = <<~DECAYCARD
    Decay vpho
    1.000 ConExc 91;
    Enddecay

    Decay psi(2S)
    0.250 pi0 pi0 J/psi PHSP;
    0.250 pi0 J/psi PHSP;
    0.250 eta J/psi PHSP;
    0.250 gamma gamma J/psi PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    Decay J/psi
    0.500 e+ e- PHOTOS VLL;
    0.500 mu+ mu- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# ISR background : e+e- -> pi+ pi- J/psi,  J/psi -> l+ l-
decay_card_isr = <<~DECAYCARD
    Decay psi(4260)
    1.000 pi+ pi- J/psi PHSP;
    Enddecay

    Decay J/psi
    0.500 e+ e- PHOTOS VLL;
    0.500 mu+ mu- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

### ------------------------------------------------------------------ ###
### 3. Exclusive MC (100k events, generated at every energy point)     ###
### ------------------------------------------------------------------ ###
# Same signal / background card run over the 16 distinct energy points
exMC_modeI = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_pipi_psi2S_modeI"
  config.events        = 100_000
  config.decay_card    = decay_card_modeI
  config.cross_section = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_pipi_psi2S_modeII"
  config.events        = 100_000
  config.decay_card    = decay_card_modeII
  config.cross_section = :default
end

exMC_isr = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_isr_pipiJpsi"
  config.events        = 100_000
  config.decay_card    = decay_card_isr
  config.cross_section = :default
end

# Persist the exclusive-MC configurations (each array element separately)
(exMC_modeI + exMC_modeII + exMC_isr).each do |m|
  m.save_to_config(format: :yaml, file_path: 'temp_for_test')
end

### ------------------------------------------------------------------ ###
### 4. Mode I algorithm  (six charged tracks)                          ###
### ------------------------------------------------------------------ ###
alg_name_I6 = "PipiPsi2S6Trk"
alg_modeI_6trk = Algorithm.new(alg_name_I6)
alg_modeI_6trk
  .set_header(["#{alg_name_I6}Alg/#{alg_name_I6}.h"])
  .set_constant({ "ECMS" => [:double, 4.260] })
  .set_alias({ "std::vector<double>" => "Vdouble" })
  .note(:energy_scan,
        "selection is applied at the 16 c.m. energy points 4.008-4.600 GeV; " \
        "the ECMS constant set here is a placeholder - the per-point beam " \
        "energy must be used for the 4C/5C constraints")

sel_modeI_6trk = Selection.new
sel_modeI_6trk
  # charged-track quality cuts (common to all modes)
  .select_track {
    cos_theta 0.93          # |cos(theta)| < 0.93
    Vz        10.0          # |Vz| < 10 cm
    Vr        1.0           # Vr < 1 cm
    nChrp     "==3"         # 3 positive tracks
    nChrn     "==3"         # 3 negative tracks
    nNet      "==0"         # net charge 0
  }
  .pid(method: :probability) {
    prob_cut 0.001
    # p > 1.0 GeV/c tracks treated as leptons (electron if EMC E > 0.6 GeV)
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    identify :pion, against: [:kaon]   # remaining tracks as pions with pi/K separation
    npip "==2"
    npim "==2"
    nlp  "==1"
    nlm  "==1"
  }
  # 4C fit to pi+ pi- pi+ pi- l+ l-, J/psi mass window, chi2 < 60
  .kinematic_fit([:pip, :pip, :pim, :pim, :lp, :lm]) {
    nominal
    invariant_mass_of(:lp, :lm).within(3.05, 3.15)
    constrain_four_momentum
    chi2_cut 60
  }
  # subsequent 5C fit, adding the J/psi nominal-mass constraint
  .kinematic_fit([:pip, :pip, :pim, :pim, :lp, :lm]) {
    invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)
    constrain_four_momentum
  }

alg_modeI_6trk.with_decay_card(decay_card_modeI).apply(sel_modeI_6trk)

### ------------------------------------------------------------------ ###
### 5. Mode I algorithm  (five charged tracks, one pion undetected)    ###
### ------------------------------------------------------------------ ###
alg_name_I5 = "PipiPsi2S5Trk"
alg_modeI_5trk = Algorithm.new(alg_name_I5)
alg_modeI_5trk
  .set_header(["#{alg_name_I5}Alg/#{alg_name_I5}.h"])
  .set_constant({ "ECMS" => [:double, 4.260] })
  .note(:energy_scan,
        "selection is applied at the 16 c.m. energy points 4.008-4.600 GeV; " \
        "the ECMS constant set here is a placeholder - the per-point beam " \
        "energy must be used for the 1C/2C constraints")

sel_modeI_5trk = Selection.new
sel_modeI_5trk
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=2"         # >= 2 positive tracks
    nChrn     ">=2"         # >= 2 negative tracks
    nTot      "==5"         # exactly five charged tracks
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    identify :pion, against: [:kaon]
    npip ">=1"
    npim ">=1"
    nlp  "==1"
    nlm  "==1"
  }
  # 1C fit with one undetected pion, J/psi mass window, chi2 < 15
  .kinematic_fit([:pip, :pip, :pim, :lp, :lm]) {
    nominal
    miss_track_of(:pim)
    invariant_mass_of(:lp, :lm).within(3.05, 3.15)
    chi2_cut 15
  }
  # subsequent 2C fit, adding the J/psi nominal-mass constraint
  .kinematic_fit([:pip, :pip, :pim, :lp, :lm]) {
    miss_track_of(:pim)
    invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)
  }

alg_modeI_5trk.with_decay_card(decay_card_modeI).apply(sel_modeI_5trk)

### ------------------------------------------------------------------ ###
### 6. Mode II algorithm  (2 pi + 2 l + >= 2 photons, no kinematic fit)###
###    psi(3686) selected through the pi+ pi- recoil mass              ###
### ------------------------------------------------------------------ ###
alg_name_II = "PipiPsi2SModeII"
alg_modeII = Algorithm.new(alg_name_II)
alg_modeII
  .set_header(["#{alg_name_II}Alg/#{alg_name_II}.h"])
  .set_constant({ "ECMS" => [:double, 4.260] })
  .note(:energy_scan,
        "selection is applied at the 16 c.m. energy points 4.008-4.600 GeV; " \
        "the ECMS constant set here is a placeholder")

sel_modeII = Selection.new
sel_modeII
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0          # EMC timing 0 - 14 (x 700 ns)
    tdc_emc_end       14
    angle_to_track    10.0       # > 10 deg isolation from tracks
    energyThreshold_b 0.025      # E > 25 MeV in barrel
    energyThreshold_e 0.050      # E > 50 MeV in endcap
    nGam ">=2"                   # at least two photons
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    identify :pion, against: [:kaon]
    npip "==1"
    npim "==1"
    nlp  "==1"
    nlm  "==1"
  }
  # No kinematic fit for mode II.  The psi(3686) is instead selected through
  # the pi+ pi- recoil mass: reconstruct the two primary pions and require the
  # recoil mass (e+e- minus the pi+pi- pair) to lie in 3.63 - 3.75 GeV/c^2.
  .partial_rec([1, 2]) {                       # recIDs of the primary pi+, pi-
    require_recoil_mass 3.63, 3.75
  }

alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)

### ------------------------------------------------------------------ ###
### 7. Execution                                                       ###
### ------------------------------------------------------------------ ###
# Real data (16 points) + inclusive MC (4.258, 4.358 GeV) + exclusive MC
datasets_modeI = data_points + [incMC_4258, incMC_4358] + exMC_modeI + exMC_isr
datasets_modeII = data_points + [incMC_4258, incMC_4358] + exMC_modeII + exMC_isr

root_files_modeI_6trk = alg_modeI_6trk.execute_on(datasets_modeI)
root_files_modeI_5trk = alg_modeI_5trk.execute_on(datasets_modeI)
root_files_modeII     = alg_modeII.execute_on(datasets_modeII)