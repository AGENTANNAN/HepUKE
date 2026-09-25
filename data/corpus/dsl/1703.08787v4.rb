# =============================================================================
# BESIII arXiv:1703.08787v4
# Measurement of e+e- -> pi+ pi- psi(3686) from 4.008 to 4.600 GeV and
# observation of a charged structure in the pi^{+-} psi(3686) mass spectrum.
#
# 16 c.m. energy points from 4.008 to 4.600 GeV, total luminosity 5.1 fb^-1.
# The psi(3686) is reconstructed in two decay modes, J/psi -> l+ l- (l = e/mu)
# in both:
#   mode I  : psi(3686) -> pi+ pi- J/psi   (6-track 4C/5C and 5-track 1C/2C
#                                           topologies -> two algorithms)
#   mode II : psi(3686) -> neutrals + J/psi (4-track topology, no kinematic
#                                           fit; psi(3686) from the pi+pi-
#                                           recoil mass -> partial_rec)
# =============================================================================

### Dataset description ###
# The 16 c.m. energy points of Table I (the sample name carries the BESIII
# energy label, which can differ slightly from the measured c.m. energy).
energy_labels = %w[703_4009 703_4090 703_4190 703_4210 703_4220 703_4230
                   703_4245 703_4260 703_4310 703_4360 703_4390 703_4420
                   703_4470 703_4530 703_4575 703_4600]
data_points = energy_labels.map { |n| DatasetManager.real_data.find(n) }

# Generic inclusive MC (EvtGen with PDG branching fractions + Lundcharm for the
# unknown modes) is generated with luminosity equivalent to data only at two
# representative energy points.
incMC_4258 = DatasetManager.inclusive_mc.find("703_4260")   # sqrt(s) = 4.258 GeV
incMC_4358 = DatasetManager.inclusive_mc.find("703_4360")   # sqrt(s) = 4.358 GeV

# ---------------------------------------------------------------------------
# Decay cards.  The process is a continuum measurement of a Born cross section
# over the 16 scan points, with the ISR and vacuum-polarisation correction
# factors (1 + delta^r) and (1 + delta^v) taken from the generator, so the
# signal is generated with the ConExc model (mode 91 = psi(2S) pi+ pi-, the
# channel e+e- -> pi+pi- psi(3686)) instead of the KKMC + psi(4260) default.
# 'Particle vpho' is omitted so that every scan point receives its own sqrt(s).
# ---------------------------------------------------------------------------

# Mode I: e+e- -> pi+ pi- psi(3686), psi(3686) -> pi+ pi- J/psi, J/psi -> l+ l-
decay_card_modeI = <<~DECAYCARD
    Decay vpho
    1 ConExc 91;
    Enddecay

    Decay vhdr
    1 pi+ pi- psi(2S) PHSP;
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

# Mode II: e+e- -> pi+ pi- psi(3686), psi(3686) -> neutrals + J/psi, where
# 'neutrals' stands for pi0 pi0, pi0, eta and gamma gamma (all four allowed
# sub-modes are generated); J/psi -> l+ l-.
decay_card_modeII = <<~DECAYCARD
    Decay vpho
    1 ConExc 91;
    Enddecay

    Decay vhdr
    1 pi+ pi- psi(2S) PHSP;
    Enddecay

    Decay psi(2S)
    0.250 pi0 pi0 J/psi PHSP;
    0.250 pi0 J/psi PHSP;
    0.250 eta J/psi PHSP;
    0.250 gamma gamma J/psi PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
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

# Dominant background of mode II: the direct (ISR) process e+e- -> pi+ pi- J/psi,
# which has the same four-charged-track topology as the mode II signal and is
# removed by the |M^{corr} - M(psi(3686))| > 8 MeV/c^2 veto.
decay_card_bkg_pipiJpsi = <<~DECAYCARD
    Decay psi(4260)
    1.000 pi+ pi- J/psi PHSP;
    Enddecay

    Decay J/psi
    0.500 e+ e- PHOTOS VLL;
    0.500 mu+ mu- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# ---------------------------------------------------------------------------
# Exclusive MC samples — generated at every c.m. energy point
# ---------------------------------------------------------------------------
exMCs_modeI = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "ee_to_pipi_psi3686_2pipijpsi"     # per-point suffix
  config.events        = 100_000
  config.decay_card    = decay_card_modeI
  config.cross_section = :default                           # required; inert for ConExc
end

exMCs_modeII = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "ee_to_pipi_psi3686_neutralsjpsi"  # per-point suffix
  config.events        = 100_000
  config.decay_card    = decay_card_modeII
  config.cross_section = :default                           # required; inert for ConExc
end

exMCs_bkg_pipiJpsi = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "ee_to_pipi_jpsi_isr_bkg"          # per-point suffix
  config.events        = 100_000
  config.decay_card    = decay_card_bkg_pipiJpsi
  config.cross_section = :default
end

# ===========================================================================
# Algorithm I : mode I, e+e- -> pi+ pi- pi+ pi- l+ l-  (six charged tracks)
# ===========================================================================
alg_name_6trk = "PiPiJpsiLL6Trk"
alg_6trk = Algorithm.new(alg_name_6trk)
alg_6trk.set_header(["#{alg_name_6trk}Alg/#{alg_name_6trk}.h"])
        .set_constant({"ECMS" => [:double, 4.260]})

sel_6trk = Selection.new
  .select_track {
      cos_theta 0.93      # |cos(theta)| < 0.93 in the MDC
      Vz        10.0      # |Vz| < 10 cm along the beam
      Vr        1.0       # Vr < 1 cm in the plane perpendicular to the beam
      nChrp     "==3"     # three positive tracks (pi+ pi+ l+)
      nChrn     "==3"     # three negative tracks (pi- pi- l-)
      nNet      "==0"     # zero net charge for the six-track topology
  }
  # Pions and leptons are kinematically well separated: p > 1.0 GeV/c tracks are
  # treated as leptons and the remaining tracks as pions.
  .pid(method: :probability) {
      prob_cut 0.001
      identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                     treat_as_electron_if_energy_above: 0.6
      identify :pion, against: [:kaon, :proton]
      npip "==2"
      npim "==2"
      nlp  "==1"
      nlm  "==1"
  }
  # 4C kinematic fit under the e+e- -> pi+pi-pi+pi-l+l- hypothesis (chi2 < 60 in
  # the paper), with the M(l+l-) window; the track combination is then re-fitted
  # with the additional J/psi mass constraint (5C), which provides the corrected
  # four-momenta used to build the psi(3686) candidate.
  .kinematic_fit([:pip, :pim, :pip, :pim, :lp, :lm]) {
      constrain_four_momentum
      invariant_mass_of(:lp, :lm).within(3.05, 3.15)   # 3.05 < M(l+l-) < 3.15 GeV/c^2
      chi2_cut 60
  }
  .kinematic_fit([:pip, :pim, :pip, :pim, :lp, :lm]) {
      nominal
      constrain_four_momentum
      invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)   # 5C
  }

alg_6trk
  .note(:jpsi_mass_window,
        "the J/psi window 3.05 < M(l+l-) < 3.15 GeV/c^2 is applied inside the 4C fit " \
        "(invariant_mass_of(:lp,:lm).within(3.05,3.15)); in the paper the window is evaluated " \
        "with the fit-corrected momenta, which is equivalent at this precision.")
  .note(:best_candidate,
        "the psi(3686) candidate is the pi+pi-J/psi combination (four per event) whose mass, " \
        "calculated with the 5C-corrected momenta, is closest to the nominal psi(3686) mass; " \
        "the best-candidate choice is a post-fit selection applied in the ROOT analysis.")
  .note(:muc_hits,
        "J/psi -> mu+mu- candidates must have at least 5 hits in the muon counter (MUC); this is " \
        "not part of the DSL lepton PID and is applied in the ROOT analysis.")
  .note(:lepton_separation,
        "electrons are separated from muons by the EMC energy: E/pc > 0.7 for electrons and " \
        "E < 0.45 GeV for muons (the threshold-based lepton path of " \
        "identify_high_momentum_leptons is used here); the E/p and E requirements themselves are " \
        "applied in the ROOT analysis.")
  .note(:kinematic_fit_chi2,
        "the paper requires chi2_4C < 60 (chi2_1C < 15); the loose default chi2 cut is applied in " \
        "BOSS and the chi2 of the 5C (2C) fit is stored so that any tighter requirement is applied " \
        "in the ROOT analysis.")
  .note(:born_cross_section,
        "sigma^B_i = N^obs_i / (L_int (1+delta^r) (1+delta^v) Br_i eps_i); the ISR correction " \
        "(1+delta^r) and the vacuum-polarisation factor (1+delta^v) are taken from the ConExc " \
        "generator log, and Br = 4.11% for mode I.")
  .note(:conexc_mode_range,
        "ConExc mode 91 (psi(2S) pi+pi-) is tabulated for sqrt(s) = 4.127-5.480 GeV; the three " \
        "lowest scan points (4.008, 4.085, 4.189 GeV) lie below that range and the ISR line shape " \
        "there should be cross-checked against the measured cross section used iteratively.")
  .with_decay_card(decay_card_modeI)
  .apply(sel_6trk)

# ===========================================================================
# Algorithm II : mode I, e+e- -> pi+ pi- pi+ pi- l+ l-  (five charged tracks,
#                one pion undetected -> 1C / 2C kinematic fits)
# ===========================================================================
alg_name_5trk = "PiPiJpsiLL5Trk"
alg_5trk = Algorithm.new(alg_name_5trk)
alg_5trk.set_header(["#{alg_name_5trk}Alg/#{alg_name_5trk}.h"])
        .set_constant({"ECMS" => [:double, 4.260]})

sel_5trk = Selection.new
  .select_track {
      cos_theta 0.93
      Vz        10.0
      Vr        1.0
      nTot      "==5"     # five charged tracks with +-1 net charge
      nChrp     ">=2"     # (3 pi+ + l+ - missing pi-) or (pi+ + l+ + 2 pi- ... )
      nChrn     ">=2"
  }
  .pid(method: :probability) {
      prob_cut 0.001
      identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                     treat_as_electron_if_energy_above: 0.6
      identify :pion, against: [:kaon, :proton]
      npip ">=1"
      npim ">=1"
      nlp  "==1"
      nlm  "==1"
  }
  # 1C kinematic fit with one missing pion, then the 2C fit that also constrains
  # M(l+l-) to the nominal J/psi mass.
  .kinematic_fit([:pip, :pim, :pip, :pim, :lp, :lm]) {
      miss_track_of :pim                                # one pion undetected (1C)
      constrain_four_momentum
      invariant_mass_of(:lp, :lm).within(3.05, 3.15)
      chi2_cut 15
  }
  .kinematic_fit([:pip, :pim, :pip, :pim, :lp, :lm]) {
      nominal
      miss_track_of :pim
      constrain_four_momentum
      invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)   # 2C
  }

alg_5trk
  .note(:charge_conjugate_topology,
        "the five-track sample contains both +1 and -1 net charge configurations, i.e. the missing " \
        "pion may be a pi- (3 pi+ + l+ measured) or a pi+ (2 pi- + l- measured). The generated " \
        "selection assumes the missing particle to be a pi- in the 1C/2C fits; the charge-conjugate " \
        "configuration is covered by the same algorithm through the nChrp/nChrn >= 2 requirements " \
        "and is separated in the ROOT analysis.")
  .note(:jpsi_mass_window,
        "the J/psi window 3.05 < M(l+l-) < 3.15 GeV/c^2 is applied inside the 1C fit.")
  .note(:best_candidate,
        "at most one psi(3686) candidate per event is retained, chosen as the pi+pi-J/psi " \
        "combination whose mass after the 2C fit is closest to the nominal psi(3686) mass " \
        "(ROOT level).")
  .note(:muc_hits,
        "J/psi -> mu+mu- candidates require at least 5 muon-counter hits (ROOT level).")
  .note(:lepton_separation,
        "electron/muon separation by EMC energy (E/pc > 0.7 for electrons, E < 0.45 GeV for muons) " \
        "is applied in the ROOT analysis.")
  .note(:kinematic_fit_chi2,
        "the paper requires chi2_1C < 15; the loose cut is applied in BOSS and the 2C chi2 is stored " \
        "for a possible tighter requirement in ROOT.")
  .note(:born_cross_section,
        "the mode I Born cross section uses (1+delta^r) and (1+delta^v) from the ConExc generator " \
        "log and Br = 4.11%; the five- and six-track samples are combined in the extraction.")
  .with_decay_card(decay_card_modeI)
  .apply(sel_5trk)

# ===========================================================================
# Algorithm III : mode II, e+e- -> pi+ pi- (neutrals) J/psi -> pi+ pi- l+ l- + N gamma
# ===========================================================================
alg_name_neut = "NeutralsJpsiLL"
alg_neut = Algorithm.new(alg_name_neut)
alg_neut.set_header(["#{alg_name_neut}Alg/#{alg_name_neut}.h"])
        .set_constant({"ECMS" => [:double, 4.260]})

sel_neut = Selection.new
  .select_track {
      cos_theta 0.93      # |cos(theta)| < 0.93
      Vz        10.0
      Vr        1.0
      nChrp     "==2"     # two positive tracks (pi+ l+)
      nChrn     "==2"     # two negative tracks (pi- l-)
      nNet      "==0"
  }
  .select_photon {
      tdc_emc_start     0     # EMC timing window (0, 700) ns
      tdc_emc_end       14
      angle_to_track    10.0  # isolated from all charged tracks by > 10 degrees
      energyThreshold_b 0.025 # E > 25 MeV in the barrel (|cos(theta)| < 0.8)
      energyThreshold_e 0.050 # E > 50 MeV in the endcap (0.86 < |cos(theta)| < 0.92)
      nGam              ">=2" # at least two photon candidates (the neutral system)
  }
  .pid(method: :probability) {
      prob_cut 0.001
      identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                     treat_as_electron_if_energy_above: 0.6
      identify :pion, against: [:kaon, :proton]
      npip "==1"
      npim "==1"
      nlp  "==1"
      nlm  "==1"
  }
  # No kinematic fit is applied in mode II: the psi(3686) signal is identified
  # through the mass recoiling against the pi+pi- system.  Reconstructing only
  # the two pions (recID 1 = pi+, recID 2 = pi-) makes the DSL recoil mass
  # P4_cms - p(pi+pi-) identical to the paper's M^{recoil}(pi+pi-).
  .partial_rec([1, 2]) {
      require_recoil_mass 3.63, 3.75   # loose window covering signal and sidebands
  }

alg_neut
  .note(:recoil_mass_window,
        "the psi(3686) signal is extracted from the M^{recoil}(pi+pi-) spectrum: " \
        "3.68 < M^{recoil} < 3.70 GeV/c^2 for the signal and 3.63-3.65 / 3.73-3.75 GeV/c^2 for the " \
        "sidebands. A loose window 3.63-3.75 GeV/c^2 is applied in BOSS " \
        "(require_recoil_mass); the signal/sideband separation is done in the ROOT analysis.")
  .note(:jpsi_mass_window,
        "the J/psi signal is selected by 3.05 < M(l+l-) < 3.15 GeV/c^2, evaluated with the measured " \
        "(uncorrected) lepton momenta since no kinematic fit is applied in mode II; handled in the " \
        "ROOT analysis.")
  .note(:background_veto,
        "|M^{corr}_{psi(3686)} - M(psi(3686))| > 8 MeV/c^2, with " \
        "M^{corr}_{psi(3686)} = M(pi+pi-l+l-) - M(l+l-) + M(J/psi), vetoes the e+e- -> pi+pi- J/psi " \
        "background; |M(gamma gamma pi+pi-) - M(eta)| > 50 MeV/c^2 (the two photons of largest " \
        "energy) vetoes e+e- -> eta J/psi (eta -> pi+pi-pi0). Both vetoes are applied in ROOT.")
  .note(:angular_cut,
        "cos(theta_{pi+pi-}) < 0.9 removes the radiative Bhabha and dimuon background in which a " \
        "photon converts into an e+e- pair misidentified as a pi+pi- pair; applied in ROOT.")
  .note(:muc_hits,
        "J/psi -> mu+mu- candidates require at least 5 muon-counter hits (ROOT level).")
  .note(:lepton_separation,
        "electron/muon separation by EMC energy (E/pc > 0.7 for electrons, E < 0.45 GeV for muons) " \
        "is applied in the ROOT analysis.")
  .note(:born_cross_section,
        "the mode II Born cross section uses (1+delta^r) and (1+delta^v) from the ConExc generator " \
        "log and Br = 2.95%.")
  .note(:conexc_mode_range,
        "ConExc mode 91 (psi(2S) pi+pi-) is tabulated for sqrt(s) = 4.127-5.480 GeV; the three " \
        "lowest scan points lie below that range and the ISR line shape should be cross-checked " \
        "against the measured cross section used iteratively.")
  .with_decay_card(decay_card_modeII)
  .apply(sel_neut)

### Execution ###
root_files_6trk = alg_6trk.execute_on(data_points + [incMC_4258, incMC_4358] + exMCs_modeI)
root_files_5trk = alg_5trk.execute_on(data_points + [incMC_4258, incMC_4358] + exMCs_modeI)
root_files_neut = alg_neut.execute_on(data_points + [incMC_4258, incMC_4358] +
                                      exMCs_modeII + exMCs_bkg_pipiJpsi)
