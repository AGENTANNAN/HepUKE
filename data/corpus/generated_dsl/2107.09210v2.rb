# =============================================================================
# BOSS side: dataset preparation + event selection
# e+e- -> pi+ pi- psi(3686) over the 4.0076 - 4.6984 GeV scan
#   Mode I  (charged): psi(3686) -> pi+ pi- J/psi , J/psi -> e+ e-
#   Mode II (neutral): psi(3686) -> neutrals  J/psi , J/psi -> e+ e-
# =============================================================================

### ------------------------- Dataset preparation ------------------------- ###
# Real data + inclusive MC at every scan point inside [4.0076, 4.6984] GeV
scan_points = %w[
  703_4009 703_4180 703_4190 703_4200 703_4210 703_4220 703_4230
  703_4237 703_4246 703_4260 703_4270 703_4280 703_4360 703_4420
  703_4600 703_4610 703_4620 703_4640 703_4660 703_4680 703_4700
]
scan_data  = scan_points.map { |p| DatasetManager.real_data.find(p) }
scan_incMC = scan_points.map { |p| DatasetManager.inclusive_mc.find(p) }

# Mode I decay card: e+e- -> pi+ pi- psi(3686), psi(3686) -> pi+ pi- J/psi, J/psi -> e+ e-
decay_card_modeI = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ pi- psi(3686) PHSP;
  Enddecay

  Decay psi(3686)
  1.0000 pi+ pi- J/psi PHSP;
  Enddecay

  Decay J/psi
  1.0000 e+ e- PHOTOS VLL;
  Enddecay

  End
DECAYCARD

# Mode II decay card: e+e- -> pi+ pi- psi(3686), psi(3686) -> neutrals (eta -> gamma gamma) J/psi, J/psi -> e+ e-
decay_card_modeII = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi+ pi- psi(3686) PHSP;
  Enddecay

  Decay psi(3686)
  1.0000 eta J/psi PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay

  Decay J/psi
  1.0000 e+ e- PHOTOS VLL;
  Enddecay

  End
DECAYCARD

# Exclusive signal MC: 100k events per energy point, for each of the two modes
exMC_modeI = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name   = "exmc_pipipsi3686_pipijpsi_ee"
  config.events        = 100_000
  config.decay_card    = decay_card_modeI
  config.cross_section = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name   = "exmc_pipipsi3686_neutraljpsi_ee"
  config.events        = 100_000
  config.decay_card    = decay_card_modeII
  config.cross_section = :default
end

### ------------------------ Mode I  (charged) ------------------------ ###
alg_name_modeI = "PiPiPsi3686ModeI"
alg_modeI = Algorithm.new(alg_name_modeI)
alg_modeI.set_header(["#{alg_name_modeI}Alg/#{alg_name_modeI}.h"])
        # ECMS is the nominal scan energy; it is set to the actual sqrt(s)
        # of each data/MC point at job level.
        .set_constant({"ECMS" => [:double, 4.26]})
        .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeI = Selection.new
  # Charged-track quality cuts: |cos(theta)| < 0.93, |Vz| < 10 cm, Vr < 1 cm,
  # at least two positive and two negative tracks (5- or 6-track topology).
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=2"
    nChrn     ">=2"
    nTot      ">=5"
  }
  # PID: probability method; high-momentum tracks are leptons
  # (electron if EMC energy > 0.7 GeV, otherwise muon), the rest are pions.
  # At least one e+ and one e- are required (same-flavour lepton pair).
  .pid(method: :probability) {
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.7
    identify :pion, against: [:kaon, :proton]
    nlp ">=1"
    nlm ">=1"
  }
  # Nominal 4C fit to pi+pi-pi+pi-e+e- (six charged tracks), chi2 < 60
  .kinematic_fit([:pip, :pim, :pip, :pim, :lp, :lm]) {
    nominal
    constrain_four_momentum
    chi2_cut 60
  }
  # Five-track 1C missing-pion fit (one pion lost), chi2 < 15
  .kinematic_fit([:pip, :pim, :pip, :pim, :lp, :lm]) {
    miss_track_of(:pim)
    constrain_four_momentum
    chi2_cut 15
  }
  # J/psi mass window [3.05, 3.15] GeV/c^2 pre-selection, then the subsequent
  # 5C fit (4C + J/psi mass constraint) after the J/psi mass constraint
  .kinematic_fit([:pip, :pim, :pip, :pim, :lp, :lm]) {
    invariant_mass_of(:lp, :lm).within(3.05, 3.15)
    invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)
    constrain_four_momentum
    chi2_cut 60
  }

alg_modeI
  .note(:missing_pion_fit, "documented five-track 1C missing-pion fit (one pion
    undetected) is performed in addition to the nominal six-track 4C fit;
    chi2 < 15 for the 1C fit")
  .note(:pion_momentum_threshold, "tracks are treated as pions only below a
    momentum threshold that depends on the scan point: p < 0.65 GeV/c for
    sqrt(s) < 4.465 GeV and p < 0.80 GeV/c for sqrt(s) > 4.465 GeV; the
    threshold is set per energy point at job level (not expressible as a single
    static cut)")
  .note(:same_flavour_lepton_pair, "same-flavour opposite-charge lepton pairs
    (e+e- from J/psi) are required; enforced on the high-momentum lepton
    candidates")
  .note(:energy_scan_ecms, "the algorithm runs over the 4.0076-4.6984 GeV scan;
    the ECMS constant must be set to the measured sqrt(s) of each data/MC point")
  .with_decay_card(decay_card_modeI)
  .apply(sel_modeI)

root_files_modeI = alg_modeI.execute_on(scan_data + scan_incMC + exMC_modeI)

### ------------------------ Mode II  (neutral) ----------------------- ###
alg_name_modeII = "PiPiPsi3686ModeII"
alg_modeII = Algorithm.new(alg_name_modeII)
alg_modeII.set_header(["#{alg_name_modeII}Alg/#{alg_name_modeII}.h"])
         .set_constant({"ECMS" => [:double, 4.26]})
         .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeII = Selection.new
  # Four charged tracks (pi+ pi- e+ e-), net charge zero, with the same
  # quality cuts (|cos(theta)| < 0.93, |Vz| < 10 cm, Vr < 1 cm)
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"
  }
  # At least two photons: EMC timing 0-14, barrel/endcap thresholds
  # 25/50 MeV, opening angle to any charged track > 10 degrees
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end   14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 10.0
    nGam ">=2"
  }
  # Same PID strategy as Mode I
  .pid(method: :probability) {
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.7
    identify :pion, against: [:kaon, :proton]
    nlp ">=1"
    nlm ">=1"
  }
  # J/psi mass window [3.05, 3.15] GeV/c^2 (no kinematic fit in this mode)
  .invariant_mass_of(:lp, :lm).within(3.05, 3.15)
  # Background veto against the eta mass region in the gamma gamma pi+ pi- system
  .invariant_mass_of(:gamma, :gamma, :pip, :pim).out_of(0.498, 0.598)

alg_modeII
  .note(:background_veto, "events with cos(theta)(pi+pi-) < 0.9 are vetoed
    (the two pions are required to be nearly back-to-back)")
  .note(:background_veto_2, "events with |M_rec(pi+pi-e+e-)| > 63 MeV/c^2 are
    vetoed; recoil mass computed against the four charged tracks")
  .note(:background_veto_3, "events with |M_corr(psi(3686)) - M(psi(3686))| >
    8 MeV/c^2 are vetoed, M_corr being the psi(3686) mass corrected for the
    missing momentum carried by the neutral system")
  .note(:no_kinematic_fit, "Mode II applies only mass windows and vetoes; no
    kinematic fit is performed in this channel")
  .note(:pion_momentum_threshold, "tracks are treated as pions only below a
    momentum threshold that depends on the scan point: p < 0.65 GeV/c for
    sqrt(s) < 4.465 GeV and p < 0.80 GeV/c for sqrt(s) > 4.465 GeV")
  .note(:same_flavour_lepton_pair, "same-flavour opposite-charge lepton pairs
    (e+e- from J/psi) are required")
  .note(:energy_scan_ecms, "ECMS must be set to the measured sqrt(s) of each
    data/MC point of the 4.0076-4.6984 GeV scan")
  .with_decay_card(decay_card_modeII)
  .apply(sel_modeII)

root_files_modeII = alg_modeII.execute_on(scan_data + scan_incMC + exMC_modeII)