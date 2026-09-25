# =============================================================================
# Inclusive prompt J/psi and psi(3686) production cross sections in e+e-
# annihilation at 49 energy points, sqrt(s) = 3.808 - 4.951 GeV
# (BESIII BOSS releases 703 / 705 / 706 / 707, ~22 fb^-1).
#
# BOSS-side scope ONLY: dataset preparation + event selection up to and
# including the final 4C kinematic fit. The post-fit yield extraction
# (M(ll) / M(pi+pi-ll) / M(gamma J/psi) window fits), the ISR gamma_ISR
# background and psi(3686)->J/psi X / chi_cJ -> gamma J/psi feed-down
# subtractions, and the iterative Born cross-section unfolding with the
# vacuum-polarisation factor belong to the ROOT analysis.
# =============================================================================

### ------------------------------- Datasets ------------------------------ ###
# Scan points between 3.808 and 4.951 GeV that have both real data and
# inclusive MC samples available (BOSS 703 / 705 / 706 / 707).
scan_points = %w[
  703_4009 703_4180 703_4190 703_4200 703_4210 703_4220 703_4230
  703_4237 703_4246 703_4260 703_4270 703_4280 703_4360 703_4420 703_4600
  705_4130 705_4160
  706_4610 706_4620 706_4640 706_4660 706_4680 706_4700
  707_4740 707_4750 707_4780 707_4840 707_4914 707_4946
]

data_samples  = scan_points.map { |name| DatasetManager.real_data.find(name) }
incMC_samples = scan_points.map { |name| DatasetManager.inclusive_mc.find(name) }

### ------------------------------ Decay cards ---------------------------- ###
# A Born cross section measured as a function of sqrt(s), requiring ISR /
# vacuum-polarisation factors from the generator, is generated with ConExc.
# The DSL auto-detects the literal token "ConExc", switches to the no-KKMC
# template and injects "Particle vpho <ECMS> 0.0" per energy point, so
# Particle vpho is deliberately omitted for this multi-energy scan.

# (1) e+e- -> J/psi X (prompt), J/psi -> mu+ mu-   [PDG 443]
decay_card_jpsi_mumu = <<~DECAYCARD
    Decay vpho
    1.000 ConExc 443;
    Enddecay

    Decay J/psi
    1.000 mu+ mu- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# (2) ISR e+e- -> gamma_ISR J/psi, J/psi -> mu+ mu-   [443 + 22]
decay_card_isr_jpsi_mumu = <<~DECAYCARD
    Decay vpho
    1.000 ConExc 443 22;
    Enddecay

    Decay J/psi
    1.000 mu+ mu- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# (3) e+e- -> psi(3686) X (prompt), psi(2S) -> pi+ pi- J/psi, J/psi -> mu+ mu-   [100443]
decay_card_psip_mumu = <<~DECAYCARD
    Decay vpho
    1.000 ConExc 100443;
    Enddecay

    Decay psi(2S)
    1.000 pi+ pi- J/psi PHSP;
    Enddecay

    Decay J/psi
    1.000 mu+ mu- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# (4) e+e- -> psi(3686) X (prompt), psi(2S) -> pi+ pi- J/psi, J/psi -> e+ e-
decay_card_psip_ee = <<~DECAYCARD
    Decay vpho
    1.000 ConExc 100443;
    Enddecay

    Decay psi(2S)
    1.000 pi+ pi- J/psi PHSP;
    Enddecay

    Decay J/psi
    1.000 e+ e- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# (5) ISR e+e- -> gamma_ISR psi(3686), psi(2S) -> pi+ pi- J/psi, J/psi -> mu+ mu-
decay_card_isr_psip_mumu = <<~DECAYCARD
    Decay vpho
    1.000 ConExc 100443 22;
    Enddecay

    Decay psi(2S)
    1.000 pi+ pi- J/psi PHSP;
    Enddecay

    Decay J/psi
    1.000 mu+ mu- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# (6) ISR e+e- -> gamma_ISR psi(3686), psi(2S) -> pi+ pi- J/psi, J/psi -> e+ e-
decay_card_isr_psip_ee = <<~DECAYCARD
    Decay vpho
    1.000 ConExc 100443 22;
    Enddecay

    Decay psi(2S)
    1.000 pi+ pi- J/psi PHSP;
    Enddecay

    Decay J/psi
    1.000 e+ e- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# (7) e+e- -> chi_c1 X, chi_c1 -> gamma J/psi, J/psi -> mu+ mu-   [20443]
decay_card_chic1_mumu = <<~DECAYCARD
    Decay vpho
    1.000 ConExc 20443;
    Enddecay

    Decay chi_c1
    1.000 gamma J/psi PHSP;
    Enddecay

    Decay J/psi
    1.000 mu+ mu- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

### --------------------- Exclusive MC (one per mode) --------------------- ###
# Same signal MC is generated at every scan energy point -> use
# create_exclusive_mc_for over the list of real datasets (100k events each).

exmc_jpsi_mumu = DatasetManager.create_exclusive_mc_for(data_samples) do |config|
  config.sample_name   = "exmc_jpsi_mumu"
  config.events        = 100_000
  config.decay_card    = decay_card_jpsi_mumu
  config.cross_section = :default
end

exmc_isr_jpsi_mumu = DatasetManager.create_exclusive_mc_for(data_samples) do |config|
  config.sample_name   = "exmc_isr_jpsi_mumu"
  config.events        = 100_000
  config.decay_card    = decay_card_isr_jpsi_mumu
  config.cross_section = :default
end

exmc_psip_mumu = DatasetManager.create_exclusive_mc_for(data_samples) do |config|
  config.sample_name   = "exmc_psip_mumu"
  config.events        = 100_000
  config.decay_card    = decay_card_psip_mumu
  config.cross_section = :default
end

exmc_psip_ee = DatasetManager.create_exclusive_mc_for(data_samples) do |config|
  config.sample_name   = "exmc_psip_ee"
  config.events        = 100_000
  config.decay_card    = decay_card_psip_ee
  config.cross_section = :default
end

exmc_isr_psip_mumu = DatasetManager.create_exclusive_mc_for(data_samples) do |config|
  config.sample_name   = "exmc_isr_psip_mumu"
  config.events        = 100_000
  config.decay_card    = decay_card_isr_psip_mumu
  config.cross_section = :default
end

exmc_isr_psip_ee = DatasetManager.create_exclusive_mc_for(data_samples) do |config|
  config.sample_name   = "exmc_isr_psip_ee"
  config.events        = 100_000
  config.decay_card    = decay_card_isr_psip_ee
  config.cross_section = :default
end

exmc_chic1_mumu = DatasetManager.create_exclusive_mc_for(data_samples) do |config|
  config.sample_name   = "exmc_chic1_mumu"
  config.events        = 100_000
  config.decay_card    = decay_card_chic1_mumu
  config.cross_section = :default
end

# =============================================================================
# Event selection (BOSS)
# =============================================================================

# -----------------------------------------------------------------------------
# (A) e+e- -> J/psi X (prompt + ISR), J/psi -> mu+ mu-  :  final state mu+ mu-
# -----------------------------------------------------------------------------
alg_jpsi = Algorithm.new("InclusiveJpsi")
alg_jpsi.set_header(["InclusiveJpsiAlg/InclusiveJpsi.h"])
        .set_constant({"ECMS" => [:double, 4.260]})   # representative scan value;
                                                       # per-point vpho energy is injected by ConExc
        .set_alias({"std::vector<double>" => "Vdouble"})

sel_jpsi = Selection.new
  .select_track {                     # charged-track quality + multiplicity
    cos_theta 0.93                    # |cos(theta)| < 0.93
    Vz        10.0                    # |Vz| < 10 cm
    Vr        1.0                     # Vr  < 1 cm
    nChrp     ">=1"                   # at least one positive track
    nChrn     ">=1"                   # at least one negative track
    nNet      "==0"                   # net charge zero
  }
  .select_photon {                    # photons are selected but not required here
    tdc_emc_start     0               # TDC window 0 - 700 ns
    tdc_emc_end       14
    angle_to_track    20.0            # > 20 deg from the nearest charged track
    energyThreshold_b 0.025           # E > 25 MeV (barrel)
    energyThreshold_e 0.050           # E > 50 MeV (endcap)
  }
  .pid(method: :probability) {        # PID: probability method, cut 0.001
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6   # p>1.0 -> lepton; EMC>0.6 -> e, else mu
    nlp ">=1"                         # at least one positive lepton
    nlm ">=1"                         # at least one negative lepton
  }
  # Kalman fit: constrain the di-lepton mass to the nominal J/psi (chi2 < 50)
  .kalman_kinematic_fit([:lp, :lm]) {
    invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)
    chi2_cut 50
    njpsi    ">=1"                    # at least one J/psi candidate
  }
  # Nominal 4C kinematic fit on mu+ mu- (chi2 < 200, loose; tight cut in ROOT)
  .kinematic_fit([:lp, :lm]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_jpsi.note(:inclusive_production,
  "the undetected hadronic system X recoiling against the J/psi is not " \
  "reconstructed; the signal is the J/psi -> l+l- candidate itself and the " \
  "Born cross section is extracted from the fitted M(ll) yield in ROOT")

alg_jpsi.with_decay_card(decay_card_jpsi_mumu).apply(sel_jpsi)

# -----------------------------------------------------------------------------
# (B) e+e- -> psi(3686) X (prompt + ISR), psi(2S) -> pi+ pi- J/psi,
#     J/psi -> mu+ mu-  OR  J/psi -> e+ e-  :  final state pi+ pi- l+ l-
# -----------------------------------------------------------------------------
alg_psip = Algorithm.new("InclusivePsi2S")
alg_psip.set_header(["InclusivePsi2SAlg/InclusivePsi2S.h"])
        .set_constant({"ECMS" => [:double, 4.260]})
        .set_alias({"std::vector<double>" => "Vdouble"})

sel_psip = Selection.new
  .select_track {                     # charged-track quality + multiplicity
    cos_theta 0.93                    # |cos(theta)| < 0.93
    Vz        10.0                    # |Vz| < 10 cm
    Vr        1.0                     # Vr  < 1 cm
    nChrp     ">=2"                   # at least two positive tracks (pi+ and l+)
    nChrn     ">=2"                   # at least two negative tracks (pi- and l-)
    nNet      "==0"                   # net charge zero
  }
  .select_photon {                    # photons selected but not required in this channel
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    20.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
  }
  .pid(method: :probability) {        # PID: probability method, cut 0.001
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    nlp ">=1"                         # at least one positive lepton (from J/psi)
    nlm ">=1"                         # at least one negative lepton (from J/psi)
  }
  .remove([:lp <= :chrgp, :lm <= :chrgn])   # remove the identified leptons
  .assign({:chrgp => :pip, :chrgn => :pim}) # remaining tracks are the pi+ pi- pair
  # Kalman fit: constrain the di-lepton mass to the nominal J/psi (chi2 < 50)
  .kalman_kinematic_fit([:lp, :lm]) {
    invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)
    chi2_cut 50
    njpsi    ">=1"
  }
  # Nominal 4C kinematic fit on pi+ pi- l+ l- (chi2 < 200, loose; tight cut in ROOT)
  .kinematic_fit([:pip, :pim, :lp, :lm]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_psip.with_decay_card(decay_card_psip_mumu).apply(sel_psip)

# -----------------------------------------------------------------------------
# (C) e+e- -> chi_c1 X, chi_c1 -> gamma J/psi, J/psi -> mu+ mu-
#     final state gamma mu+ mu-
# -----------------------------------------------------------------------------
alg_chic1 = Algorithm.new("InclusiveChiC1")
alg_chic1.set_header(["InclusiveChiC1Alg/InclusiveChiC1.h"])
         .set_constant({"ECMS" => [:double, 4.260]})
         .set_alias({"std::vector<double>" => "Vdouble"})

sel_chic1 = Selection.new
  .select_track {                     # charged-track quality + multiplicity
    cos_theta 0.93                    # |cos(theta)| < 0.93
    Vz        10.0                    # |Vz| < 10 cm
    Vr        1.0                     # Vr  < 1 cm
    nChrp     ">=1"                   # at least one positive track
    nChrn     ">=1"                   # at least one negative track
    nNet      "==0"                   # net charge zero
  }
  .select_photon {                    # radiative photon is required in this channel
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    20.0            # > 20 deg from the nearest charged track
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=1"           # at least one photon
  }
  .pid(method: :probability) {        # PID: probability method, cut 0.001
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    nlp ">=1"                         # at least one positive lepton
    nlm ">=1"                         # at least one negative lepton
  }
  # Kalman fit: constrain the di-lepton mass to the nominal J/psi (chi2 < 50)
  .kalman_kinematic_fit([:lp, :lm]) {
    invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)
    chi2_cut 50
    njpsi    ">=1"
  }
  # Nominal 4C kinematic fit on gamma mu+ mu- (chi2 < 200, loose; tight cut in ROOT)
  .kinematic_fit([:gamma, :lp, :lm]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
  # Alternative 4C fit: same tracks, additionally constraining the di-muon to
  # the J/psi mass. No chi2_cut / not nominal -> its chi2 is stored for the
  # ROOT-level competing-hypothesis veto.
  .kinematic_fit([:gamma, :lp, :lm]) {
    constrain_four_momentum
    invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)
  }

alg_chic1.note(:photon_isolation,
  "the chi_c1 channel requires the radiative photon to be isolated by more " \
  "than 20 deg from BOTH muon tracks; select_photon applies the 20 deg cut to " \
  "the nearest charged track (a generic isolation primitive for muon tracks is " \
  "not available in the DSL), the explicit two-muon isolation is enforced at " \
  "the ROOT level")

alg_chic1.with_decay_card(decay_card_chic1_mumu).apply(sel_chic1)

### ------------------------------ Execution ------------------------------ ###
# Real data + inclusive MC at every energy point, plus the mode-specific
# exclusive MC samples.
alg_jpsi.execute_on(data_samples + incMC_samples +
                    exmc_jpsi_mumu + exmc_isr_jpsi_mumu)

alg_psip.execute_on(data_samples + incMC_samples +
                    exmc_psip_mumu + exmc_psip_ee +
                    exmc_isr_psip_mumu + exmc_isr_psip_ee)

alg_chic1.execute_on(data_samples + incMC_samples + exmc_chic1_mumu)