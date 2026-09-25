# ============================================================================
#  BOSS part : psi(4260) -> gamma X(3872) on the 12 BOSS-703 energy-scan points
#              (4.178 - 4.288 GeV), real data + matching inclusive MC
#              X(3872) -> pi0 chi_c0 / pi+pi- chi_c0 / pi0pi0 chi_c0
#              chi_c0  -> pi+pi-, K+K-, pi+pi-pi+pi-, pi+pi-K+K-, pi+pi-pi0pi0
#              normalisation : X(3872) -> pi+pi- J/psi, J/psi -> l+l-
# ============================================================================

### ------------------------------- datasets ------------------------------- ###
# 12 BOSS 703 scan points between 4.178 and 4.288 GeV
scan_point_names = %w[
  703_4180 703_4190 703_4200 703_4210 703_4220 703_4230
  703_4237 703_4245 703_4246 703_4260 703_4270 703_4280
]
scan_data  = scan_point_names.map { |n| DatasetManager.real_data.find(n) }     # real data
scan_incMC = scan_point_names.map { |n| DatasetManager.inclusive_mc.find(n) } # matching inclusive MC

### ---------------------------- decay cards (EvtGen) ---------------------- ###
# (1) X(3872) -> pi0 chi_c0, chi_c0 -> pi+ pi-
dc_pi0chic0_pipi = <<~DECAYCARD
  Decay psi(4260)
  1.000 gamma X(3872) PHSP;
  Enddecay

  Decay X(3872)
  1.000 pi0 chi_c0 PHSP;
  Enddecay

  Decay chi_c0
  1.000 pi+ pi- PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# (2) X(3872) -> pi0 chi_c0, chi_c0 -> K+ K-
dc_pi0chic0_KK = <<~DECAYCARD
  Decay psi(4260)
  1.000 gamma X(3872) PHSP;
  Enddecay

  Decay X(3872)
  1.000 pi0 chi_c0 PHSP;
  Enddecay

  Decay chi_c0
  1.000 K+ K- PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# (3) X(3872) -> pi0 chi_c0, chi_c0 -> pi+ pi- pi+ pi-
dc_pi0chic0_4pi = <<~DECAYCARD
  Decay psi(4260)
  1.000 gamma X(3872) PHSP;
  Enddecay

  Decay X(3872)
  1.000 pi0 chi_c0 PHSP;
  Enddecay

  Decay chi_c0
  1.000 pi+ pi- pi+ pi- PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# (4) X(3872) -> pi0 chi_c0, chi_c0 -> pi+ pi- K+ K-
dc_pi0chic0_pipiKK = <<~DECAYCARD
  Decay psi(4260)
  1.000 gamma X(3872) PHSP;
  Enddecay

  Decay X(3872)
  1.000 pi0 chi_c0 PHSP;
  Enddecay

  Decay chi_c0
  1.000 pi+ pi- K+ K- PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# (5) X(3872) -> pi0 chi_c0, chi_c0 -> pi+ pi- pi0 pi0
dc_pi0chic0_pipi2pi0 = <<~DECAYCARD
  Decay psi(4260)
  1.000 gamma X(3872) PHSP;
  Enddecay

  Decay X(3872)
  1.000 pi0 chi_c0 PHSP;
  Enddecay

  Decay chi_c0
  1.000 pi+ pi- pi0 pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# (6) X(3872) -> pi+ pi- chi_c0, chi_c0 -> pi+ pi-
dc_pipichic0_pipi = <<~DECAYCARD
  Decay psi(4260)
  1.000 gamma X(3872) PHSP;
  Enddecay

  Decay X(3872)
  1.000 pi+ pi- chi_c0 PHSP;
  Enddecay

  Decay chi_c0
  1.000 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# (7) X(3872) -> pi+ pi- chi_c0, chi_c0 -> K+ K-
dc_pipichic0_KK = <<~DECAYCARD
  Decay psi(4260)
  1.000 gamma X(3872) PHSP;
  Enddecay

  Decay X(3872)
  1.000 pi+ pi- chi_c0 PHSP;
  Enddecay

  Decay chi_c0
  1.000 K+ K- PHSP;
  Enddecay

  End
DECAYCARD

# (8) X(3872) -> pi+ pi- chi_c0, chi_c0 -> pi+ pi- pi+ pi-
dc_pipichic0_4pi = <<~DECAYCARD
  Decay psi(4260)
  1.000 gamma X(3872) PHSP;
  Enddecay

  Decay X(3872)
  1.000 pi+ pi- chi_c0 PHSP;
  Enddecay

  Decay chi_c0
  1.000 pi+ pi- pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# (9) X(3872) -> pi+ pi- chi_c0, chi_c0 -> pi+ pi- K+ K-
dc_pipichic0_pipiKK = <<~DECAYCARD
  Decay psi(4260)
  1.000 gamma X(3872) PHSP;
  Enddecay

  Decay X(3872)
  1.000 pi+ pi- chi_c0 PHSP;
  Enddecay

  Decay chi_c0
  1.000 pi+ pi- K+ K- PHSP;
  Enddecay

  End
DECAYCARD

# (10) X(3872) -> pi+ pi- chi_c0, chi_c0 -> pi+ pi- pi0 pi0
dc_pipichic0_pipi2pi0 = <<~DECAYCARD
  Decay psi(4260)
  1.000 gamma X(3872) PHSP;
  Enddecay

  Decay X(3872)
  1.000 pi+ pi- chi_c0 PHSP;
  Enddecay

  Decay chi_c0
  1.000 pi+ pi- pi0 pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# (11) X(3872) -> pi0 pi0 chi_c0, chi_c0 -> pi+ pi-
dc_2pi0chic0_pipi = <<~DECAYCARD
  Decay psi(4260)
  1.000 gamma X(3872) PHSP;
  Enddecay

  Decay X(3872)
  1.000 pi0 pi0 chi_c0 PHSP;
  Enddecay

  Decay chi_c0
  1.000 pi+ pi- PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# (12) X(3872) -> pi0 pi0 chi_c0, chi_c0 -> K+ K-
dc_2pi0chic0_KK = <<~DECAYCARD
  Decay psi(4260)
  1.000 gamma X(3872) PHSP;
  Enddecay

  Decay X(3872)
  1.000 pi0 pi0 chi_c0 PHSP;
  Enddecay

  Decay chi_c0
  1.000 K+ K- PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# (13) X(3872) -> pi0 pi0 chi_c0, chi_c0 -> pi+ pi- pi+ pi-
dc_2pi0chic0_4pi = <<~DECAYCARD
  Decay psi(4260)
  1.000 gamma X(3872) PHSP;
  Enddecay

  Decay X(3872)
  1.000 pi0 pi0 chi_c0 PHSP;
  Enddecay

  Decay chi_c0
  1.000 pi+ pi- pi+ pi- PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# (14) X(3872) -> pi0 pi0 chi_c0, chi_c0 -> pi+ pi- K+ K-
dc_2pi0chic0_pipiKK = <<~DECAYCARD
  Decay psi(4260)
  1.000 gamma X(3872) PHSP;
  Enddecay

  Decay X(3872)
  1.000 pi0 pi0 chi_c0 PHSP;
  Enddecay

  Decay chi_c0
  1.000 pi+ pi- K+ K- PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# (15) X(3872) -> pi0 pi0 chi_c0, chi_c0 -> pi+ pi- pi0 pi0
dc_2pi0chic0_pipi2pi0 = <<~DECAYCARD
  Decay psi(4260)
  1.000 gamma X(3872) PHSP;
  Enddecay

  Decay X(3872)
  1.000 pi0 pi0 chi_c0 PHSP;
  Enddecay

  Decay chi_c0
  1.000 pi+ pi- pi0 pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# (16) normalisation channel : X(3872) -> pi+ pi- J/psi, J/psi -> e+ e-
#      (the l+l- channel is handled by the combined lepton lists :lp/:lm,
#       the muon mode uses the same card with mu+ mu- as J/psi daughters)
dc_pipiJpsi_ll = <<~DECAYCARD
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

### ------------------------ exclusive MC samples -------------------------- ###
# 100k events of psi(4260) -> gamma X(3872), X(3872) -> pi0 chi_c0,
# chi_c0 -> pi+ pi-, pi0 -> gamma gamma, generated at each of the scan points
exMC_pi0chic0_pipi = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name   = "exmc_gamX3872_pi0chic0_pipi"  # per-point names like ..._703_4180
  config.events        = 100_000
  config.decay_card    = dc_pi0chic0_pipi
  config.cross_section = :default
end

### --------------------------- event selection ---------------------------- ###
# Single ECMS constant per algorithm (nominal psi(4260) energy); the individual
# scan points lie between 4.178 and 4.288 GeV.
ecms = 4.260

# Common charged-track + photon selection; the three arguments are the
# multiplicity requirements (positive tracks, negative tracks, photons).
def build_selection(n_chrp, n_chrn, n_gam)
  Selection.new
    .select_track {
      cos_theta 0.93        # |cos(theta)| < 0.93
      Vz        10.0        # |Vz| < 10 cm
      Vr        1.0         # Vr < 1 cm
      nChrp     n_chrp
      nChrn     n_chrn
      nNet      "==0"       # zero net charge
    }
    .select_photon {
      tdc_emc_start     0     # EMC timing window
      tdc_emc_end       14
      angle_to_track    10.0  # > 10 degrees from any charged track
      energyThreshold_b 0.025 # 25 MeV in the barrel
      energyThreshold_e 0.050 # 50 MeV in the endcap
      nGam              n_gam
    }
end

# Kaon-mode PID: identify K+/K- against pions, remove them from the charged
# lists and treat the remaining tracks as pions.
def kaon_pid(selection)
  selection
    .pid(method: :probability) {
      prob_cut 0.001
      identify :kaon, against: [:pion, :proton]
      nkp ">=1"
      nkm ">=1"
    }
    .remove([:kp <= :chrgp, :km <= :chrgn])       # drop identified kaons
    .assign({:chrgp => :pip, :chrgn => :pim})     # remaining tracks -> pions
end

# ---------------------------------------------------------------------------
# (1) X(3872) -> pi0 chi_c0, chi_c0 -> pi+ pi-
alg_pi0chic0_pipi = Algorithm.new("X3872Pi0Chic0ToPiPi")
alg_pi0chic0_pipi.set_header(["X3872Pi0Chic0ToPiPiAlg/X3872Pi0Chic0ToPiPi.h"])
                 .set_constant({"ECMS" => [:double, ecms]})
sel_pi0chic0_pipi = build_selection(">=1", ">=1", ">=2")
  .assign({:chrgp => :pip, :chrgn => :pim})                    # pions without PID
  .kalman_kinematic_fit([:gamma, :gamma]) {                    # 1-C : gamma gamma -> pi0
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:gamma, :pi0, :pip, :pim]) {                 # 4C fit
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_pi0chic0_pipi.with_decay_card(dc_pi0chic0_pipi).apply(sel_pi0chic0_pipi)

# ---------------------------------------------------------------------------
# (2) X(3872) -> pi0 chi_c0, chi_c0 -> K+ K-
alg_pi0chic0_KK = Algorithm.new("X3872Pi0Chic0ToKK")
alg_pi0chic0_KK.set_header(["X3872Pi0Chic0ToKKAlg/X3872Pi0Chic0ToKK.h"])
               .set_constant({"ECMS" => [:double, ecms]})
sel_pi0chic0_KK = kaon_pid(build_selection(">=1", ">=1", ">=2"))
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:gamma, :pi0, :kp, :km]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_pi0chic0_KK.with_decay_card(dc_pi0chic0_KK).apply(sel_pi0chic0_KK)

# ---------------------------------------------------------------------------
# (3) X(3872) -> pi0 chi_c0, chi_c0 -> pi+ pi- pi+ pi-
alg_pi0chic0_4pi = Algorithm.new("X3872Pi0Chic0To4Pi")
alg_pi0chic0_4pi.set_header(["X3872Pi0Chic0To4PiAlg/X3872Pi0Chic0To4Pi.h"])
                .set_constant({"ECMS" => [:double, ecms]})
sel_pi0chic0_4pi = build_selection(">=2", ">=2", ">=2")
  .assign({:chrgp => :pip, :chrgn => :pim})
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:gamma, :pi0, :pip, :pip, :pim, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_pi0chic0_4pi.with_decay_card(dc_pi0chic0_4pi).apply(sel_pi0chic0_4pi)

# ---------------------------------------------------------------------------
# (4) X(3872) -> pi0 chi_c0, chi_c0 -> pi+ pi- K+ K-
alg_pi0chic0_pipiKK = Algorithm.new("X3872Pi0Chic0ToPiPiKK")
alg_pi0chic0_pipiKK.set_header(["X3872Pi0Chic0ToPiPiKKAlg/X3872Pi0Chic0ToPiPiKK.h"])
                   .set_constant({"ECMS" => [:double, ecms]})
sel_pi0chic0_pipiKK = kaon_pid(build_selection(">=2", ">=2", ">=2"))
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:gamma, :pi0, :pip, :pim, :kp, :km]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_pi0chic0_pipiKK.with_decay_card(dc_pi0chic0_pipiKK).apply(sel_pi0chic0_pipiKK)

# ---------------------------------------------------------------------------
# (5) X(3872) -> pi0 chi_c0, chi_c0 -> pi+ pi- pi0 pi0   (3 pi0 -> >=6 photons)
alg_pi0chic0_pipi2pi0 = Algorithm.new("X3872Pi0Chic0ToPiPi2Pi0")
alg_pi0chic0_pipi2pi0.set_header(["X3872Pi0Chic0ToPiPi2Pi0Alg/X3872Pi0Chic0ToPiPi2Pi0.h"])
                     .set_constant({"ECMS" => [:double, ecms]})
sel_pi0chic0_pipi2pi0 = build_selection(">=1", ">=1", ">=6")
  .assign({:chrgp => :pip, :chrgn => :pim})
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=3"
  }
  .kinematic_fit([:gamma, :pi0, :pi0, :pi0, :pip, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_pi0chic0_pipi2pi0.with_decay_card(dc_pi0chic0_pipi2pi0).apply(sel_pi0chic0_pipi2pi0)

# ---------------------------------------------------------------------------
# (6) X(3872) -> pi+ pi- chi_c0, chi_c0 -> pi+ pi-   (only the prompt photon)
alg_pipichic0_pipi = Algorithm.new("X3872PiPiChic0ToPiPi")
alg_pipichic0_pipi.set_header(["X3872PiPiChic0ToPiPiAlg/X3872PiPiChic0ToPiPi.h"])
                  .set_constant({"ECMS" => [:double, ecms]})
sel_pipichic0_pipi = build_selection(">=2", ">=2", ">=1")
  .assign({:chrgp => :pip, :chrgn => :pim})
  .kinematic_fit([:gamma, :pip, :pip, :pim, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_pipichic0_pipi.with_decay_card(dc_pipichic0_pipi).apply(sel_pipichic0_pipi)

# ---------------------------------------------------------------------------
# (7) X(3872) -> pi+ pi- chi_c0, chi_c0 -> K+ K-
alg_pipichic0_KK = Algorithm.new("X3872PiPiChic0ToKK")
alg_pipichic0_KK.set_header(["X3872PiPiChic0ToKKAlg/X3872PiPiChic0ToKK.h"])
                .set_constant({"ECMS" => [:double, ecms]})
sel_pipichic0_KK = kaon_pid(build_selection(">=2", ">=2", ">=1"))
  .kinematic_fit([:gamma, :pip, :pim, :kp, :km]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_pipichic0_KK.with_decay_card(dc_pipichic0_KK).apply(sel_pipichic0_KK)

# ---------------------------------------------------------------------------
# (8) X(3872) -> pi+ pi- chi_c0, chi_c0 -> pi+ pi- pi+ pi-  (six-pion mode)
alg_pipichic0_4pi = Algorithm.new("X3872PiPiChic0To4Pi")
alg_pipichic0_4pi.set_header(["X3872PiPiChic0To4PiAlg/X3872PiPiChic0To4Pi.h"])
                 .set_constant({"ECMS" => [:double, ecms]})
sel_pipichic0_4pi = build_selection(">=3", ">=3", ">=1")
  .assign({:chrgp => :pip, :chrgn => :pim})
  .kinematic_fit([:gamma, :pip, :pip, :pip, :pim, :pim, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_pipichic0_4pi.with_decay_card(dc_pipichic0_4pi).apply(sel_pipichic0_4pi)

# ---------------------------------------------------------------------------
# (9) X(3872) -> pi+ pi- chi_c0, chi_c0 -> pi+ pi- K+ K-
alg_pipichic0_pipiKK = Algorithm.new("X3872PiPiChic0ToPiPiKK")
alg_pipichic0_pipiKK.set_header(["X3872PiPiChic0ToPiPiKKAlg/X3872PiPiChic0ToPiPiKK.h"])
                    .set_constant({"ECMS" => [:double, ecms]})
sel_pipichic0_pipiKK = kaon_pid(build_selection(">=3", ">=3", ">=1"))
  .kinematic_fit([:gamma, :pip, :pip, :pim, :pim, :kp, :km]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_pipichic0_pipiKK.with_decay_card(dc_pipichic0_pipiKK).apply(sel_pipichic0_pipiKK)

# ---------------------------------------------------------------------------
# (10) X(3872) -> pi+ pi- chi_c0, chi_c0 -> pi+ pi- pi0 pi0   (2 pi0 -> >=4 photons)
alg_pipichic0_pipi2pi0 = Algorithm.new("X3872PiPiChic0ToPiPi2Pi0")
alg_pipichic0_pipi2pi0.set_header(["X3872PiPiChic0ToPiPi2Pi0Alg/X3872PiPiChic0ToPiPi2Pi0.h"])
                      .set_constant({"ECMS" => [:double, ecms]})
sel_pipichic0_pipi2pi0 = build_selection(">=2", ">=2", ">=4")
  .assign({:chrgp => :pip, :chrgn => :pim})
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=2"
  }
  .kinematic_fit([:gamma, :pi0, :pi0, :pip, :pip, :pim, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_pipichic0_pipi2pi0.with_decay_card(dc_pipichic0_pipi2pi0).apply(sel_pipichic0_pipi2pi0)

# ---------------------------------------------------------------------------
# (11) X(3872) -> pi0 pi0 chi_c0, chi_c0 -> pi+ pi-   (2 pi0 -> >=4 photons)
alg_2pi0chic0_pipi = Algorithm.new("X38722Pi0Chic0ToPiPi")
alg_2pi0chic0_pipi.set_header(["X38722Pi0Chic0ToPiPiAlg/X38722Pi0Chic0ToPiPi.h"])
                  .set_constant({"ECMS" => [:double, ecms]})
sel_2pi0chic0_pipi = build_selection(">=1", ">=1", ">=4")
  .assign({:chrgp => :pip, :chrgn => :pim})
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=2"
  }
  .kinematic_fit([:gamma, :pi0, :pi0, :pip, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_2pi0chic0_pipi.with_decay_card(dc_2pi0chic0_pipi).apply(sel_2pi0chic0_pipi)

# ---------------------------------------------------------------------------
# (12) X(3872) -> pi0 pi0 chi_c0, chi_c0 -> K+ K-
alg_2pi0chic0_KK = Algorithm.new("X38722Pi0Chic0ToKK")
alg_2pi0chic0_KK.set_header(["X38722Pi0Chic0ToKKAlg/X38722Pi0Chic0ToKK.h"])
                .set_constant({"ECMS" => [:double, ecms]})
sel_2pi0chic0_KK = kaon_pid(build_selection(">=1", ">=1", ">=4"))
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=2"
  }
  .kinematic_fit([:gamma, :pi0, :pi0, :kp, :km]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_2pi0chic0_KK.with_decay_card(dc_2pi0chic0_KK).apply(sel_2pi0chic0_KK)

# ---------------------------------------------------------------------------
# (13) X(3872) -> pi0 pi0 chi_c0, chi_c0 -> pi+ pi- pi+ pi-
alg_2pi0chic0_4pi = Algorithm.new("X38722Pi0Chic0To4Pi")
alg_2pi0chic0_4pi.set_header(["X38722Pi0Chic0To4PiAlg/X38722Pi0Chic0To4Pi.h"])
                 .set_constant({"ECMS" => [:double, ecms]})
sel_2pi0chic0_4pi = build_selection(">=2", ">=2", ">=4")
  .assign({:chrgp => :pip, :chrgn => :pim})
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=2"
  }
  .kinematic_fit([:gamma, :pi0, :pi0, :pip, :pip, :pim, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_2pi0chic0_4pi.with_decay_card(dc_2pi0chic0_4pi).apply(sel_2pi0chic0_4pi)

# ---------------------------------------------------------------------------
# (14) X(3872) -> pi0 pi0 chi_c0, chi_c0 -> pi+ pi- K+ K-
alg_2pi0chic0_pipiKK = Algorithm.new("X38722Pi0Chic0ToPiPiKK")
alg_2pi0chic0_pipiKK.set_header(["X38722Pi0Chic0ToPiPiKKAlg/X38722Pi0Chic0ToPiPiKK.h"])
                    .set_constant({"ECMS" => [:double, ecms]})
sel_2pi0chic0_pipiKK = kaon_pid(build_selection(">=2", ">=2", ">=4"))
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=2"
  }
  .kinematic_fit([:gamma, :pi0, :pi0, :pip, :pim, :kp, :km]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_2pi0chic0_pipiKK.with_decay_card(dc_2pi0chic0_pipiKK).apply(sel_2pi0chic0_pipiKK)

# ---------------------------------------------------------------------------
# (15) X(3872) -> pi0 pi0 chi_c0, chi_c0 -> pi+ pi- pi0 pi0  (4 pi0 -> >=8 photons)
alg_2pi0chic0_pipi2pi0 = Algorithm.new("X38722Pi0Chic0ToPiPi2Pi0")
alg_2pi0chic0_pipi2pi0.set_header(["X38722Pi0Chic0ToPiPi2Pi0Alg/X38722Pi0Chic0ToPiPi2Pi0.h"])
                      .set_constant({"ECMS" => [:double, ecms]})
sel_2pi0chic0_pipi2pi0 = build_selection(">=1", ">=1", ">=8")
  .assign({:chrgp => :pip, :chrgn => :pim})
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=4"
  }
  .kinematic_fit([:gamma, :pi0, :pi0, :pi0, :pi0, :pip, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_2pi0chic0_pipi2pi0.with_decay_card(dc_2pi0chic0_pipi2pi0).apply(sel_2pi0chic0_pipi2pi0)

# ---------------------------------------------------------------------------
# (16) normalisation : X(3872) -> pi+ pi- J/psi, J/psi -> l+ l-
alg_pipiJpsi_ll = Algorithm.new("X3872PiPiJpsiToLL")
alg_pipiJpsi_ll.set_header(["X3872PiPiJpsiToLLAlg/X3872PiPiJpsiToLL.h"])
               .set_constant({"ECMS" => [:double, ecms]})
sel_pipiJpsi_ll = build_selection(">=2", ">=2", ">=1")
  .pid(method: :probability) {
    # high-momentum tracks are leptons; E/p (EMC energy > 0.85 GeV) separates e from mu
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.85
    identify :pion, against: [:kaon]        # the two pions recoiling against the J/psi
    nlp  "==1"                              # one l+
    nlm  "==1"                              # one l-
    npip ">=1"                              # at least one pi+
    npim ">=1"                              # at least one pi-
  }
  .kinematic_fit([:gamma, :lp, :lm, :pip, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
alg_pipiJpsi_ll.with_decay_card(dc_pipiJpsi_ll).apply(sel_pipiJpsi_ll)

### ------------------------------ execution -------------------------------- ###
all_algorithms = [
  alg_pi0chic0_pipi, alg_pi0chic0_KK, alg_pi0chic0_4pi,
  alg_pi0chic0_pipiKK, alg_pi0chic0_pipi2pi0,
  alg_pipichic0_pipi, alg_pipichic0_KK, alg_pipichic0_4pi,
  alg_pipichic0_pipiKK, alg_pipichic0_pipi2pi0,
  alg_2pi0chic0_pipi, alg_2pi0chic0_KK, alg_2pi0chic0_4pi,
  alg_2pi0chic0_pipiKK, alg_2pi0chic0_pipi2pi0,
  alg_pipiJpsi_ll
]

# ECMS is a single constant in the generated algorithm, while the data are taken
# at twelve different scan-point energies (4.178 - 4.288 GeV); the exact
# centre-of-mass energy of each job must be taken run by run.
all_algorithms.each do |a|
  a.note(:ecms_per_scan_point,
         "the analysis uses twelve BOSS 703 scan points between 4.178 and 4.288 GeV; " \
         "the ECMS constant of the generated algorithm is set to the nominal psi(4260) " \
         "value 4.260 GeV, while the per-point beam energy must be used for the 4C fits")
end

all_datasets = scan_data + scan_incMC + exMC_pi0chic0_pipi

alg_pi0chic0_pipi.execute_on(all_datasets)
alg_pi0chic0_KK.execute_on(all_datasets)
alg_pi0chic0_4pi.execute_on(all_datasets)
alg_pi0chic0_pipiKK.execute_on(all_datasets)
alg_pi0chic0_pipi2pi0.execute_on(all_datasets)
alg_pipichic0_pipi.execute_on(all_datasets)
alg_pipichic0_KK.execute_on(all_datasets)
alg_pipichic0_4pi.execute_on(all_datasets)
alg_pipichic0_pipiKK.execute_on(all_datasets)
alg_pipichic0_pipi2pi0.execute_on(all_datasets)
alg_2pi0chic0_pipi.execute_on(all_datasets)
alg_2pi0chic0_KK.execute_on(all_datasets)
alg_2pi0chic0_4pi.execute_on(all_datasets)
alg_2pi0chic0_pipiKK.execute_on(all_datasets)
alg_2pi0chic0_pipi2pi0.execute_on(all_datasets)
alg_pipiJpsi_ll.execute_on(all_datasets)