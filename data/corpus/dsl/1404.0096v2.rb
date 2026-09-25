# Observation of eta' -> pi+ pi- pi+ pi- and eta' -> pi+ pi- pi0 pi0
# from J/psi -> gamma eta'   [arXiv:1404.0096]
#
# Two independent reconstruction modes with different final states and
# different kinematic-fit hypotheses, therefore one Algorithm each (Rule T1).

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")   # 1.3 x 10^9 J/psi events (2009 + 2012)
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# ---------------------------------------------------------------- decay cards
# Mode A: J/psi -> gamma eta', eta' -> pi+ pi- pi+ pi-
decay_card_modeA = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta'  PHSP;
  Enddecay

  Decay eta'
  1.0000 pi+ pi- pi+ pi-  PHSP;
  Enddecay

  End
DECAYCARD

# Mode B: J/psi -> gamma eta', eta' -> pi+ pi- pi0 pi0, pi0 -> gamma gamma
decay_card_modeB = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta'  PHSP;
  Enddecay

  Decay eta'
  1.0000 pi+ pi- pi0 pi0  PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma  PHSP;
  Enddecay

  End
DECAYCARD

# Dominant background for Mode A: eta' -> pi+ pi- eta with eta -> gamma pi+ pi-,
# which produces the enhancement just below the eta' peak.
decay_card_bkgA = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta'  PHSP;
  Enddecay

  Decay eta'
  1.0000 pi+ pi- eta  PHSP;
  Enddecay

  Decay eta
  1.0000 gamma pi+ pi-  PHSP;
  Enddecay

  End
DECAYCARD

# Dominant backgrounds for Mode B: eta' -> pi+ pi- eta (eta -> pi0 pi0 pi0) and
# eta' -> pi0 pi0 eta (eta -> gamma pi+ pi-), plus eta' -> gamma omega.
decay_card_bkgB = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta'  PHSP;
  Enddecay

  Decay eta'
  1.0000 pi+ pi- eta  PHSP;
  Enddecay

  Decay eta
  1.0000 pi0 pi0 pi0  PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma  PHSP;
  Enddecay

  End
DECAYCARD

# Non-resonant / f2(1270) continuum background under J/psi -> gamma 4pi
decay_card_bkg_f2 = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma f2(1270)  PHSP;
  Enddecay

  Decay f2(1270)
  1.0000 pi+ pi- pi+ pi-  PHSP;
  Enddecay

  End
DECAYCARD

exMC_modeA   = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "jpsi_gamma_etap_4pi_charged"; c.related_dataset = jpsi_data
  c.events = 500_000; c.decay_card = decay_card_modeA; c.cross_section = :default
end
exMC_modeB   = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "jpsi_gamma_etap_pipipi0pi0"; c.related_dataset = jpsi_data
  c.events = 500_000; c.decay_card = decay_card_modeB; c.cross_section = :default
end
exMC_bkgA    = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "jpsi_gamma_etap_pipieta"; c.related_dataset = jpsi_data
  c.events = 500_000; c.decay_card = decay_card_bkgA; c.cross_section = :default
end
exMC_bkgB    = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "jpsi_gamma_etap_etapipi0eta"; c.related_dataset = jpsi_data
  c.events = 500_000; c.decay_card = decay_card_bkgB; c.cross_section = :default
end
exMC_bkg_f2  = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "jpsi_gamma_f2_4pi"; c.related_dataset = jpsi_data
  c.events = 500_000; c.decay_card = decay_card_bkg_f2; c.cross_section = :default
end

### Event selection (BOSS) ###

# ------------------------------------------------------------------ Mode A
# J/psi -> gamma eta', eta' -> pi+ pi- pi+ pi-   (4 charged tracks + >=1 photon)
alg_name = "JpsiGammaEtapTo4PiCharged"
algA = Algorithm.new(alg_name)
algA.set_header(["#{alg_name}Alg/#{alg_name}.h"])
    .set_constant({"ECMS" => [:double, 3.097]})

selA = Selection.new
selA.select_track do
      cos_theta 0.93    # |cos(theta)| < 0.93
      Vz        10.0    # within 10 cm of the IP along the beam direction
      Vr        1.0     # within 1 cm in the radial direction
      nChrp     "==2"   # 2 pi+
      nChrn     "==2"   # 2 pi-
      nNet      "==0"   # four charged tracks with zero net charge
    end
    .select_photon do
      tdc_emc_start     0      # EMC cluster timing (suppress noise / unrelated deposits)
      tdc_emc_end       14
      angle_to_track    10.0   # shower at least 10 degrees from the nearest charged track
      energyThreshold_b 0.025  # E > 25 MeV in the barrel (|cos(theta)| < 0.8)
      energyThreshold_e 0.050  # E > 50 MeV in the end cap (0.86 < |cos(theta)| < 0.92)
      nGam              ">=1"  # at least one photon from J/psi -> gamma eta'
    end
    .pid(method: :probability) do
      prob_cut 0.001
      identify :pion, against: [:kaon, :proton]
      npip ">=1"        # at least two oppositely charged tracks identified as pions
      npim ">=1"
    end
    .assign({chrgp: :pip, chrgn: :pim})  # remaining tracks assumed pions (no PID requirement)
    .kinematic_fit([:gamma, :pip, :pim, :pip, :pim]) do
      nominal
      vertex_fit([1, 2, 3, 4])  # vertex fit on the four charged tracks (loose chi2 at the IP)
      constrain_four_momentum   # 4C fit under the gamma pi+ pi- pi+ pi- hypothesis
      chi2_cut 200
    end
    # Competing hypothesis with an extra photon: chi2_4C(gamma 4pi) must be smaller
    # than chi2_4C(gamma gamma 4pi). No chi2_cut and no nominal -> the value is stored
    # for the ROOT-level veto chi2_4c_signal < chi2_4c_2gamma.
    .kinematic_fit([:gamma, :gamma, :pip, :pim, :pip, :pim]) do
      vertex_fit([2, 3, 4, 5])
      constrain_four_momentum
    end

algA.note(:kinematic_fit_chi2,
          "The published analysis applies chi2_4C < 35; the BOSS first pass keeps the " \
          "loose chi2_cut 200 and the tight cut is applied in the ROOT analysis. When more " \
          "than one photon candidate is present, the combination with the smallest chi2_4C " \
          "is retained, as done automatically by kinematic_fit.")
    .note(:background_veto,
          "Backgrounds that could fake the eta' peak were studied with a 1-billion-event " \
          "Lund-model J/psi inclusive sample: the enhancement below the eta' peak comes " \
          "from eta' -> pi+ pi- eta with eta -> gamma pi+ pi-, the region above 1 GeV/c^2 " \
          "from eta' -> pi+ pi- e+ e-, plus J/psi -> gamma f2(1270) with " \
          "f2(1270) -> pi+ pi- pi+ pi- and non-resonant J/psi -> gamma 4pi. None of them " \
          "peaks at the eta' mass.")
    .note(:signal_mc_model,
          "In addition to phase-space events, the signal MC is modelled with the ChPT+VMD " \
          "decay amplitudes of the paper's Ref. [8]; the ChPT+VMD model describes the " \
          "background-subtracted M(pi+pi-) distribution better than phase space and is " \
          "therefore used to determine the detection efficiency.")
    .note(:signal_yield_extraction,
          "The yield is obtained from an extended unbinned maximum likelihood fit to the " \
          "pi+ pi- pi+ pi- invariant mass. The signal PDF is the MC shape convoluted with " \
          "a Gaussian; the background PDFs comprise the dedicated-background MC shapes, a " \
          "Breit-Wigner f2(1270) tail convoluted with a Gaussian, and the phase-space shape.")

algA.with_decay_card(decay_card_modeA).apply(selA)
algA.execute_on([jpsi_data, jpsi_incMC, exMC_modeA, exMC_bkgA, exMC_bkg_f2])

# ------------------------------------------------------------------ Mode B
# J/psi -> gamma eta', eta' -> pi+ pi- pi0 pi0 (2 charged tracks + >=5 photons)
alg_name = "JpsiGammaEtapToPipipi0pi0"
algB = Algorithm.new(alg_name)
algB.set_header(["#{alg_name}Alg/#{alg_name}.h"])
    .set_constant({"ECMS" => [:double, 3.097]})

selB = Selection.new
selB.select_track do
      cos_theta 0.93
      Vz        10.0
      Vr        1.0
      nChrp     "==1"   # pi+
      nChrn     "==1"   # pi-
      nNet      "==0"   # two charged tracks with zero net charge
    end
    .select_photon do
      tdc_emc_start     0
      tdc_emc_end       14
      angle_to_track    10.0
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      nGam              ">=5"  # one prompt photon + four photons from the two pi0
    end
    .pid(method: :probability) do
      prob_cut 0.001
      identify :pion, against: [:kaon, :proton]
      npip "==1"        # both charged tracks identified as pions
      npim "==1"
    end
    # 1C fits on the pi0 candidates from photon pairs, mass-constrained to the pi0 mass
    .kalman_kinematic_fit([:gamma, :gamma]) do
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 50      # chi2_1C(gamma gamma) < 50
      npi0 ">=2"       # two pi0 candidates required
    end
    # 6C kinematic fit under the J/psi -> gamma pi+ pi- pi0 pi0 hypothesis
    # (four-momentum conservation plus the two pi0 mass constraints already
    # carried by the Kalman-reconstructed pi0 four-momenta).
    .kinematic_fit([:gamma, :pip, :pim, :pi0, :pi0]) do
      nominal
      vertex_fit([1, 2])
      constrain_four_momentum
      chi2_cut 200
    end
    # Competing hypothesis with an extra photon: chi2_6C must be smaller than that of
    # the gamma gamma pi+ pi- pi0 pi0 hypothesis, rejecting events with six photons.
    .kinematic_fit([:gamma, :gamma, :pip, :pim, :pi0, :pi0]) do
      vertex_fit([2, 3])
      constrain_four_momentum
    end

algB.note(:kinematic_fit_chi2,
          "The published analysis applies chi2 < 35 (a loose criterion to exclude " \
          "kinematically incompatible events); the BOSS first pass keeps chi2_cut 200 " \
          "and the tight cut is applied in the ROOT analysis. When more than two pi0 " \
          "candidates exist, the combination with the smallest chi2_6C is retained, as " \
          "done automatically by kinematic_fit.")
    .note(:background_veto,
          "Backgrounds from eta and omega are suppressed by rejecting events with " \
          "|M(pi+ pi- pi0) - m_eta| < 0.02 GeV/c^2 or " \
          "|M(pi+ pi- pi0) - m_omega| < 0.02 GeV/c^2, where M(pi+ pi- pi0) is the " \
          "combination closest to the nominal eta or omega mass. This three-body veto is " \
          "applied on kinematically-fitted momenta and has no dedicated DSL expression.")
    .note(:background_study,
          "The main background channels found with a 1-billion-event J/psi inclusive " \
          "sample are: (1) eta' -> pi+ pi- eta with eta -> pi0 pi0 pi0, (2) " \
          "eta' -> pi0 pi0 eta with eta -> gamma pi+ pi-, (3) eta' -> gamma omega with " \
          "omega -> pi+ pi- pi0, (4) J/psi -> gamma f2(1270) with " \
          "f2(1270) -> pi+ pi- pi0 pi0, and (5) non-resonant J/psi -> gamma 4pi. None of " \
          "them peaks at the eta' mass.")
    .note(:signal_yield_extraction,
          "The yield is obtained from an extended unbinned maximum likelihood fit to the " \
          "pi+ pi- pi0 pi0 invariant mass, with the signal described by the MC shape " \
          "convoluted with a Gaussian and the backgrounds by the dedicated MC shapes, a " \
          "Breit-Wigner f2(1270) tail and the phase-space shape.")

algB.with_decay_card(decay_card_modeB).apply(selB)
algB.execute_on([jpsi_data, jpsi_incMC, exMC_modeB, exMC_bkgB, exMC_bkg_f2])
