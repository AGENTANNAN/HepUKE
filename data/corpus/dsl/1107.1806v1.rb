# ============================================================
# J/psi -> omega eta pi+ pi-   (BESIII, 1107.1806v1)
# omega -> pi+ pi- pi0 ,  eta -> gamma gamma
# Study of f1(1285), eta(1405) and X(1870) in the eta pi+ pi- system
# ============================================================

### Dataset description ###
# J/psi data and the corresponding inclusive MC (2 x 10^8 J/psi, PDG decay table)
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Decay card for the signal process J/psi -> omega eta pi+ pi-
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 omega eta pi+ pi- PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for the non-resonant phase-space control sample
# J/psi -> pi+ pi- pi0 eta pi+ pi- (same final state, no omega resonance)
decay_card_phsp = <<~DECAYCARD
    Decay J/psi
    1.0000 pi+ pi- pi0 eta pi+ pi- PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for the peaking background J/psi -> b1(1235) a0(980),
# b1(1235) -> omega pi and a0(980) -> eta pi
decay_card_bkg_b1a0 = <<~DECAYCARD
    Decay J/psi
    1.0000 b1(1235)+ a0(980)- PHSP;
    Enddecay

    Decay b1(1235)+
    1.000 omega pi+ PHSP;
    Enddecay

    Decay a0(980)-
    1.000 eta pi- PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for the signal channel
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_omega_eta_pip_pim_signal"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

# Exclusive MC for the phase-space process used to derive the phase-space
# correction curve (two million events as quoted in the paper)
exMC_phsp = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_pip_pim_pi0_eta_pip_pim_phsp"
  config.related_dataset = jpsi_data
  config.events          = 2000000
  config.decay_card      = decay_card_phsp
  config.cross_section   = :default
end

# Exclusive MC for the b1(1235) a0(980) background channel
exMC_bkg_b1a0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_b1_1235_a0_980_bkg"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_bkg_b1a0
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "JpsiOmegaEtaPiPi"
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})   # J/psi center-of-mass energy (GeV)
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {
                  # Four charged tracks with net charge zero, all assumed to be pions
                  cos_theta   0.93    # |cos(theta)| < 0.93
                  Vz          20.0    # |Vz| < 20 cm from the IP along the beam direction
                  Vr          2.0     # |Vr| < 2 cm from the beam line in the transverse plane
                  nChrp       "==2"
                  nChrn       "==2"
                  nNet        "==0"
                }
               .select_photon {
                  # EMC showers: barrel E > 25 MeV, endcap 50 MeV, >= 10 deg from nearest track
                  tdc_emc_start      0
                  tdc_emc_end        14
                  angle_to_track     10.0
                  energyThreshold_b  0.025
                  energyThreshold_e  0.050
                  nGam               ">=4"   # two photons for eta and two for pi0
                }
               .assign({chrgp: :pip, chrgn: :pim})   # all tracks treated as pions (no PID required)
               .kalman_kinematic_fit([:gamma, :gamma]) {
                  # Reconstruct pi0 (from omega -> pi+ pi- pi0) as a gamma gamma pair
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                  chi2_cut 25
                  npi0 ">=1"
                }
               .kalman_kinematic_fit([:gamma, :gamma]) {
                  # Reconstruct eta (from the eta pi+ pi- system) as a gamma gamma pair
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
                  chi2_cut 25
                  neta ">=1"
                }
                # Nominal 4C kinematic fit: constrain the total four-momentum of all
                # final-state particles to the J/psi mass and the initial e+e- three-momentum.
                # The four charged tracks are also required to share a common vertex.
               .kinematic_fit([:pip, :pim, :pip, :pim, :pi0, :eta]) {
                  nominal
                  vertex_fit([0, 1, 2, 3])   # IP-constrained vertex fit on the four charged tracks
                  constrain_four_momentum
                  chi2_cut 50                # paper quotes chi2_4C < 50
                }

my_Algorithm
  .note(:vertex_fit_chi2, "the vertex fit of the four final-state tracks is required to satisfy chi2_V < 100; the vertex constraint is applied inside the nominal 4C kinematic fit via vertex_fit([0,1,2,3]) and the chi2_V < 100 quality requirement cannot be expressed directly in the DSL")
  .note(:eta_pi0_mass_windows, "eta (524, 572) MeV/c^2 and pi0 (122, 148) MeV/c^2 gamma-gamma invariant-mass windows are applied on the four-momenta returned by the kinematic fit, hence they are post-fit selections handled at the ROOT level")
  .note(:omega_mass_window, "the omega candidate is chosen among the selected charged pion pairs by minimizing |M(pi+ pi- pi0) - m_omega|, and required to be less than 28 MeV/c^2; this uses the kinematic-fit-corrected four-momenta and is applied at the ROOT level")
  .note(:best_candidate_selection, "when several combinations satisfy the requirements, the eta pi0 4pi combination with minimum chi2_4C is retained; if more than one four-photon combination falls in the eta/pi0 mass ranges the assignment with the lowest chi2_{eta pi0} = sqrt(P_eta^2 + P_pi0^2) is used, with pulls P = (M_{gamma gamma} - m_{eta/pi0}) / sigma_{gamma gamma}")
  .note(:background_sidebands, "non-omega and/or non-a0(980) background is estimated from the two-dimensional mass sidebands of omega and a0(980), weighted by the horizontal (0.48), vertical (1.58) and diagonal (0.76) factors obtained from a two-dimensional fit to M_omega(pi+ pi- pi0) versus M_a0(980)(eta pi); this is a ROOT-level fit procedure")
  .note(:peaking_background_veto, "the peaking background J/psi -> b1(1235) a0(980) (b1(1235) -> omega pi, a0(980) -> eta pi) is estimated from a dedicated exclusive MC sample and from a two-dimensional fit to M(omega pi) versus M(eta pi)")

my_Algorithm.with_decay_card(decay_card_signal).apply(event_selection)

root_files = my_Algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal, exMC_phsp, exMC_bkg_b1a0])
