# Analysis: Measurement of the branching fraction for psi(3686) -> omega K+ K-
# with omega -> pi+ pi- pi0, pi0 -> gamma gamma
# (BESIII, 1.06e8 psi(3686) events; continuum subtracted with 44.49 pb-1 at 3.65 GeV)

### Datasets ###
psip_data    = DatasetManager.real_data.find("709_3686")     # psi(3686) data at 3.686 GeV
psip_incMC   = DatasetManager.inclusive_mc.find("709_3686")  # inclusive psi(3686) MC (same size as data)
cont_data    = DatasetManager.real_data.find("709_3650")     # 44.49 pb-1 e+e- data at 3.65 GeV (continuum)
cont_incMC   = DatasetManager.inclusive_mc.find("709_3650")  # continuum inclusive MC

# Decay card for the signal process psi(3686) -> omega K+ K-
# NOTE: the description states the exclusive MC is produced with the data-driven BODY3
# generator (based on EvtGen). BODY3 is not expressible in the EvtGen decay-card syntax,
# so the closest expressible card is given here (omega -> pi+ pi- pi0 via OMEGA_DALITZ)
# and the BODY3 treatment is recorded below with algorithm.note(:efficiency_curve, ...).
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000  omega K+ K-                       PHSP;
    Enddecay

    Decay omega
    1.0000  pi+ pi- pi0                       OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.0000  gamma gamma                       PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for the peaking background psi(3686) -> gamma eta_c(2S), eta_c(2S) -> omega K+ K-
# (studied with a dedicated exclusive MC sample per the description)
decay_card_etac2S = <<~DECAYCARD
    Decay psi(2S)
    1.0000  gamma eta_c(2S)                   PHSP;
    Enddecay

    Decay eta_c(2S)
    1.0000  omega K+ K-                       PHSP;
    Enddecay

    Decay omega
    1.0000  pi+ pi- pi0                       OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.0000  gamma gamma                       PHSP;
    Enddecay

    End
DECAYCARD

# Signal exclusive MC at the psi(3686) energy point (3.686 GeV)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_omegaKK_signal"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

# Signal exclusive MC at the continuum energy point (3.65 GeV), used for the continuum efficiency
exMC_signal_365 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_omegaKK_signal_3650"
  config.related_dataset = cont_data
  config.events          = 50000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

# Exclusive MC for the peaking background psi(3686) -> gamma eta_c(2S) -> gamma omega K+ K-
exMC_etac2S = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_gamma_etac2S_omegaKK"
  config.related_dataset = psip_data
  config.events          = 50000
  config.decay_card      = decay_card_etac2S
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "PsippOmegaKK"  # psi(3686) -> omega K+ K-, omega -> pi+ pi- pi0
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})

event_selection = Selection.new
event_selection.select_track {          # four good charged tracks, zero net charge
                  cos_theta  0.93       # |cos(theta)| < 0.93 (MDC angular coverage)
                  Vz         10.0       # |Vz| < 10 cm (interaction point in beam direction)
                  Vr         1.0        # |Vr| < 1 cm (transverse plane)
                  nChrp      "==2"      # two positively charged tracks
                  nChrn      "==2"      # two negatively charged tracks
                  nNet       "==0"      # net charge zero
                }
               .select_photon {         # at least two good photons (from pi0 -> gamma gamma)
                  tdc_emc_start     0
                  tdc_emc_end       14    # average EMC hit time in [0, 700] ns
                  angle_to_track    20.0  # angle to any charged track > 20 degrees
                  energyThreshold_b 0.025 # E > 25 MeV in the barrel (|cos(theta)| < 0.8)
                  energyThreshold_e 0.050 # E > 50 MeV in the end-cap (0.86 < |cos(theta)| < 0.92)
                  nGam              ">=2" # at least two photon candidates
                }
               .pid(method: :probability) {  # TOF + dE/dx combined PID probability
                  prob_cut 0.001             # P_i > 0.001
                  identify :kaon, against: [:pion, :proton]  # P_K > 0.001 and P_K > P_pi
                  identify :pion, against: [:kaon, :proton]  # P_pi > 0.001 and P_pi > P_K
                  nkp ">=1"
                  nkm ">=1"
                  npip ">=1"
                  npim ">=1"
                }
               .kalman_kinematic_fit([:gamma, :gamma]) {   # reconstruct pi0 from photon pairs
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                  chi2_cut 25
                  npi0 ">=1"
                }
               .assign({:chrgp => :kp, :chrgn => :km})     # K+ K- hypothesis for the kaon pair
               # Vertex fit (all charged tracks from the IP) plus the 4C energy-momentum
               # conserving kinematic fit. If more than two photon candidates exist the
               # combination with the minimum 4C chi2 is kept, which is the DSL default
               # (smallest-chi2 combination).
               .kinematic_fit([:kp, :km, :pip, :pim, :pi0]) {
                  nominal                    # nominal fit: corrected four-momenta used downstream
                  vertex_fit([0, 1, 2, 3])   # vertex fit assuming all tracks originate from the IP
                  constrain_four_momentum    # 4C energy-momentum conservation
                  invariant_mass_of(:pip, :pim, :pi0).within(0.752, 0.812)  # omega mass window
                  chi2_cut 200               # loose BOSS cut: the paper's 5C chi2 < 90 applied in ROOT
               }
               # pi0 mass-constrained (5C) fit: the two photons already enter through the
               # reconstructed :pi0, so the mass constraint on :pi0 is added here.
               .kinematic_fit([:kp, :km, :pip, :pim, :pi0]) {
                  constrain_four_momentum
                  invariant_mass_of(:gamma, :gamma).within(0.11, 0.15)  # 0.11 < M(gamma gamma) < 0.15 GeV/c2
                  chi2_cut 200   # loose BOSS cut: the paper's 5C chi2 < 90 (FOM-optimised) applied in ROOT
               }

my_Algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# BOSS-side procedures that cannot be expressed with the current DSL capabilities
my_Algorithm
  .note(:efficiency_curve, "Signal efficiency determined with the data-driven BODY3 generator \
    (EvtGen based): two BODY3 generators are applied in sequence, one for the three-body decay \
    psi(3686) -> omega K+ K- and one for omega -> pi+ pi- pi0, using the measured Dalitz plot \
    and the K+/K- angular distributions in the psi(3686) CMS as input. Separate BODY3 inputs are \
    used for the 3.686 GeV and 3.65 GeV data samples; eps = 16.9% (3.686 GeV) and 20.7% (3.65 GeV).")
  .note(:background_veto, "Peaking background psi(3686) -> gamma eta_c(2S), eta_c(2S) -> omega K+ K- \
    (not present in the inclusive MC) studied with a dedicated exclusive MC under the assumption \
    B(eta_c(2S) -> omega K+ K-) = 1e-3; its contribution is ~0.1% of the observed omega K+ K- \
    candidates and is neglected.")
  .note(:continuum_subtraction, "Continuum contribution determined from 44.49 pb-1 of e+e- data at \
    3.65 GeV: B = [N_3.686/eps_3.686 - f_c * N_3.65/eps_3.65] / [B(omega -> pi+ pi- pi0) * \
    B(pi0 -> gamma gamma) * N_psi(3686)], with the scaling factor f_c = 3.677 from the luminosities \
    and continuum hadronic cross sections. Interference between psi(3686) decay and continuum \
    production assumed absent.")

my_Algorithm.execute_on([psip_data, psip_incMC, cont_data, cont_incMC,
                         exMC_signal, exMC_signal_365, exMC_etac2S])
