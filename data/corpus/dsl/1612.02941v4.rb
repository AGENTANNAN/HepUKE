# ============================================================================
# BESIII: improved measurement of Br(eta_c -> phi phi) and search for
# eta_c -> omega phi
# arXiv:1612.02941v4,  (223.7 +- 1.4) x 10^6 J/psi events at sqrt(s) = 3.097 GeV
#
# Signal chain:
#   Mode I  : J/psi -> gamma eta_c, eta_c -> phi phi, phi -> K+ K-
#             final state gamma 2(K+ K-)
#   Mode II : J/psi -> gamma eta_c, eta_c -> omega phi,
#             omega -> pi+ pi- pi0, phi -> K+ K-, pi0 -> gamma gamma
#             final state 3 gamma K+ K- pi+ pi-
# The two modes have different final states, photon multiplicities, PID
# requirements and kinematic-fit hypotheses -> two Algorithm objects (Rule T1).
#
# The eta_c yield is obtained from a ROOT-level amplitude analysis of the
# M(phi phi) / M(omega phi) spectra (helicity-covariant amplitudes, unbinned
# maximum likelihood fit with MINUIT, Bayesian upper limit for omega phi);
# that part is outside the BOSS scope.
# ============================================================================

### ---------------------------- Dataset preparation ---------------------------- ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # J/psi (3.097 GeV) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # inclusive J/psi MC (background study)

# ---------------------------------------------------------------------------
# Mode I : J/psi -> gamma eta_c, eta_c -> phi phi -> gamma 2(K+ K-)
# ---------------------------------------------------------------------------
decay_card_phiphi = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta_c           PHSP;
  Enddecay

  Decay eta_c
  1.0000 phi phi               PHSP;
  Enddecay

  Decay phi
  1.0000 K+ K-                 VSS;
  Enddecay

  End
DECAYCARD

exMC_phiphi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_gamma_etac_phiphi"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_phiphi
  config.cross_section   = :default
end
exMC_phiphi.save_to_config(format: :yaml, file_path: 'temp_for_test')

# Peaking / non-peaking backgrounds with the same 2(K+ K-) final state:
# J/psi -> gamma phi K+ K- and J/psi -> gamma K+ K- K+ K- (with or without an
# eta_c intermediate state), and the pi0 2(K+ K-) channels J/psi -> phi f_1
# (f_1 -> K+ K- pi0) and J/psi -> phi K*(892)+- K-+ (K*(892)+- -> K+- pi0).
decay_card_bkg_gphiKK = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma phi K+ K-       PHSP;
  Enddecay

  Decay phi
  1.0000 K+ K-                 VSS;
  Enddecay

  End
DECAYCARD

decay_card_bkg_gKKKK = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma K+ K- K+ K-     PHSP;
  Enddecay

  End
DECAYCARD

exMC_bkg_gphiKK = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_gamma_phi_KK"
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_bkg_gphiKK
  config.cross_section   = :default
end

exMC_bkg_gKKKK = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_gamma_KKKK"
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_bkg_gKKKK
  config.cross_section   = :default
end

# ---------------------------------------------------------------------------
# Mode II : J/psi -> gamma eta_c, eta_c -> omega phi,
#           omega -> pi+ pi- pi0, phi -> K+ K-, pi0 -> gamma gamma
# ---------------------------------------------------------------------------
decay_card_omegaphi = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta_c           PHSP;
  Enddecay

  Decay eta_c
  1.0000 omega phi             PHSP;
  Enddecay

  Decay omega
  1.0000 pi+ pi- pi0           OMEGA_DALITZ;
  Enddecay

  Decay phi
  1.0000 K+ K-                 VSS;
  Enddecay

  Decay pi0
  1.0000 gamma gamma           PHSP;
  Enddecay

  End
DECAYCARD

exMC_omegaphi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_gamma_etac_omegaphi"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_omegaphi
  config.cross_section   = :default
end

# Dominant background for the omega phi search: J/psi -> eta' phi with
# eta' -> gamma omega; a small amount comes from J/psi -> f_0(980) omega ->
# K+ K- omega and J/psi -> f_X omega -> pi0 K+ K- omega (f_X = f_1(1285),
# f_1(1420)).  These have the same 3 gamma K+ K- pi+ pi- final state.
decay_card_bkg_etapPhi = <<~DECAYCARD
  Decay J/psi
  1.0000 eta' phi              PHSP;
  Enddecay

  Decay eta'
  1.0000 gamma omega           PHSP;
  Enddecay

  Decay omega
  1.0000 pi+ pi- pi0           OMEGA_DALITZ;
  Enddecay

  Decay phi
  1.0000 K+ K-                 VSS;
  Enddecay

  Decay pi0
  1.0000 gamma gamma           PHSP;
  Enddecay

  End
DECAYCARD

exMC_bkg_etapPhi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_etap_phi_etap2gammaomega"
  config.related_dataset = jpsi_data
  config.events          = 200_000
  config.decay_card      = decay_card_bkg_etapPhi
  config.cross_section   = :default
end

# ============================================================================
# ALGORITHM I: J/psi -> gamma eta_c, eta_c -> phi phi -> gamma 2(K+ K-)
# Topology: 1 radiative photon + 4 charged tracks (2 K+, 2 K-); no PID.
# ============================================================================
alg_name_phiphi = "JpsiGammaEtacPhiPhi"
alg_phiphi = Algorithm.new(alg_name_phiphi)
alg_phiphi.set_header(["#{alg_name_phiphi}Alg/#{alg_name_phiphi}.h"])
          .set_constant({"ECMS" => [:double, 3.097]})

sel_phiphi = Selection.new
sel_phiphi.select_track {                 # four charged tracks with net charge zero
             cos_theta 0.93               # |cos(theta)| < 0.93
             Vz        10.0               # within +-10 cm of the IP along the beam
             Vr         1.0               # within 1 cm in the plane perpendicular to the beam
             nChrp     "==2"
             nChrn     "==2"
             nNet      "==0"
           }
           .select_photon {               # one radiative photon from J/psi -> gamma eta_c
             tdc_emc_start     0          # EMC timing suppresses electronic noise
             tdc_emc_end       14         # (0, 700) ns
             angle_to_track    10.0       # > 10 degrees from the nearest charged track
             energyThreshold_b 0.025      # barrel (|cos(theta)| < 0.8):  E > 25 MeV
             energyThreshold_e 0.050      # end-cap (0.86 < |cos(theta)| < 0.92): E > 50 MeV
             nGam              ">=1"
           }
           # No PID: the four tracks are direct K+/K- candidates (the paper
           # applies no PID for the 2(K+ K-) channel)
           .assign({:chrgp => :kp, :chrgn => :km})
           # 4C kinematic fit under the J/psi -> gamma 2(K+ K-) hypothesis.
           # The four-momenta of the four tracks and the photon are constrained
           # to the initial e+e- (J/psi) four-momentum.  When more than one
           # photon or K+ K- pairing is possible, the combination with the
           # smallest chi2_4C is retained.  The paper requires chi2_4C < 100
           # (optimised on S/sqrt(S+B)); the loose BOSS default is used here and
           # the tight value is applied in ROOT (Rule T3).
           # 4C kinematic fit under the J/psi -> gamma 2(K+ K-) hypothesis.
           # The four-momenta of the four tracks and the photon are constrained
           # to the initial e+e- (J/psi) four-momentum.  When more than one
           # photon or K+ K- pairing is possible, the combination with the
           # smallest chi2_4C is retained.  No mass window is imposed inside the
           # fit: the two phi candidates are defined on the fit-corrected
           # M(K+ K-) afterwards in ROOT (see the phi_candidate_pairing and
           # phi_mass_window notes).  The paper requires chi2_4C < 100
           # (optimised on S/sqrt(S+B)); the loose BOSS default is used here and
           # the tight value is applied in ROOT (Rule T3).
           .kinematic_fit([:kp, :km, :kp, :km, :gamma]) {
             nominal
             constrain_four_momentum
             chi2_cut 200
           }

alg_phiphi
  .note(:phi_candidate_pairing,
        "The two phi candidates are reconstructed from the selected 2(K+ K-) tracks; only the " \
        "K+ K- pairing that minimises |M(K+K-)^(1) - M_phi|^2 + |M(K+K-)^(2) - M_phi|^2 is " \
        "retained.  No dedicated DSL primitive expresses this two-pair combinatorial criterion " \
        "(nested for_each is experimental); the pairing is fixed at the ROOT level.")
  .note(:phi_mass_window,
        "Both phi candidates must satisfy |M(K+ K-) - M_phi| < 0.02 GeV/c^2 (the window is " \
        "optimised on S/sqrt(S+B)); applied on the kinematic-fit-corrected invariant masses in ROOT.")
  .note(:pid_omission,
        "No PID requirement is applied to the four charged tracks in the eta_c -> phi phi " \
        "channel, as in the paper: any track combination satisfying the track-quality and " \
        "multiplicity requirements enters the 4C fit.")
  .note(:background_veto,
        "Dominant backgrounds with the same 2(K+ K-) final state are J/psi -> gamma phi K+ K- " \
        "and J/psi -> gamma K+ K- K+ K-, with or without an eta_c intermediate state; they are " \
        "peaking and non-peaking backgrounds in the M(2(K+ K-)) distribution with expected " \
        "yields of 26 and 75 events from the inclusive MC.  In addition J/psi -> phi f_1(1420) / " \
        "f_1(1285) (f_1 -> K+ K- pi0) and J/psi -> phi K*(892)+- K-+ (K*(892)+- -> K+- pi0) give " \
        "a pi0 2(K+ K-) final state similar to signal, with a very low detection efficiency " \
        "(< 0.1%) and no peaking in the eta_c signal range; 43 such events are found in data.")
  .note(:amplitude_analysis,
        "The eta_c -> phi phi yield (549 +- 65) is extracted by an amplitude analysis of the 1276 " \
        "selected candidates, assuming J/psi -> gamma phi phi with or without an eta_c " \
        "intermediate state.  Helicity-covariant amplitudes combine the eta_c Breit-Wigner " \
        "(multiplied by an E_gamma damping factor with beta = 0.065 GeV) with nonresonant " \
        "J/psi -> gamma phi phi components of J^P = 0-, 0+ and 2+ in the phi phi system; the " \
        "eta_c mass and width are fixed to M = 2.984 GeV/c^2 and Gamma = 0.032 GeV.  An unbinned " \
        "maximum-likelihood fit is minimised with MINUIT, with the background log-likelihood " \
        "(normalised to 101 peaking + non-peaking events from MC) subtracted.  The mass " \
        "resolution, the non-eta_c components and the fit range contribute to the systematic " \
        "uncertainty.  ROOT-level procedure.")
  .note(:signal_mc_model,
        "The detection efficiency (epsilon = 24%) is determined from an amplitude-weighted MC " \
        "sample generated with the amplitude model and parameters fixed to the fit results, not " \
        "from a phase-space sample; the phase-space card is used only to build the efficiency " \
        "map needed by the amplitude fit.")
  .note(:branching_fraction,
        "Br(J/psi -> gamma eta_c) Br(eta_c -> phi phi) = N_sig / (N_J/psi epsilon Br^2(phi -> " \
        "K+ K-)) = (4.3 +- 0.5 +0.5/-1.2) x 10^-5 with N_J/psi = 223.7 x 10^6; using " \
        "Br(J/psi -> gamma eta_c) = (1.7 +- 0.4)% this gives Br(eta_c -> phi phi) = " \
        "(2.5 +- 0.3 +0.3/-0.7 +- 0.6) x 10^-3.")
  .with_decay_card(decay_card_phiphi)
  .apply(sel_phiphi)

# ============================================================================
# ALGORITHM II: J/psi -> gamma eta_c, eta_c -> omega phi -> 3 gamma K+ K- pi+ pi-
# Topology: 3 photons (2 from pi0 -> gamma gamma + 1 radiative) + 4 charged
# tracks (2 K+ K- from phi, pi+ pi- from omega).
# ============================================================================
alg_name_omegaphi = "JpsiGammaEtacOmegaPhi"
alg_omegaphi = Algorithm.new(alg_name_omegaphi)
alg_omegaphi.set_header(["#{alg_name_omegaphi}Alg/#{alg_name_omegaphi}.h"])
           .set_constant({"ECMS" => [:double, 3.097]})

sel_omegaphi = Selection.new
sel_omegaphi.select_track {                # four charged tracks with net charge zero
              cos_theta 0.93
              Vz        10.0
              Vr         1.0
              nChrp     "==2"
              nChrn     "==2"
              nNet      "==0"
            }
            .select_photon {               # three photons: pi0 -> gamma gamma + radiative photon
              tdc_emc_start     0
              tdc_emc_end       14
              angle_to_track    10.0
              energyThreshold_b 0.025
              energyThreshold_e 0.050
              nGam              ">=3"
            }
            .pid(method: :probability) {   # dE/dx + TOF combined confidence levels
              prob_cut 0.001
              identify :kaon, against: [:pion, :proton]
              nkp "==1"
              nkm "==1"
            }
            .remove([:kp <= :chrgp, :km <= :chrgn])   # keep the two kaons in the kp/km lists
            .assign({:chrgp => :pip, :chrgn => :pim}) # remaining tracks are pi+ and pi-
            # pi0 -> gamma gamma (1C mass constraint); when more than one photon
            # pair is possible the combination with the mass closest to the
            # nominal pi0 mass is chosen (smallest chi2 of the Kalman fit).
            .kalman_kinematic_fit([:gamma, :gamma]) {
              invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
              chi2_cut 25
              npi0 ">=1"
            }
            # 4C kinematic fit under the J/psi -> gamma K+ K- pi+ pi- pi0
            # hypothesis (the pi0 enters through the reconstructed :pi0, so its
            # daughter photons are not listed separately).  All photon
            # combinations with the four charged tracks are tried and the
            # smallest chi2_4C combination is retained.  The paper requires
            # chi2_4C < 40; the loose BOSS default is used here and the tight
            # value is applied in ROOT (Rule T3).
            .kinematic_fit([:gamma, :kp, :km, :pip, :pim, :pi0]) {
              nominal
              constrain_four_momentum
              invariant_mass_of(:kp, :km).within(1.011, 1.027)   # |M(K+K-) - M_phi| < 0.008 GeV/c^2
              invariant_mass_of(:pip, :pim, :pi0).within(0.752, 0.812) # |M(3pi) - M_omega| < 0.03 GeV/c^2
              chi2_cut 200
            }

alg_omegaphi
  .note(:pi0_candidate_selection,
        "The photon pair whose invariant mass is closest to the nominal pi0 mass is chosen as the " \
        "pi0 candidate and |M(gamma gamma) - M_pi0| < 0.02 GeV/c^2 is required; the condition is " \
        "expressed here through the 1C mass-constrained Kalman fit on the photon pair, with the " \
        "minimum-chi2 combination retained by default.  The explicit 0.02 GeV/c^2 window on the " \
        "unconstrained gamma gamma mass is applied in ROOT.")
  .note(:kaon_momentum_window,
        "The two identified kaons are required to lie in the momentum range 0.3-0.9 GeV/c (the " \
        "average kaon PID efficiency over this range is about 8%); this momentum-dependent " \
        "selection is applied in ROOT on the reconstructed track momenta.")
  .note(:phi_omega_mass_windows,
        "The phi and omega candidates are selected with |M(K+ K-) - M_phi| < 0.008 GeV/c^2 and " \
        "|M(pi+ pi- pi0) - M_omega| < 0.03 GeV/c^2 (windows optimised on S/sqrt(S+B)); they are " \
        "applied on the kinematic-fit-corrected invariant masses in ROOT.")
  .note(:background_veto,
        "The dominant background for the eta_c -> omega phi search is J/psi -> eta' phi with " \
        "eta' -> gamma omega; a small contribution comes from J/psi -> f_0(980) omega -> K+ K- " \
        "omega and J/psi -> f_X omega -> pi0 K+ K- omega (f_X = f_1(1285), f_1(1420)).  The total " \
        "background from the inclusive MC is small compared with the number of selected " \
        "candidates and appears flat in M(omega phi).")
  .note(:upper_limit,
        "No significant eta_c signal is observed in the M(omega phi) distribution between 2.70 and " \
        "3.05 GeV/c^2.  The signal yield upper limit at the 90% C.L. (N_up = 18) is obtained with a " \
        "Bayesian method from the normalised likelihood distribution versus the signal yield, with " \
        "the eta_c line shape taken from MC (M and Gamma from the BESIII measurement), the known " \
        "MC background fixed in shape and magnitude and the remaining background described by a " \
        "second-order Chebyshev polynomial.  This gives Br(eta_c -> omega phi) < N_up / (N_J/psi " \
        "epsilon Br (1 - sigma_sys)) = 2.5 x 10^-4 with epsilon = 5.9%, sigma_sys = 25.8% and Br " \
        "the product of Br(J/psi -> gamma eta_c), Br(phi -> K+ K-) and Br(omega -> pi+ pi- pi0).  " \
        "ROOT-level procedure.")
  .note(:signal_mc_model,
        "The omega -> pi+ pi- pi0 decay is generated with the OMEGA_DALITZ model (PDG-driven " \
        "Dalitz distribution); the eta_c -> omega phi decay is generated with a phase-space model.")
  .with_decay_card(decay_card_omegaphi)
  .apply(sel_omegaphi)

### --------------------------------- Execution --------------------------------- ###
alg_phiphi.execute_on([jpsi_data, jpsi_incMC, exMC_phiphi,
                       exMC_bkg_gphiKK, exMC_bkg_gKKKK])

alg_omegaphi.execute_on([jpsi_data, jpsi_incMC, exMC_omegaphi, exMC_bkg_etapPhi])
