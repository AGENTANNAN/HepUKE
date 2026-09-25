# =============================================================================
# BESIII : Observation of a charged charmoniumlike structure in
#          e+e- -> (D* Dbar*)± pi∓ at sqrt(s) = 4.26 GeV
# arXiv:1308.2760v2
#
# Partial-reconstruction technique:
#   reconstructed  : bachelor pi-, D+ (D*+ -> D+ pi0, D+ -> K- pi+ pi+),
#                    and at least one soft pi0 (D*+ -> D+ pi0 or
#                    D*0bar -> D0bar pi0)
#   untagged recoil: D*0bar (used as the D+ pi- recoil mass)
# Data : 827 pb^-1 at 4.260 GeV
# BOSS part only : dataset preparation + event selection up to partial_rec.
# =============================================================================

### ---------------------------------------------------------------------------
### Datasets
### ---------------------------------------------------------------------------
data_4260  = DatasetManager.real_data.find("703_4260")     # 827 pb^-1 @ 4.260 GeV
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")  # Y(4260) + ISR + QED

### ---------------------------------------------------------------------------
### Decay card : e+e- -> D*+ D*0bar pi-  (three-body, PHSP)
### RecID order (DecayCardResolver.rec_id_list, DFS order):
###   0  psi(4260)        1  D*+            2  anti-D*0        3  pi-  (bachelor)
###   4  D+               5  pi0 (from D*+) 6  anti-D0
###   7  pi0 (from D*0bar)                  8  K-   9  pi+_1  10 pi+_2
###  11  K+              12  pi-
### ===========================================================================
decay_card_phsp = <<~DECAYCARD
  Decay psi(4260)
  1.0000 D*+ anti-D*0 pi-   PHSP;
  Enddecay

  Decay D*+
  1.0000 D+ pi0   PHSP;
  Enddecay

  Decay anti-D*0
  1.0000 anti-D0 pi0   PHSP;
  Enddecay

  Decay D+
  1.0000 K- pi+ pi+   PHSP;
  Enddecay

  Decay anti-D0
  1.0000 K+ pi-   PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma   PHSP;
  Enddecay

  End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_4260_DstarDstarpi_phsp"
  config.related_dataset = data_4260
  config.events          = 100_000
  config.decay_card      = decay_card_phsp
  config.cross_section   = :default   # ISR uses the measured sqrt(s) dependence of sigma(D*+D*0 pi-)
end

### ---------------------------------------------------------------------------
### Event selection (BOSS)
### ---------------------------------------------------------------------------
alg_name = "DstarDstarPi"
algorithm = Algorithm.new(alg_name)
algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 4.260]})

event_selection = Selection.new
event_selection.select_track {
                  cos_theta 0.93    # |cos(theta)| < 0.93 in the MDC
                  Vz        10.0    # PCA to the IP within +/-10 cm along the beam
                  Vr        1.0     # PCA within 1 cm in the plane perpendicular to the beam
                  nTot      ">=4"   # at least four charged tracks
                  nChrp     ">=2"   # two pi+
                  nChrn     ">=2"   # one K-, one pi-
                }
               .pid(method: :probability) {
                  prob_cut 0.001
                  identify :kaon, against: [:pion]    # L(K) > L(pi)
                  identify :pion, against: [:kaon]    # L(pi) > L(K)
                  nkm  ">=1"    # at least one K-
                  npip ">=2"    # at least two pi+
                  npim ">=1"    # at least one pi-
                }
               .select_photon {
                  tdc_emc_start     0
                  tdc_emc_end       14
                  energyThreshold_b 0.025
                  energyThreshold_e 0.050
                  angle_to_track    10.0
                  nGam              ">=2"
                }
               # D+ -> K- pi+ pi+ : vertex fit of the three tracks to a common vertex;
               # the fit quality requirement is imposed by minimising the vertex fit chi2
               .secondary_vertex_fit([:km, :pip, :pip]) {
                  build_virtual_particle(:D_plus).by_minimizing_verfit_chi2
                  remove_used_particle_from_candidate_list
                }
               # soft pi0 -> gamma gamma (mass window 0.120 - 0.145 GeV/c^2)
               .kalman_kinematic_fit([:gamma, :gamma]) {
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                  chi2_cut 25
                  npi0 ">=1"   # at least one soft pi0 in the final state
                }
               # partial reconstruction: tag the bachelor pi- and the D+,
               # the D*0bar is inferred from the recoil four-momentum.
               # Recoil mass window = signal region of RM(D+ pi-) (Fig. 3(a))
               .partial_rec([3, 4]) {
                  best_combination_by_mass :D_plus, 1.86962   # M(K- pi+ pi+) closest to m(D+)
                  require_recoil_mass 2.135, 2.175            # RM(D+ pi-) signal region
                }

algorithm.note(:dstar_plus_mass_requirement,
               "The pi0 from D*+ -> D+ pi0 is selected with 2.008 < M(D+ pi0) - M(D+) + m(D+) - M(pi0) + m(pi0) < 2.013 GeV/c^2; the pi0 from D*0bar -> D0bar pi0 is selected by its momentum in the D+ pi- recoil system, P*(pi0) in (0.03, 0.05) GeV/c. Events with at least one pi0 satisfying either requirement are retained.")
       .note(:dstar0_recoil_window,
               "Two-body backgrounds e+e- -> D(*)Dbar(*) are reduced by RM(D+) + M(D+) - m(D+) > 2.3 GeV/c^2; events with RM(pi-) > 4.1 GeV/c^2 (soft pi- from D*- decays) are excluded.")
       .note(:dplus_mass_window,
               "All K- pi+ pi+ combinations with invariant mass in (1.854, 1.884) GeV/c^2 are identified as D+ candidates; the vertex fit quality requirement suppresses non-D+ decays.")
       .note(:resonant_signal_mc,
               "Besides the PHSP D*+ D*0bar pi- sample, the resonant process e+e- -> Zc+(4025) pi- -> D*+ D*0bar pi- is simulated as a cascade with J^P = 1+ and angular distributions from the corresponding matrix element (not expressible as a PHSP EvtGen card).")
       .note(:dstarstar_background,
               "Backgrounds from e+e- -> D** Dbar(*), D** -> D(*) pi(pi) (D0*(2400), D1(2420), D1(2430), D2*(2460)) are simulated separately and subtracted in the fit; they cannot be described by the signal decay card.")
       .note(:combinatorial_background,
               "Combinatorial background is estimated from wrong-sign (WS) events (a reconstructed D+ combined with a pion of the wrong charge), scaled by 1.9 and modelled with a kernel estimate in the RM(pi-) fit.")

algorithm.with_decay_card(decay_card_phsp).apply(event_selection)
algorithm.execute_on([data_4260, incMC_4260, exMC_signal])
