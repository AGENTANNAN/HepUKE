### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # 10 billion J/psi events
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # J/psi inclusive MC

# Signal decay card: J/psi -> Sigma0 Sigmabar0, with Sigma0 -> Lambda e+ e- (Dalitz)
# and Sigmabar0 -> Lambdabar gamma (tag side).  Lambda(bar) -> p pi- (charge conj.).
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000  Sigma0  anti-Sigma0                   HELAMP 1.0 0.0 1.0 0.0;
    Enddecay

    Decay Sigma0
    1.0000  Lambda0  e+   e-                      PHSP;
    Enddecay

    Decay anti-Sigma0
    1.0000  anti-Lambda0  gamma                   PHSP;
    Enddecay

    Decay Lambda0
    1.0000  p+   pi-                              HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000  anti-p-   pi+                         HypWK;
    Enddecay

    End
DECAYCARD

# Background MC: peaking background from Sigma0 -> Lambda gamma (double radiative)
decay_card_bkg = <<~DECAYCARD
    Decay J/psi
    1.0000  Sigma0  anti-Sigma0                   HELAMP 1.0 0.0 1.0 0.0;
    Enddecay

    Decay Sigma0
    1.0000  Lambda0  gamma                        PHSP;
    Enddecay

    Decay anti-Sigma0
    1.0000  anti-Lambda0  gamma                   PHSP;
    Enddecay

    Decay Lambda0
    1.0000  p+   pi-                              HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000  anti-p-   pi+                         HypWK;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_Sigma0Sigmabar_LambdaEE"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

exMC_bkg = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_Sigma0Sigmabar_LambdaGamma"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_bkg
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "Sigma0Dalitz"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {
                  cos_theta 0.93         # |cos(theta)| < 0.93
                  Vz        30.0         # |Vz| < 30 cm; no Vr cut (long Lambda decay length)
                  nChrp     ">=3"        # DT: p (from Lambda), pi+ (from Lambdabar), e+
                  nChrn     ">=3"        # DT: pbar, pi-, e-
                }
               .select_photon {
                  energyThreshold_b 0.025 # >25 MeV in barrel  (|cos theta|<0.80)
                  energyThreshold_e 0.050 # >50 MeV in end cap (0.86<|cos theta|<0.92)
                  tdc_emc_start     0
                  tdc_emc_end       14    # shower time in [0, 700] ns
                  angle_to_track    10.0
                  nGam              ">=1" # radiative photon from Sigmabar0 -> Lambdabar gamma
                }
               .pid(method: :probability) {
                  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.5,
                                                 treat_as_electron_if_energy_above: 0.6
                  prob_cut 0.001
                  identify :proton, against: [:kaon, :pion]  # p and pbar
                  nprp "==1"
                  nprm "==1"
                  nlp  "==1"
                  nlm  "==1"
                }
               # Assign the remaining charged tracks (soft pions from Lambda / Lambdabar)
               .remove([:prp <= :chrgp])
               .remove([:prm <= :chrgn])
               .remove([:lp  <= :chrgp])
               .remove([:lm  <= :chrgn])
               .assign({:chrgp => :pip, :chrgn => :pim})
               # Secondary vertex fit: Lambda -> p pi-
               .secondary_vertex_fit([:prp, :pim]) {
                  build_virtual_particle(:Lambda).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
                }
               # Secondary vertex fit: Lambdabar -> pbar pi+
               .secondary_vertex_fit([:prm, :pip]) {
                  build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
                }
               # 4C global kinematic fit on J/psi -> p pbar pi+ pi- gamma e+ e-
               .kinematic_fit([:Lambda, :Lambda_bar, :gamma, :lp, :lm]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 200            # loose cut in BOSS; tight cut applied in ROOT
                }

my_algorithm
  .note(:double_tag_strategy,
        "Absolute BF measurement uses a double-tag (DT) method: single-tag " \
        "reconstructs Sigmabar0 -> Lambdabar gamma with |M(Lambdabar gamma) - " \
        "M(Sigma0)| < 15 MeV/c^2; DT then reconstructs Sigma0 -> Lambda e+ e- " \
        "from remaining tracks recoiling against the ST Sigmabar0.")
  .note(:lambda_selection,
        "Lambda(bar) reconstruction: |M(p pi) - M(Lambda)| < 8 MeV/c^2, " \
        "L_decay > 0, secondary vertex chi^2 < 200. Enforced in ROOT " \
        "post-selection variables.")
  .note(:pbar_pip_momentum_cuts,
        "In the ST branch, negatively-charged proton candidate is required to " \
        "have p > 0.55 GeV/c and positively-charged pion candidate p < 0.35 GeV/c " \
        "based on kinematic separation from MC (Lambdabar side).")
  .note(:electron_pid,
        "Signal e+e- pair: |Vz|<10 cm, |Vr|<1 cm, and PID requirement " \
        "P(e)/(P(e)+P(pi)+P(K)) > 0.8 applied per candidate.")
  .note(:gamma_conversion_veto,
        "Gamma-conversion background from beam pipe (r~31.5 mm) and MDC inner " \
        "wall (r~54 mm) suppressed using an M_ee^BP vs Rxy and Phi_ee vs Rxy " \
        "veto (procedure per Ref. [30] of the paper). Applied post-selection.")
  .note(:mass_recoil_tag,
        "M_rec^corr = sqrt((E_cms - E_Lambdabar - E_gamma)^2 - " \
        "|P_Lambdabar + P_gamma|^2) + M(Lambdabar gamma) - M(Sigma0) used to " \
        "identify J/psi -> Sigma0 Sigmabar0 in the ST fit (ROOT).")

my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal, exMC_bkg])
