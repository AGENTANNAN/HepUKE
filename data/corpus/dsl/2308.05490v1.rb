# Paper: 2308.05490v1
# J/psi → phi eta, phi → K+K- (reference) and phi → pi+pi+e-e- (signal, LNV search)
# J/psi data (708_3097), KKMC generator

### Dataset preparation ###
jpsi_data = DatasetManager.load_real_data.find("708_3097")
jpsi_incMC = DatasetManager.load_inclusive_mc.find("708_3097")

# Decay card for reference mode: J/psi → phi eta, phi → K+K-, eta → gamma gamma
decay_card_ref = <<~DECAYCARD
    Decay J/psi
    1.0000 phi eta    HELAMP 1.0 0.0 1.0 0.0 0.0 0.0;
    Enddecay

    Decay phi
    1.0000 K+ K-     VSS;
    Enddecay

    Decay eta
    1.0000 gamma gamma  PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for signal mode: J/psi → phi eta, phi → pi+pi+e-e-, eta → gamma gamma
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 phi eta    HELAMP 1.0 0.0 1.0 0.0 0.0 0.0;
    Enddecay

    Decay phi
    1.0000 pi+ pi+ e- e-   PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma  PHSP;
    Enddecay

    End
DECAYCARD

exMC_ref = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "jpsi_phi_eta_KK_gammagamma"
  config.related_dataset = jpsi_data
  config.events = 1000000
  config.decay_card = decay_card_ref
  config.cross_section = :default
end

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "jpsi_phi_eta_pipiee_gammagamma"
  config.related_dataset = jpsi_data
  config.events = 500000
  config.decay_card = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
# Two independent analysis paths: reference (phi→K+K-) and signal (phi→pi+pi+e-e-)
# Both share J/psi → phi eta with eta → gamma gamma

# ===== Algorithm 1: Reference mode phi → K+K- =====
alg_ref = Algorithm.new("JpsiPhiEta_KK_Ref")
alg_ref.set_header(["JpsiPhiEta_KK_RefAlg/JpsiPhiEta_KK_Ref.h"])
        .set_constant({"ECMS" => [:double, 3.097]})
        .set_alias({"std::vector<double>" => "Vdouble"})

sel_ref = Selection.new
sel_ref.select_track {
             cos_theta   0.93
             Vz   10.0
             Vr   1.0
             nChrp   "==1"     # exactly one positive track
             nChrn   "==1"     # exactly one negative track
             nTot    "==2"     # exactly 2 good charged tracks
             nNet    "==0"     # net charge zero
           }
          .select_photon {
             tdc_emc_start   0
             tdc_emc_end     700
             angle_to_track   20.0   # >20 degrees from nearest track
             energyThreshold_b   0.025  # 25 MeV in barrel
             energyThreshold_e   0.050  # 50 MeV in endcap
             nGam   ">=2"      # At least 2 photons for eta
           }
          .pid(method: :probability) {
             prob_cut   0.001
             identify :kaon, against: [:pion]   # kaon PID: CL_K > 0.001, CL_K > CL_pi
             nkp   "==1"
             nkm   "==1"
           }
          # 4C kinematic fit: e+e- → K+K- gamma gamma
          .kinematic_fit([:kp, :km, :gamma, :gamma]) {
             nominal
             constrain_four_momentum
             chi2_cut 200   # loose cut; optimal in ROOT
             # In ROOT: eta mass cut 0.45 < M_gammagamma < 0.65 GeV
             # phi mass window around 1.0195 GeV
          }

alg_ref.with_decay_card(decay_card_ref).apply(sel_ref)
alg_ref.note(:eta_signal_region, "eta signal region: M_gammagamma in [0.525, 0.565] GeV; sidebands [0.452,0.492] and [0.598,0.638]")
alg_ref.note(:phi_fit, "phi signal extracted from fit to M_KK with MC-convolved double Gaussian signal + inverted ARGUS × polynomial background")
alg_ref.note(:branching_fraction_ref, "B(phi → K+K-) = 49.2%, used to normalize the signal mode")

# ===== Algorithm 2: Signal mode (LNV) phi → pi+pi+e-e- =====
alg_sig = Algorithm.new("JpsiPhiEta_pipiee_LNV")
alg_sig.set_header(["JpsiPhiEta_pipiee_LNVAlg/JpsiPhiEta_pipiee_LNV.h"])
        .set_constant({"ECMS" => [:double, 3.097]})
        .set_alias({"std::vector<double>" => "Vdouble"})

sel_sig = Selection.new
sel_sig.select_track {
             cos_theta   0.93
             Vz   10.0
             Vr   1.0
             nChrp   ">=2"     # at least 2 positive (pi+, e+)
             nChrn   ">=2"     # at least 2 negative (pi-, e-)
             nTot    "==4"     # exactly 4 charged tracks
             nNet    "==0"     # net charge zero
           }
          .select_photon {
             tdc_emc_start   0
             tdc_emc_end     700
             angle_to_track   20.0
             energyThreshold_b   0.025
             energyThreshold_e   0.050
             nGam   ">=2"
           }
          .pid(method: :probability) {
             prob_cut   0.001
             identify :electron, against: [:kaon, :pion]
             identify :pion, against: [:kaon]
             identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                            treat_as_electron_if_energy_above: 0.6
             nep   ">=1"
             nem   ">=1"
             npip  ">=2"   # two pi+
             # In ROOT: CL_e / (CL_e + CL_K + CL_pi) > 0.8 for electrons
             # In ROOT: E/p > 0.8 for p_e >= 0.5 GeV for electrons
           }
          # 4C kinematic fit: e+e- → pi+pi+e-e- gamma gamma
          .kinematic_fit([:pip, :pip, :ep, :ep, :gamma, :gamma]) {
             nominal
             constrain_four_momentum
             chi2_cut 30   # optimized via Punzi significance
           }
          # Competing hypothesis veto: KK pipi gamma gamma (mis-ID kaons→pions)
          .kinematic_fit([:kp, :km, :pip, :pim, :gamma, :gamma]) {
             constrain_four_momentum
             chi2_cut 200
             # In ROOT: require the pipiee hypothesis gives minimum chi2
           }
          # Competing hypothesis veto: pipipipi gamma gamma
          .kinematic_fit([:pip, :pip, :pim, :pim, :gamma, :gamma]) {
             constrain_four_momentum
             chi2_cut 200
           }
          # Competing hypothesis veto: KK ee gamma gamma
          .kinematic_fit([:kp, :km, :ep, :em, :gamma, :gamma]) {
             constrain_four_momentum
             chi2_cut 200
           }
          # Competing hypothesis veto: KK e pi gamma gamma
          .kinematic_fit([:kp, :km, :ep, :pip, :gamma, :gamma]) {
             constrain_four_momentum
             chi2_cut 200
           }
          # Competing hypothesis veto: KK pi e e gamma gamma
          .kinematic_fit([:kp, :km, :pip, :ep, :ep, :gamma, :gamma]) {
             constrain_four_momentum
             chi2_cut 200
           }
          # Competing hypothesis veto: KK pi e pi gamma gamma
          .kinematic_fit([:kp, :km, :pip, :ep, :pip, :gamma, :gamma]) {
             constrain_four_momentum
             chi2_cut 200
           }

alg_sig.with_decay_card(decay_card_signal).apply(sel_sig)
alg_sig.note(:gamma_conversion_veto, "Opening angle theta_pie > 8 degrees between any pion and electron to suppress gamma conversion")
alg_sig.note(:ep_requirement, "Electrons: CL_e > 0.001, CL_e/(CL_e+CL_K+CL_pi) > 0.8; additional E/p > 0.8 for p_e >= 0.5 GeV/c")
alg_sig.note(:signal_window, "Signal region: M_pipiee in [0.99, 1.04] GeV and M_gammagamma in [0.52, 0.57] GeV")
alg_sig.note(:competing_hypothesis, "Events must have minimum chi2 among 6 competing 4C hypotheses to be accepted")
alg_sig.note(:upper_limit, "No signal events observed; upper limit B(phi → pi+pi+e-e-) < 9.7e-6 at 90% CL")

root_files_ref = alg_ref.execute_on([jpsi_data, jpsi_incMC, exMC_ref])
root_files_sig = alg_sig.execute_on([jpsi_data, jpsi_incMC, exMC_signal])