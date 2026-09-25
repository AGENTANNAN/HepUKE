### Dataset preparation ###
# J/psi -> gamma pi0 eta' analysis at sqrt(s)=3.097 GeV
# 1310.6 million J/psi events under BOSS 708

data_Jpsi = DatasetManager.real_data.find("708_3097")
incMC_Jpsi = DatasetManager.inclusive_mc.find("708_3097")

# Decay card: J/psi -> gamma pi0 eta', eta' -> pi+ pi- eta, eta -> gamma gamma, pi0 -> gamma gamma
decay_card_gpi0etap = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma pi0 eta'      PHSP;
  Enddecay
  Decay eta'
  1.0000 pi+ pi- eta          PHSP;
  Enddecay
  Decay eta
  1.0000 gamma gamma           PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma           PHSP;
  Enddecay
  End
DECAYCARD

# Decay card for omega eta' (omega -> gamma pi0) - used for omega eta' BF measurement
decay_card_omega_etap = <<~DECAYCARD
  Decay J/psi
  1.0000 omega eta'           VSS;
  Enddecay
  Decay omega
  1.0000 gamma pi0             VSP_PWAVE;
  Enddecay
  Decay eta'
  1.0000 pi+ pi- eta           PHSP;
  Enddecay
  Decay eta
  1.0000 gamma gamma            PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma            PHSP;
  Enddecay
  End
DECAYCARD

# Decay card for J/psi -> gamma eta_c (eta_c -> pi0 eta') - CPV search
decay_card_gammac = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta_c           VSP_PWAVE;
  Enddecay
  Decay eta_c
  1.0000 pi0 eta'              PHSP;
  Enddecay
  Decay eta'
  1.0000 pi+ pi- eta           PHSP;
  Enddecay
  Decay eta
  1.0000 gamma gamma            PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma            PHSP;
  Enddecay
  End
DECAYCARD

# Exclusive MC for gamma pi0 eta' (phase space)
exMC_gpi0etap = DatasetManager.create_exclusive_mc do |config|
  config.sample_name   = "sig_Jpsi_gpi0etap"
  config.events        = 200_000
  config.decay_card    = decay_card_gpi0etap
  config.related_dataset = data_Jpsi
  config.cross_section = :default
end

# Exclusive MC for omega eta'
exMC_omega_etap = DatasetManager.create_exclusive_mc do |config|
  config.sample_name   = "sig_Jpsi_omega_etap"
  config.events        = 200_000
  config.decay_card    = decay_card_omega_etap
  config.related_dataset = data_Jpsi
  config.cross_section = :default
end

# Exclusive MC for gamma eta_c -> pi0 eta'
exMC_gammac = DatasetManager.create_exclusive_mc do |config|
  config.sample_name   = "sig_Jpsi_gammac_pi0etap"
  config.events        = 200_000
  config.decay_card    = decay_card_gammac
  config.related_dataset = data_Jpsi
  config.cross_section = :default
end

### Algorithm 1: J/psi -> gamma pi0 eta' (common selection for all three sub-analyses) ###
alg_gpi0etap = Algorithm.new("Jpsi_gpi0etap")
alg_gpi0etap.set_header(["Jpsi_gpi0etapAlg/Jpsi_gpi0etap.h"])
             .set_constant({"ECMS" => [:double, 3.097]})

sel_gpi0etap = Selection.new
sel_gpi0etap.select_track {
               cos_theta 0.93
               Vz 10.0
               Vr 1.0
               nChrp ">=1"
               nChrn ">=1"
               nNet "==0"
             }
             .select_photon {
               tdc_emc_start 0
               tdc_emc_end 14
               energyThreshold_b 0.025
               energyThreshold_e 0.050
               angle_to_track 10.0
               nGam ">=5"
             }
             .pid(method: :probability) {
               prob_cut 0.001
               identify :pion, against: [:kaon]
               npip ">=1"
               npim ">=1"
             }
             # 4C kinematic fit: pi+ pi- 5gamma
             .kinematic_fit([:pip, :pim, :gamma, :gamma, :gamma, :gamma, :gamma]) {
               nominal
               constrain_four_momentum
               chi2_cut 200
             }

# Competing hypothesis: pi+ pi- 6gamma (to suppress 6-gamma backgrounds)
sel_gpi0etap
  .kinematic_fit([:pip, :pim, :gamma, :gamma, :gamma, :gamma, :gamma, :gamma]) {
    constrain_four_momentum
  }

# Reconstruct pi0 and eta from photon pairs via Kalman fit (mass-constrained)
sel_gpi0etap
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200
    npi0 ">=1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 200
    neta ">=1"
  }
  # 5C kinematic fit: pi+ pi- gamma gamma gamma eta (eta mass constrained)
  .kinematic_fit([:pip, :pim, :gamma, :gamma, :gamma, :eta]) {
    constrain_four_momentum
    chi2_cut 30
  }

alg_gpi0etap
  .note(:mass_windows, "pi0 mass window: |M(gamma gamma)-m_pi0| < 15 MeV; eta' mass window: |M(pi+pi-eta)-m_eta'| < 15 MeV; pi0 veto: M(radiative_gamma + photon_from_eta) outside [0.115, 0.155] GeV; omega region (0.75-0.82 GeV) excluded in U' search.")
  .note(:chi2_selection, "4C fit chi2(pi+pi-5gamma) required < chi2(pi+pi-6gamma) to suppress 6-gamma backgrounds; pi0/eta photon pair assignment via chi2_pi0eta minimization; 5C fit chi2 < 30 applied.")
  .note(:sub_analyses, "Three sub-analyses using same selection: (1) eta_c -> pi0 eta' CPV search via M(pi0 eta') spectrum; (2) Dark photon U' -> gamma pi0 search via M(gamma pi0) spectrum scanning 0.2-2.1 GeV; (3) omega eta' BF measurement via 2D fit to M(pi+pi-eta) vs M(gamma pi0). All post-fit analyses performed at ROOT level.")
  .note(:dark_photon, "U' search scans M(gamma pi0) in 10 MeV steps from 0.2 to 2.1 GeV excluding omega region [0.75, 0.82] GeV. Signal shape from MC; background modeled with Chebychev polynomial. Mass resolution varies 3.6-10.4 MeV depending on U' mass.")
  .with_decay_card(decay_card_gpi0etap)
  .apply(sel_gpi0etap)

# Execute on J/psi data
alg_gpi0etap.execute_on([data_Jpsi, incMC_Jpsi, exMC_gpi0etap, exMC_omega_etap, exMC_gammac])