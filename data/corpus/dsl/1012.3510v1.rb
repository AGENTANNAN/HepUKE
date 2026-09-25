# Dataset preparation
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# ============================================================
# Mode 1: J/psi -> gamma pi+ pi- eta', eta' -> gamma rho (rho -> pi+ pi-)
# ============================================================
decay_card_gpipietap_grho = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma X(1835)                    PHSP;
  Enddecay

  Decay X(1835)
  1.0000 pi+ pi- eta'                     PHSP;
  Enddecay

  Decay eta'
  1.0000 gamma rho0                       PHSP;
  Enddecay

  Decay rho0
  1.0000 pi+ pi-                          VSS;
  Enddecay

  End
DECAYCARD

exMC_grho = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_gamma_pipietap_grho"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_gpipietap_grho
  config.cross_section   = :default
end

alg_grho = Algorithm.new("JpsiGammaPiPiEtapGRho")
alg_grho.set_header(["JpsiGammaPiPiEtapGRhoAlg/JpsiGammaPiPiEtapGRho.h"])
        .set_constant({"ECMS" => [:double, 3.097]})

sel_grho = Selection.new
sel_grho.select_track {
           cos_theta 0.93
           Vz        20.0
           Vr        2.0
           nChrp     "==2"
           nChrn     "==2"
           nNet      "==0"
         }
        .select_photon {
           tdc_emc_start     0
           tdc_emc_end       14
           angle_to_track    5.0
           energyThreshold_b 0.100
           energyThreshold_e 0.100
           nGam              ">=2"
         }
        .pid(method: :probability) {
           prob_cut 0.001
           identify :pion, against: [:kaon, :proton]
           npip ">=2"
           npim ">=1"  # >=3 charged tracks identified as pions in total
         }
        .kinematic_fit([:gamma, :gamma, :pip, :pip, :pim, :pim]) {
           nominal
           constrain_four_momentum
           chi2_cut 40
         }

alg_grho
  .note(:min_pion_id_count,
        "At least 3 out of 4 charged tracks required to be identified as " \
        "pions.")
  .note(:pi0_veto,
        "Veto |M(gamma gamma) - m_pi0| < 0.04 GeV/c^2 to suppress pi0 pi+pi- " \
        "pi+ pi- background.")
  .note(:eta_veto_gg,
        "Veto |M(gamma gamma) - m_eta| < 0.03 GeV/c^2 to suppress " \
        "eta pi+ pi- pi+ pi- background.")
  .note(:omega_veto,
        "Veto 0.72 GeV/c^2 < M(gamma gamma) < 0.82 GeV/c^2 to suppress " \
        "omega(-> gamma pi0) pi+ pi- pi+ pi- background.")
  .note(:eta_gpipi_veto,
        "Veto |M(gamma pi+ pi-) - m_eta| < 0.007 GeV/c^2 to suppress " \
        "gamma pi+ pi- eta (eta -> gamma pi+ pi-) background.")
  .note(:rho_selection,
        "Select rho candidate: |M(pi+ pi-) - m_rho| < 0.2 GeV/c^2 (ROOT).")
  .note(:etap_selection,
        "Select eta' candidate: |M(gamma pi+ pi-) - m_etap| < 0.015 GeV/c^2 " \
        "and choose combination closest to m_etap if ambiguous (ROOT).")

alg_grho.with_decay_card(decay_card_gpipietap_grho).apply(sel_grho)
alg_grho.execute_on([jpsi_data, jpsi_incMC, exMC_grho])

# ============================================================
# Mode 2: J/psi -> gamma pi+ pi- eta', eta' -> pi+ pi- eta, eta -> gamma gamma
# ============================================================
decay_card_pipieta = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma X(1835)                    PHSP;
  Enddecay

  Decay X(1835)
  1.0000 pi+ pi- eta'                     PHSP;
  Enddecay

  Decay eta'
  1.0000 pi+ pi- eta                      ETAPRIME_DALITZ;
  Enddecay

  Decay eta
  1.0000 gamma gamma                      PHSP;
  Enddecay

  End
DECAYCARD

exMC_pipieta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_gamma_pipietap_pipieta"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_pipieta
  config.cross_section   = :default
end

alg_pipieta = Algorithm.new("JpsiGammaPiPiEtapPiPiEta")
alg_pipieta.set_header(["JpsiGammaPiPiEtapPiPiEtaAlg/JpsiGammaPiPiEtapPiPiEta.h"])
           .set_constant({"ECMS" => [:double, 3.097]})

sel_pipieta = Selection.new
sel_pipieta.select_track {
              cos_theta 0.93
              Vz        20.0
              Vr        2.0
              nChrp     "==2"
              nChrn     "==2"
              nNet      "==0"
            }
           .select_photon {
              tdc_emc_start     0
              tdc_emc_end       14
              angle_to_track    5.0
              energyThreshold_b 0.100
              energyThreshold_e 0.100
              nGam              ">=3"
            }
           .pid(method: :probability) {
              prob_cut 0.001
              identify :pion, against: [:kaon, :proton]
              npip ">=2"
              npim ">=1"
            }
           .kinematic_fit([:gamma, :gamma, :gamma, :pip, :pip, :pim, :pim]) {
              nominal
              constrain_four_momentum
              chi2_cut 40
            }
           .kalman_kinematic_fit([:gamma, :gamma]) {
              invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
              chi2_cut 40      # paper's 5C chi^2 < 40
              neta ">=1"
            }

alg_pipieta
  .note(:pi0_pairs_veto,
        "Veto |M(gamma gamma) - m_pi0| < 0.04 GeV/c^2 for all photon pairs " \
        "to suppress pi0 combinatorial background.")
  .note(:eta_preselection,
        "Preselect eta candidate by |M(gamma gamma) - m_eta| < 0.03 GeV/c^2 " \
        "before performing the 5C mass-constrained fit.")
  .note(:etap_selection,
        "Select eta' candidate: |M(pi+ pi- eta) - m_etap| < 0.01 GeV/c^2, " \
        "choose combination closest to m_etap (ROOT).")

alg_pipieta.with_decay_card(decay_card_pipieta).apply(sel_pipieta)
alg_pipieta.execute_on([jpsi_data, jpsi_incMC, exMC_pipieta])
