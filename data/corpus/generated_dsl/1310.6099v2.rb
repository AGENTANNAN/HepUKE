### Dataset description ###
data_3686      = DatasetManager.real_data.find("709_3686")        # psi(2S) real data @ 3.686 GeV
incMC_3686     = DatasetManager.inclusive_mc.find("709_3686")     # matching inclusive MC
continuum_3650 = DatasetManager.real_data.find("709_3650")        # 44 pb^-1 continuum data @ 3.65 GeV

# ------------------------------------------------------------------
# Decay cards (EvtGen format); all sub-decays to p pbar in phase space
# ------------------------------------------------------------------
# psi(2S) -> gamma eta_c(2S), eta_c(2S) -> p pbar
decay_card_etac2s = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma eta_c(2S) PHSP;
    Enddecay
    Decay eta_c(2S)
    1.000 p+ anti-p- PHSP;
    Enddecay
    End
DECAYCARD

# psi(2S) -> gamma chi_c0, chi_c0 -> p pbar
decay_card_chic0 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c0 PHSP;
    Enddecay
    Decay chi_c0
    1.000 p+ anti-p- PHSP;
    Enddecay
    End
DECAYCARD

# psi(2S) -> gamma chi_c1, chi_c1 -> p pbar
decay_card_chic1 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c1 PHSP;
    Enddecay
    Decay chi_c1
    1.000 p+ anti-p- PHSP;
    Enddecay
    End
DECAYCARD

# psi(2S) -> gamma chi_c2, chi_c2 -> p pbar
decay_card_chic2 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c2 PHSP;
    Enddecay
    Decay chi_c2
    1.000 p+ anti-p- PHSP;
    Enddecay
    End
DECAYCARD

# psi(2S) -> pi0 h_c, h_c -> p pbar, pi0 -> gamma gamma
decay_card_hc = <<~DECAYCARD
    Decay psi(2S)
    1.000 pi0 h_c PHSP;
    Enddecay
    Decay h_c
    1.000 p+ anti-p- PHSP;
    Enddecay
    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

# psi(2S) -> pi0 p pbar : irreducible background to the pi0 h_c channel
decay_card_pi0ppbar = <<~DECAYCARD
    Decay psi(2S)
    1.000 pi0 p+ anti-p- PHSP;
    Enddecay
    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

# ------------------------------------------------------------------
# Exclusive MC samples (100k events each)
# ------------------------------------------------------------------
exMC_etac2s = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_psip_gamma_etac2s_ppbar"
    config.related_dataset = data_3686
    config.events          = 100000
    config.decay_card      = decay_card_etac2s
    config.cross_section   = :default
end

exMC_chic0 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_psip_gamma_chic0_ppbar"
    config.related_dataset = data_3686
    config.events          = 100000
    config.decay_card      = decay_card_chic0
    config.cross_section   = :default
end

exMC_chic1 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_psip_gamma_chic1_ppbar"
    config.related_dataset = data_3686
    config.events          = 100000
    config.decay_card      = decay_card_chic1
    config.cross_section   = :default
end

exMC_chic2 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_psip_gamma_chic2_ppbar"
    config.related_dataset = data_3686
    config.events          = 100000
    config.decay_card      = decay_card_chic2
    config.cross_section   = :default
end

exMC_hc = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_psip_pi0_hc_ppbar"
    config.related_dataset = data_3686
    config.events          = 100000
    config.decay_card      = decay_card_hc
    config.cross_section   = :default
end

exMC_pi0ppbar = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_psip_pi0_ppbar"
    config.related_dataset = data_3686
    config.events          = 100000
    config.decay_card      = decay_card_pi0ppbar
    config.cross_section   = :default
end

# ==================================================================
# Channel 1 : psi(2S) -> gamma eta_c(2S) / gamma chi_cJ (J=0,1,2),
#             all -> gamma p pbar  (identical final state -> one Algorithm)
# ==================================================================
alg_gamma_name = "GammaPPbar"
alg_gamma = Algorithm.new(alg_gamma_name)
alg_gamma.set_header(["#{alg_gamma_name}Alg/#{alg_gamma_name}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})
         .note(:floating_photon_momentum_3c_fit,
               "for gamma p pbar signal-yield extraction a second 3C kinematic fit is " \
               "performed in which the photon momentum is left floating; expressed as an " \
               "additional 3-momentum-constrained fit whose chi2 is stored for the ROOT-level yield fit")

event_selection_gamma = Selection.new
    .select_track {
        cos_theta 0.93                                   # |cos(theta)| < 0.93
        Vz        10.0                                   # |Vz| < 10 cm
        Vr        1.0                                    # Vr < 1 cm
        nChrp     "==1"                                  # exactly one positive track
        nChrn     "==1"                                  # exactly one negative track
        nNet      "==0"                                  # net charge zero
    }
    .select_photon {
        tdc_emc_start     0                              # EMC timing window 0-14
        tdc_emc_end       14
        energyThreshold_b 0.025                          # > 25 MeV in the barrel
        energyThreshold_e 0.025                          # > 25 MeV in the endcap
        angle_to_track    10.0                           # >= 10 deg from any charged track
        nGam              ">=1"                          # at least one photon
    }
    .pid(method: :probability) {
        prob_cut 0.001                                   # PID probability > 0.001
        identify :proton, against: [:kaon, :pion]        # p and pbar vs K and pi
        nprp     "==1"                                   # exactly one proton
        nprm     "==1"                                   # exactly one anti-proton
    }
    .select_isolated_photon {
        angle_to_prp_track 15.0                          # >= 15 deg from the proton
        angle_to_prm_track 25.0                          # >= 25 deg from the anti-proton
        nGam               ">=1"                         # at least one isolated photon
    }
    .kinematic_fit([:gamma, :prp, :prm]) {               # 4C fit to gamma p pbar
        nominal
        constrain_four_momentum
        chi2_cut 40                                      # chi2 < 40
    }
    .kinematic_fit([:gamma, :prp, :prm]) {               # second 3C fit (floating photon momentum)
        constrain_three_momentum
    }

alg_gamma.with_decay_card(decay_card_etac2s).apply(event_selection_gamma)
root_files_gamma = alg_gamma.execute_on([data_3686, incMC_3686, continuum_3650,
                                         exMC_etac2s, exMC_chic0, exMC_chic1, exMC_chic2])

# ==================================================================
# Channel 2 : psi(2S) -> pi0 h_c -> gamma gamma p pbar
# ==================================================================
alg_pi0hc_name = "Pi0HC"
alg_pi0hc = Algorithm.new(alg_pi0hc_name)
alg_pi0hc.set_header(["#{alg_pi0hc_name}Alg/#{alg_pi0hc_name}.h"])
        .set_constant({"ECMS" => [:double, 3.686]})
        .note(:background_veto,
              "events with M(p pbar gamma_high) < 3.66 GeV/c^2 vetoed to remove the " \
              "gamma chi_cJ (J=1,2) background, where gamma_high is the higher-energy of the two photons")

event_selection_pi0hc = Selection.new
    .select_track {
        cos_theta 0.93                                   # same track criteria
        Vz        10.0
        Vr        1.0
        nChrp     "==1"
        nChrn     "==1"
        nNet      "==0"
    }
    .select_photon {
        tdc_emc_start     0                              # same photon criteria
        tdc_emc_end       14
        energyThreshold_b 0.025
        energyThreshold_e 0.025
        angle_to_track    10.0
        nGam              ">=2"                          # but at least two photons
    }
    .pid(method: :probability) {
        prob_cut 0.001                                   # same PID
        identify :proton, against: [:kaon, :pion]
        nprp     ">=1"                                   # at least one proton
    }
    .select_isolated_photon {
        angle_to_prp_track 15.0
        angle_to_prm_track 25.0
        nGam               ">=2"
    }
    .kinematic_fit([:gamma, :gamma, :prp, :prm]) {       # 4C fit to gamma gamma p pbar
        nominal
        constrain_four_momentum
        invariant_mass_of(:gamma, :gamma).within(0.11, 0.15)  # 0.11 < M(gg) < 0.15 GeV/c^2
        chi2_cut 40                                      # chi2 < 40
    }

alg_pi0hc.with_decay_card(decay_card_hc).apply(event_selection_pi0hc)
root_files_pi0hc = alg_pi0hc.execute_on([data_3686, incMC_3686, continuum_3650,
                                         exMC_hc, exMC_pi0ppbar])