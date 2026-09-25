### Dataset preparation ###
jpsi_data   = DatasetManager.real_data.find("708_3097")
jpsi_incMC  = DatasetManager.inclusive_mc.find("708_3097")
psip_data   = DatasetManager.real_data.find("709_3686")
psip_incMC  = DatasetManager.inclusive_mc.find("709_3686")

# Decay card: J/psi -> Sigma(1385)0 anti-Sigma(1385)0
decay_card_jpsi_s1385 = <<~DECAYCARD
    Decay J/psi
    1.0000 Sigma(1385)0 anti-Sigma(1385)0        PHSP;
    Enddecay

    Decay Sigma(1385)0
    1.0000 Lambda0 pi0                            PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-                                 HypWK;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                            PHSP;
    Enddecay

    End
DECAYCARD

# Decay card: J/psi -> Xi0 anti-Xi0
decay_card_jpsi_xi0 = <<~DECAYCARD
    Decay J/psi
    1.0000 Xi0 anti-Xi0                           PHSP;
    Enddecay

    Decay Xi0
    1.0000 Lambda0 pi0                            PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-                                 HypWK;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                            PHSP;
    Enddecay

    End
DECAYCARD

# Decay card: psi(3686) -> Sigma(1385)0 anti-Sigma(1385)0
decay_card_psip_s1385 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 Sigma(1385)0 anti-Sigma(1385)0        PHSP;
    Enddecay

    Decay Sigma(1385)0
    1.0000 Lambda0 pi0                            PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-                                 HypWK;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                            PHSP;
    Enddecay

    End
DECAYCARD

# Decay card: psi(3686) -> Xi0 anti-Xi0
decay_card_psip_xi0 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 Xi0 anti-Xi0                           PHSP;
    Enddecay

    Decay Xi0
    1.0000 Lambda0 pi0                            PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-                                 HypWK;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                            PHSP;
    Enddecay

    End
DECAYCARD

exMC_jpsi_s1385 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "Jpsi_Sigma1385_exclusive_mc"
  config.related_dataset = jpsi_data
  config.events          = 1000000
  config.decay_card      = decay_card_jpsi_s1385
  config.cross_section   = :default
end

exMC_jpsi_xi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "Jpsi_Xi0_exclusive_mc"
  config.related_dataset = jpsi_data
  config.events          = 1000000
  config.decay_card      = decay_card_jpsi_xi0
  config.cross_section   = :default
end

exMC_psip_s1385 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_Sigma1385_exclusive_mc"
  config.related_dataset = psip_data
  config.events          = 1000000
  config.decay_card      = decay_card_psip_s1385
  config.cross_section   = :default
end

exMC_psip_xi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_Xi0_exclusive_mc"
  config.related_dataset = psip_data
  config.events          = 1000000
  config.decay_card      = decay_card_psip_xi0
  config.cross_section   = :default
end

########################################################################
# Common selection template: single baryon tag via pi0 Lambda
# reconstruct Lambda -> p pi-, pi0 -> gamma gamma, then Sigma(1385)0 or Xi0
# -> pi0 Lambda; the anti-baryon inferred by recoil mass.
########################################################################
def build_selection
  sel = Selection.new
  sel.select_track {
        cos_theta 0.93
        Vz        20.0
        Vr        20.0
        nChrp     ">=1"
        nChrn     ">=1"
      }
     .select_photon {
        tdc_emc_start     0
        tdc_emc_end       14
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam              ">=2"
     }
     .pid(method: :probability) {
        prob_cut 0.001
        identify :proton, against: [:kaon, :pion]
        identify :pion,   against: [:kaon, :proton]
        nprp ">=1"
        npim ">=1"
     }
     .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 20
        npi0 ">=1"
     }
     .secondary_vertex_fit([:prp, :pim]) {
        build_virtual_particle(:Lambda).by_minimizing_mass_difference
        decay_length ">0.0"
        remove_used_particle_from_candidate_list
     }
     .partial_rec([:Lambda, :pi0]) {
        constrain_four_momentum
     }
  sel
end

########################################################################
# J/psi -> Sigma(1385)0 anti-Sigma(1385)0
########################################################################
alg_jpsi_s1385 = Algorithm.new("JpsiToSigma1385")
alg_jpsi_s1385.set_header(["JpsiToSigma1385Alg/JpsiToSigma1385.h"])
              .set_constant({"ECMS" => [:double, 3.097]})

alg_jpsi_s1385.note(:lambda_mass_window,
                    "|M(p pi-) - M(Lambda)| < 5 MeV/c^2 required; window from FOM optimization")
              .note(:sigma1385_signal_window,
                    "|M(pi0 Lambda) - M(Sigma(1385)0)| < 34 MeV/c^2 for J/psi (from FOM); tag Sigma(1385)0 by minimizing |M(pi0 Lambda) - M(Sigma(1385)0)|")
              .note(:anti_baryon_recoil,
                    "anti-Sigma(1385)0 identified via recoil mass of pi0 Lambda system within +/-80 MeV/c^2 of nominal mass")
              .with_decay_card(decay_card_jpsi_s1385)
              .apply(build_selection)

alg_jpsi_s1385.execute_on([jpsi_data, jpsi_incMC, exMC_jpsi_s1385])

########################################################################
# J/psi -> Xi0 anti-Xi0
########################################################################
alg_jpsi_xi0 = Algorithm.new("JpsiToXi0")
alg_jpsi_xi0.set_header(["JpsiToXi0Alg/JpsiToXi0.h"])
            .set_constant({"ECMS" => [:double, 3.097]})

alg_jpsi_xi0.note(:lambda_mass_window,
                  "|M(p pi-) - M(Lambda)| < 5 MeV/c^2 required")
            .note(:xi0_signal_window,
                  "|M(pi0 Lambda) - M(Xi0)| < 10 MeV/c^2 for J/psi (from FOM); tag Xi0 by minimizing |M(pi0 Lambda) - M(Xi0)|")
            .note(:anti_baryon_recoil,
                  "anti-Xi0 identified via recoil mass of pi0 Lambda system within +/-50 MeV/c^2 of nominal Xi0 mass")
            .with_decay_card(decay_card_jpsi_xi0)
            .apply(build_selection)

alg_jpsi_xi0.execute_on([jpsi_data, jpsi_incMC, exMC_jpsi_xi0])

########################################################################
# psi(3686) -> Sigma(1385)0 anti-Sigma(1385)0
########################################################################
alg_psip_s1385 = Algorithm.new("PsipToSigma1385")
alg_psip_s1385.set_header(["PsipToSigma1385Alg/PsipToSigma1385.h"])
              .set_constant({"ECMS" => [:double, 3.686]})

alg_psip_s1385.note(:sigma1385_signal_window,
                    "|M(pi0 Lambda) - M(Sigma(1385)0)| < 35 MeV/c^2 for psi(3686); FOM optimized")
              .note(:jpsi_pipi_transition_veto,
                    "|M_recoil(pi+ pi-) - M(J/psi)| > 5 MeV/c^2 to suppress psi(3686) -> pi+ pi- J/psi")
              .note(:jpsi_pi0pi0_transition_veto,
                    "|M_recoil(pi0 pi0) - M(J/psi)| > 15 MeV/c^2 to suppress psi(3686) -> pi0 pi0 J/psi")
              .with_decay_card(decay_card_psip_s1385)
              .apply(build_selection)

alg_psip_s1385.execute_on([psip_data, psip_incMC, exMC_psip_s1385])

########################################################################
# psi(3686) -> Xi0 anti-Xi0
########################################################################
alg_psip_xi0 = Algorithm.new("PsipToXi0")
alg_psip_xi0.set_header(["PsipToXi0Alg/PsipToXi0.h"])
            .set_constant({"ECMS" => [:double, 3.686]})

alg_psip_xi0.note(:xi0_signal_window,
                  "|M(pi0 Lambda) - M(Xi0)| < 11 MeV/c^2 for psi(3686); FOM optimized")
            .note(:jpsi_pipi_transition_veto,
                  "|M_recoil(pi+ pi-) - M(J/psi)| > 5 MeV/c^2 to suppress psi(3686) -> pi+ pi- J/psi")
            .note(:jpsi_pi0pi0_transition_veto,
                  "|M_recoil(pi0 pi0) - M(J/psi)| > 15 MeV/c^2 to suppress psi(3686) -> pi0 pi0 J/psi")
            .with_decay_card(decay_card_psip_xi0)
            .apply(build_selection)

alg_psip_xi0.execute_on([psip_data, psip_incMC, exMC_psip_xi0])
