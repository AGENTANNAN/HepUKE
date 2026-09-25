# =============================================================================
# BESIII analysis at the J/psi(3097):
#   e+e- -> gamma eta pi0,  eta -> gamma gamma,  pi0 -> gamma gamma
#   (five-photon, all-neutral final state)
# BOSS part: dataset preparation + event selection up to the 4C kinematic fit.
# =============================================================================

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # J/psi(3097) real data (BOSS 7.0.8)
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # matching inclusive MC sample

### Decay cards (EvtGen format) — one per considered process ###

# Mode 1 — J/psi -> gamma eta pi0 (phase space)
decay_card_phsp = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta pi0 PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode 2 — J/psi -> gamma a0(980), a0(980) -> eta pi0
decay_card_a0 = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma a_0(980)0 PHSP;
    Enddecay

    Decay a_0(980)0
    1.0000 eta pi0 PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode 3 — J/psi -> gamma a2(1320), a2(1320) -> eta pi0
decay_card_a2 = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma a_2(1320)0 PHSP;
    Enddecay

    Decay a_2(1320)0
    1.0000 eta pi0 PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode 4 — J/psi -> eta omega, omega -> pi0 gamma (all-neutral chain)
decay_card_eta_omega = <<~DECAYCARD
    Decay J/psi
    1.0000 eta omega PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay omega
    1.0000 pi0 gamma PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode 5 — J/psi -> eta phi, phi -> pi0 gamma (all-neutral chain)
decay_card_eta_phi = <<~DECAYCARD
    Decay J/psi
    1.0000 eta phi PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay phi
    1.0000 pi0 gamma PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode 6 — J/psi -> gamma eta', eta' -> eta pi0 pi0
decay_card_etap_eta2pi0 = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta' PHSP;
    Enddecay

    Decay eta'
    1.0000 eta pi0 pi0 PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode 7 — J/psi -> gamma eta', eta' -> gamma omega, omega -> pi0 gamma
decay_card_etap_gamma_omega = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta' PHSP;
    Enddecay

    Decay eta'
    1.0000 gamma omega PHSP;
    Enddecay

    Decay omega
    1.0000 pi0 gamma PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive MC samples (100k events each) ###
exMC_phsp = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3097_gamma_eta_pi0"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_phsp
  config.cross_section   = :default
end

exMC_a0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3097_gamma_a0_eta_pi0"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_a0
  config.cross_section   = :default
end

exMC_a2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3097_gamma_a2_eta_pi0"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_a2
  config.cross_section   = :default
end

exMC_eta_omega = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3097_eta_omega"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_eta_omega
  config.cross_section   = :default
end

exMC_eta_phi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3097_eta_phi"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_eta_phi
  config.cross_section   = :default
end

exMC_etap_eta2pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3097_gamma_etap_eta_2pi0"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_etap_eta2pi0
  config.cross_section   = :default
end

exMC_etap_gamma_omega = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3097_gamma_etap_gamma_omega"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_etap_gamma_omega
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# Every considered process ends in a purely neutral, five-photon final state,
# therefore no charged track is required and no particle identification is
# applied (npip = npim = 0). One common selection is shared by all channels.
event_selection = Selection.new
    .select_track {                    # charged-track counting: require an all-neutral event
        cos_theta 0.93                 # |cos(theta)| < 0.93
        Vz        10.0                 # |Vz| < 10 cm (beam direction)
        Vr        1.0                  # Vr < 1 cm (transverse plane)
        nChrp     "==0"                # zero positively charged tracks
        nChrn     "==0"                # zero negatively charged tracks
        nNet      "==0"                # net charge zero
    }
    .select_photon {                   # photon (EMC shower) selection
        tdc_emc_start     0            # EMC timing window starts at 0
        tdc_emc_end       14           # ... ends at 14 (unit 50 ns)
        energyThreshold_b 0.025        # barrel energy threshold 25 MeV
        energyThreshold_e 0.050        # endcap energy threshold 50 MeV
        angle_to_track    10.0         # > 10 degrees to the nearest charged track
        nGam              ">=5"        # at least five photons
    }
    # No PID block: the final state is neutral (no identified hadrons/leptons).
    .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :gamma]) {  # 4C fit to the five photons
        nominal                        # nominal fit: its corrected four-momenta are saved
        constrain_four_momentum        # 4C energy-momentum conservation
        chi2_cut 200                   # loose BOSS cut; the published chi2 < 30 is applied later in ROOT
    }

### One Algorithm per decay card (one decay card per Algorithm instance) ###

# Mode 1 — gamma eta pi0
alg_phsp = Algorithm.new("GammaEtaPi0")
alg_phsp.set_header(["GammaEtaPi0Alg/GammaEtaPi0.h"])
        .set_constant({"ECMS" => [:double, 3.097]})
alg_phsp.with_decay_card(decay_card_phsp).apply(event_selection.dup)
root_files_phsp = alg_phsp.execute_on([jpsi_data, jpsi_incMC, exMC_phsp])

# Mode 2 — gamma a0(980) -> gamma eta pi0
alg_a0 = Algorithm.new("GammaA0EtaPi0")
alg_a0.set_header(["GammaA0EtaPi0Alg/GammaA0EtaPi0.h"])
      .set_constant({"ECMS" => [:double, 3.097]})
alg_a0.with_decay_card(decay_card_a0).apply(event_selection.dup)
root_files_a0 = alg_a0.execute_on([jpsi_data, jpsi_incMC, exMC_a0])

# Mode 3 — gamma a2(1320) -> gamma eta pi0
alg_a2 = Algorithm.new("GammaA2EtaPi0")
alg_a2.set_header(["GammaA2EtaPi0Alg/GammaA2EtaPi0.h"])
      .set_constant({"ECMS" => [:double, 3.097]})
alg_a2.with_decay_card(decay_card_a2).apply(event_selection.dup)
root_files_a2 = alg_a2.execute_on([jpsi_data, jpsi_incMC, exMC_a2])

# Mode 4 — eta omega
alg_eta_omega = Algorithm.new("EtaOmega")
alg_eta_omega.set_header(["EtaOmegaAlg/EtaOmega.h"])
             .set_constant({"ECMS" => [:double, 3.097]})
alg_eta_omega.with_decay_card(decay_card_eta_omega).apply(event_selection.dup)
root_files_eta_omega = alg_eta_omega.execute_on([jpsi_data, jpsi_incMC, exMC_eta_omega])

# Mode 5 — eta phi
alg_eta_phi = Algorithm.new("EtaPhi")
alg_eta_phi.set_header(["EtaPhiAlg/EtaPhi.h"])
           .set_constant({"ECMS" => [:double, 3.097]})
alg_eta_phi.with_decay_card(decay_card_eta_phi).apply(event_selection.dup)
root_files_eta_phi = alg_eta_phi.execute_on([jpsi_data, jpsi_incMC, exMC_eta_phi])

# Mode 6 — gamma eta', eta' -> eta pi0 pi0
alg_etap_eta2pi0 = Algorithm.new("GammaEtapEta2Pi0")
alg_etap_eta2pi0.set_header(["GammaEtapEta2Pi0Alg/GammaEtapEta2Pi0.h"])
                .set_constant({"ECMS" => [:double, 3.097]})
alg_etap_eta2pi0.with_decay_card(decay_card_etap_eta2pi0).apply(event_selection.dup)
root_files_etap_eta2pi0 = alg_etap_eta2pi0.execute_on([jpsi_data, jpsi_incMC, exMC_etap_eta2pi0])

# Mode 7 — gamma eta', eta' -> gamma omega
alg_etap_gamma_omega = Algorithm.new("GammaEtapGammaOmega")
alg_etap_gamma_omega.set_header(["GammaEtapGammaOmegaAlg/GammaEtapGammaOmega.h"])
                    .set_constant({"ECMS" => [:double, 3.097]})
alg_etap_gamma_omega.with_decay_card(decay_card_etap_gamma_omega).apply(event_selection.dup)
root_files_etap_gamma_omega = alg_etap_gamma_omega.execute_on([jpsi_data, jpsi_incMC, exMC_etap_gamma_omega])