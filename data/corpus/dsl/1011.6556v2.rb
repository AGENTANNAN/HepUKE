# Dataset preparation
psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# psi(2S) -> gamma chi_cJ, chi_cJ -> 4 pi0 (J = 0, 1, 2 share the same
# final state and selection; use a single Algorithm object per Rule Special Case)
decay_card_chic0 = <<~DECAYCARD
  Decay psi(2S)
  1.0000 gamma chi_c0        PHSP;
  Enddecay

  Decay chi_c0
  1.0000 pi0 pi0 pi0 pi0    PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma         PHSP;
  Enddecay

  End
DECAYCARD

decay_card_chic1 = <<~DECAYCARD
  Decay psi(2S)
  1.0000 gamma chi_c1        PHSP;
  Enddecay

  Decay chi_c1
  1.0000 pi0 pi0 pi0 pi0    PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma         PHSP;
  Enddecay

  End
DECAYCARD

decay_card_chic2 = <<~DECAYCARD
  Decay psi(2S)
  1.0000 gamma chi_c2        PHSP;
  Enddecay

  Decay chi_c2
  1.0000 pi0 pi0 pi0 pi0    PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma         PHSP;
  Enddecay

  End
DECAYCARD

exMC_chic0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_gamma_chic0_4pi0"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_chic0
  config.cross_section   = :default
end

exMC_chic1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_gamma_chic1_4pi0"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_chic1
  config.cross_section   = :default
end

exMC_chic2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_gamma_chic2_4pi0"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_chic2
  config.cross_section   = :default
end

# Event selection (BOSS): single algorithm because all three chi_cJ modes
# share identical final state (gamma + 4 pi0) and selection criteria.
alg_name = "PsipGammaChicJ4Pi0"
algorithm = Algorithm.new(alg_name)
algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})

event_selection = Selection.new
event_selection.select_track {
                  cos_theta 0.93
                  Vz        10.0
                  Vr        1.0
                  nChrp     "==0"
                  nChrn     "==0"
                }
               .select_photon {
                  tdc_emc_start     0
                  tdc_emc_end       14
                  energyThreshold_b 0.025
                  energyThreshold_e 0.050
                  nGam              "==9"     # exactly 9 photons required
                }
               .kalman_kinematic_fit([:gamma, :gamma]) {
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                  chi2_cut 25
                  npi0 ">=4"
                }
               .kinematic_fit([:gamma, :pi0, :pi0, :pi0, :pi0]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 200
                }

algorithm
  .note(:total_photon_energy,
        "Sum of energies of the nine photons required to lie in " \
        "3.45 - 3.80 GeV.")
  .note(:pi0_mass_window,
        "pi0 candidate M(gamma gamma) window (110, 150) MeV/c^2.")
  .note(:best_pairing,
        "Best combination of the nine photons into one radiative gamma + " \
        "four pi0 chosen by minimizing chi^2_4pi = Sum_i (m_gg,i - m_pi0)^2 " \
        "/ sigma^2 with sigma = 6.5 MeV/c^2, requiring chi^2_4pi < 15.")
  .note(:jpsi_dipi0_veto,
        "Veto |m_R(pi0 pi0) - m_J/psi| < 100 MeV/c^2 to suppress " \
        "psi' -> pi0 pi0 J/psi transitions (recoil mass of any di-pi0).")
  .note(:kskS_veto,
        "Exclude chi_c0/chi_c2 -> K_S K_S region: reject events with " \
        "sqrt((m12 - m_KS)^2 + (m34 - m_KS)^2) < 100 MeV/c^2 for any " \
        "assignment of the four pi0 to two di-pi0 pairs.")

algorithm.with_decay_card(decay_card_chic0).apply(event_selection)
algorithm.execute_on([psip_data, psip_incMC, exMC_chic0, exMC_chic1, exMC_chic2])
