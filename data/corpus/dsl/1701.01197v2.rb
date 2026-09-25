### Dataset preparation ###
psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Decay card: psi(3686) -> gamma chi_cJ -> gamma gamma J/psi -> gamma gamma l+ l-
decay_card_chic = <<~DECAYCARD
    Decay psi(2S)
    0.3333 gamma chi_c0                       HELAMP 1.0 0.0 1.0 0.0;
    0.3333 gamma chi_c1                       HELAMP 1.0 0.0 1.0 0.0 -1.0 0.0 -1.0 0.0;
    0.3334 gamma chi_c2                       HELAMP 1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.0;
    Enddecay

    Decay chi_c0
    1.0000 gamma J/psi                        PHSP;
    Enddecay

    Decay chi_c1
    1.0000 gamma J/psi                        HELAMP 1.0 0.0 1.0 0.0 -1.0 0.0 -1.0 0.0;
    Enddecay

    Decay chi_c2
    1.0000 gamma J/psi                        HELAMP 1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.0;
    Enddecay

    Decay J/psi
    0.5000 e+ e-                              PHOTOS VLL;
    0.5000 mu+ mu-                            PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# Decay card: psi(3686) -> gamma eta_c(2S) -> gamma gamma J/psi
decay_card_etac2s = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma eta_c(2S)                    HELAMP 1.0 0.0 -1.0 0.0;
    Enddecay

    Decay eta_c(2S)
    1.0000 gamma J/psi                        HELAMP 1.0 0.0 -1.0 0.0;
    Enddecay

    Decay J/psi
    0.5000 e+ e-                              PHOTOS VLL;
    0.5000 mu+ mu-                            PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# Background reference: psi(3686) -> pi0 J/psi (JPIPI)
decay_card_pi0Jpsi = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi0 J/psi                          PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                        PHSP;
    Enddecay

    Decay J/psi
    0.5000 e+ e-                              PHOTOS VLL;
    0.5000 mu+ mu-                            PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# Background reference: psi(3686) -> pi0 pi0 J/psi
decay_card_pi0pi0Jpsi = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi0 pi0 J/psi                      PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma                        PHSP;
    Enddecay

    Decay J/psi
    0.5000 e+ e-                              PHOTOS VLL;
    0.5000 mu+ mu-                            PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# Background reference: psi(3686) -> gamma gamma J/psi (non-resonant)
decay_card_ggJpsi = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma gamma J/psi                  PHSP;
    Enddecay

    Decay J/psi
    0.5000 e+ e-                              PHOTOS VLL;
    0.5000 mu+ mu-                            PHOTOS VLL;
    Enddecay

    End
DECAYCARD

exMC_chic = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_gamma_chicJ_gamma_gamma_Jpsi_ll_exclusive_mc"
  config.related_dataset = psip_data
  config.events          = 1000000
  config.decay_card      = decay_card_chic
  config.cross_section   = :default
end

exMC_etac2s = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_gamma_etac2s_gamma_gamma_Jpsi_ll_exclusive_mc"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_etac2s
  config.cross_section   = :default
end

exMC_pi0Jpsi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_pi0_Jpsi_exclusive_mc"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_pi0Jpsi
  config.cross_section   = :default
end

exMC_pi0pi0Jpsi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_pi0pi0_Jpsi_exclusive_mc"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_pi0pi0Jpsi
  config.cross_section   = :default
end

exMC_ggJpsi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_gamma_gamma_Jpsi_nonresonant_exclusive_mc"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_ggJpsi
  config.cross_section   = :default
end

########################################################################
# Signal channel: psi(3686) -> gamma1 chi_cJ / eta_c(2S) -> gamma1 gamma2 J/psi
# J/psi -> e+ e- or mu+ mu-. Exactly two oppositely charged tracks and 2-4 photons.
########################################################################
alg = Algorithm.new("PsipGGll")
alg.set_header(["PsipGGllAlg/PsipGGll.h"])
   .set_constant({"ECMS" => [:double, 3.686]})

sel = Selection.new
sel.select_track {
      cos_theta 0.93
      Vz        10.0
      Vr        1.0
      nChrp     "==1"
      nChrn     "==1"
      nNet      "==0"
      p         ">1.0"        # momentum > 1 GeV/c
    }
   .select_photon {
      tdc_emc_start     0
      tdc_emc_end       14
      angle_to_track    10.0
      energyThreshold_b 0.025
      energyThreshold_e 0.025
      nGam              ">=2"
   }
   .pid(method: :probability) {
      prob_cut 0.001
      identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                     treat_as_electron_if_energy_above: 0.6
      nlp ">=1"
      nlm ">=1"
   }
   # 4C kinematic fit constraining to the psi(3686) 4-momentum;
   # choose the two-photon combination with smallest chi2_4C.
   .kinematic_fit([:gamma, :gamma, :lp, :lm]) {
      nominal
      constrain_four_momentum
      chi2_cut 60
      invariant_mass_of(:lp, :lm).within(3.08, 3.12)              # J/psi mass window
      invariant_mass_of(:gamma, :gamma).out_of(0.11, 0.15)        # pi0 veto
      invariant_mass_of(:gamma, :gamma).out_of(0.51, 0.57)        # eta veto (M(gg) > 0.51 rejected)
   }

alg.note(:pi0_veto,
         "reject events with M4C(gamma gamma) in (0.11, 0.15) GeV/c^2 to remove psi(3686) -> pi0 J/psi")
   .note(:eta_veto,
         "reject events with M4C(gamma gamma) > 0.51 GeV/c^2 to remove psi(3686) -> eta J/psi")
   .note(:bhabha_polar_cut,
         "for J/psi -> e+ e-: cos(theta_e+) < 0.3 and cos(theta_e-) > -0.3 to suppress radiative Bhabha")
   .note(:kinematic_fit_3c,
         "downstream, a 3C kinematic fit with the soft photon energy left free is used for the M(gamma2 l+ l-) spectrum")
   .note(:nphoton_upper,
         "N_gamma <= 4 candidates required")
   .with_decay_card(decay_card_chic)
   .apply(sel)

alg.execute_on([psip_data, psip_incMC, exMC_chic, exMC_etac2s,
                exMC_pi0Jpsi, exMC_pi0pi0Jpsi, exMC_ggJpsi])
