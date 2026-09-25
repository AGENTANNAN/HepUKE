# Dataset preparation
jpsi_data   = DatasetManager.real_data.find("708_3097")     # J/psi data (2009+2012, 1.311e9 events)
jpsi_incMC  = DatasetManager.inclusive_mc.find("708_3097")  # 1.2e9 inclusive J/psi MC

# ---------------------------------------------------------------------------
# Decay cards
# ---------------------------------------------------------------------------

# Signal: J/psi -> phi pi0 f0(980), f0(980) -> pi+ pi-, phi -> K+ K-
decay_card_phi_pi0_f0_charged = <<~DECAYCARD
    Decay J/psi
    1.000 phi pi0 f_0                                PHSP;
    Enddecay

    Decay phi
    1.000 K+ K-                                      VSS;
    Enddecay

    Decay f_0
    1.000 pi+ pi-                                    PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma                                PHSP;
    Enddecay

    End
DECAYCARD

# Signal: J/psi -> phi pi0 f0(980), f0(980) -> pi0 pi0, phi -> K+ K-
decay_card_phi_pi0_f0_neutral = <<~DECAYCARD
    Decay J/psi
    1.000 phi pi0 f_0                                PHSP;
    Enddecay

    Decay phi
    1.000 K+ K-                                      VSS;
    Enddecay

    Decay f_0
    1.000 pi0 pi0                                    PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma                                PHSP;
    Enddecay

    End
DECAYCARD

# Signal: J/psi -> phi f1(1285), f1 -> pi0 f0(980), f0 -> pi+ pi-
decay_card_phi_f1_charged = <<~DECAYCARD
    Decay J/psi
    1.000 phi f_1                                    SVV_HELAMP 1.0 0.0 1.0 0.0 1.0 0.0;
    Enddecay

    Decay phi
    1.000 K+ K-                                      VSS;
    Enddecay

    Decay f_1
    1.000 pi0 f_0                                    PHSP;
    Enddecay

    Decay f_0
    1.000 pi+ pi-                                    PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma                                PHSP;
    Enddecay

    End
DECAYCARD

# Signal: J/psi -> phi f1(1285), f1 -> pi0 f0(980), f0 -> pi0 pi0
decay_card_phi_f1_neutral = <<~DECAYCARD
    Decay J/psi
    1.000 phi f_1                                    SVV_HELAMP 1.0 0.0 1.0 0.0 1.0 0.0;
    Enddecay

    Decay phi
    1.000 K+ K-                                      VSS;
    Enddecay

    Decay f_1
    1.000 pi0 f_0                                    PHSP;
    Enddecay

    Decay f_0
    1.000 pi0 pi0                                    PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma                                PHSP;
    Enddecay

    End
DECAYCARD

# Signal: J/psi -> phi eta', eta' -> pi+ pi- pi0
decay_card_phi_etap_charged = <<~DECAYCARD
    Decay J/psi
    1.000 phi eta'                                   SVS;
    Enddecay

    Decay phi
    1.000 K+ K-                                      VSS;
    Enddecay

    Decay eta'
    1.000 pi+ pi- pi0                                PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma                                PHSP;
    Enddecay

    End
DECAYCARD

# Signal: J/psi -> phi eta', eta' -> pi0 pi0 pi0
decay_card_phi_etap_neutral = <<~DECAYCARD
    Decay J/psi
    1.000 phi eta'                                   SVS;
    Enddecay

    Decay phi
    1.000 K+ K-                                      VSS;
    Enddecay

    Decay eta'
    1.000 pi0 pi0 pi0                                PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma                                PHSP;
    Enddecay

    End
DECAYCARD

# ---------------------------------------------------------------------------
# Exclusive MC samples
# ---------------------------------------------------------------------------
exMC_phi_pi0_f0_ch = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "phi_pi0_f0_pipi"
  c.related_dataset = jpsi_data
  c.events          = 200000
  c.decay_card      = decay_card_phi_pi0_f0_charged
  c.cross_section   = :default
end

exMC_phi_pi0_f0_ne = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "phi_pi0_f0_pi0pi0"
  c.related_dataset = jpsi_data
  c.events          = 200000
  c.decay_card      = decay_card_phi_pi0_f0_neutral
  c.cross_section   = :default
end

exMC_phi_f1_ch = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "phi_f1_pi0_f0_pipi"
  c.related_dataset = jpsi_data
  c.events          = 200000
  c.decay_card      = decay_card_phi_f1_charged
  c.cross_section   = :default
end

exMC_phi_f1_ne = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "phi_f1_pi0_f0_pi0pi0"
  c.related_dataset = jpsi_data
  c.events          = 200000
  c.decay_card      = decay_card_phi_f1_neutral
  c.cross_section   = :default
end

exMC_phi_etap_ch = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "phi_etap_pippimpi0"
  c.related_dataset = jpsi_data
  c.events          = 200000
  c.decay_card      = decay_card_phi_etap_charged
  c.cross_section   = :default
end

exMC_phi_etap_ne = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "phi_etap_3pi0"
  c.related_dataset = jpsi_data
  c.events          = 200000
  c.decay_card      = decay_card_phi_etap_neutral
  c.cross_section   = :default
end

# ===========================================================================
# Channel 1: J/psi -> K+ K- pi+ pi- pi0  (5C kinematic fit)
# ===========================================================================
alg_charged = Algorithm.new("PhiPi0F0Charged")
alg_charged.set_header(["PhiPi0F0ChargedAlg/PhiPi0F0Charged.h"])
           .set_constant({ "ECMS" => [:double, 3.097] })

sel_charged = Selection.new
sel_charged.select_track {
              cos_theta 0.93
              Vz        10.0
              Vr        1.0
              nChrp     "==2"
              nChrn     "==2"
              nNet      "==0"
            }
           .select_photon {
              tdc_emc_start     0
              tdc_emc_end       14
              angle_to_track    10.0
              energyThreshold_b 0.025
              energyThreshold_e 0.050
              nGam              ">=2"
            }
           .pid(method: :probability) {
              prob_cut 0.0
              identify :kaon, against: [:pion]
              nkp ">=1"
              nkm ">=1"
            }
           .remove([:kp <= :chrgp, :km <= :chrgn])
           .assign({ :chrgp => :pip, :chrgn => :pim })
           .kinematic_fit([:kp, :km, :pip, :pim, :gamma, :gamma]) {
              nominal
              constrain_four_momentum
              invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
              chi2_cut 100
            }

alg_charged.note(:kstar_veto,
                 "events with |M(K+-pi-+) - M(K*0)| < 0.050 GeV/c^2 are rejected " \
                 "to suppress K*0 / K*0bar intermediate-state backgrounds; " \
                 "applied at ROOT level using kinematic-fit-corrected momenta")
           .note(:kaon_assignment,
                 "the two oppositely charged tracks with invariant mass closest to " \
                 "the nominal phi mass are assigned as kaons (chi2 minimisation of the " \
                 "5C fit performs the effective assignment)")
           .with_decay_card(decay_card_phi_pi0_f0_charged)
           .apply(sel_charged)

alg_charged.execute_on([jpsi_data, jpsi_incMC,
                        exMC_phi_pi0_f0_ch, exMC_phi_f1_ch, exMC_phi_etap_ch])

# ===========================================================================
# Channel 2: J/psi -> K+ K- pi0 pi0 pi0  (7C kinematic fit)
# ===========================================================================
alg_neutral = Algorithm.new("PhiPi0F0Neutral")
alg_neutral.set_header(["PhiPi0F0NeutralAlg/PhiPi0F0Neutral.h"])
           .set_constant({ "ECMS" => [:double, 3.097] })

sel_neutral = Selection.new
sel_neutral.select_track {
              cos_theta 0.93
              Vz        10.0
              Vr        1.0
              nChrp     "==1"
              nChrn     "==1"
              nNet      "==0"
            }
           .select_photon {
              tdc_emc_start     0
              tdc_emc_end       14
              angle_to_track    10.0
              energyThreshold_b 0.025
              energyThreshold_e 0.050
              nGam              ">=6"
            }
           .pid(method: :probability) {
              prob_cut 0.0
              identify :kaon, against: [:pion]
              nkp ">=1"
              nkm ">=1"
            }
           .kinematic_fit([:kp, :km, :gamma, :gamma, :gamma, :gamma, :gamma, :gamma]) {
              nominal
              constrain_four_momentum
              invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
              invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
              invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
              chi2_cut 90
            }

alg_neutral.note(:photon_pairing,
                 "the six selected photons are paired into three pi0 candidates by " \
                 "minimising sum_i ((M(gamma_i gamma_j) - M_pi0)^2 / sigma_pi0^2); " \
                 "the 7C kinematic fit's chi2 minimisation over photon combinations " \
                 "performs the equivalent selection")
           .with_decay_card(decay_card_phi_pi0_f0_neutral)
           .apply(sel_neutral)

alg_neutral.execute_on([jpsi_data, jpsi_incMC,
                        exMC_phi_pi0_f0_ne, exMC_phi_f1_ne, exMC_phi_etap_ne])
