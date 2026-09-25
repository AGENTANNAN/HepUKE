### Dataset description ###
# J/psi peak: sqrt(s) = 3.097 GeV (1.31e9 events) -> BOSS sample 708_3097
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

### Decay cards (EvtGen syntax) ###
# Invisible mode: J/psi -> omega eta, omega invisible, eta -> pi+ pi- pi0
decay_card_invisible = <<~DECAYCARD
    Decay J/psi
    1.000 omega eta PHSP;
    Enddecay

    Decay omega
    1.000 nu_mu anti-nu_mu PHSP;
    Enddecay

    Decay eta
    1.000 pi+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Visible reference: J/psi -> omega eta, omega -> pi+ pi- pi0, eta -> pi+ pi- pi0
decay_card_omega = <<~DECAYCARD
    Decay J/psi
    1.000 omega eta PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay eta
    1.000 pi+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Visible reference: J/psi -> phi eta, phi -> K+ K-, eta -> pi+ pi- pi0
decay_card_phi = <<~DECAYCARD
    Decay J/psi
    1.000 phi eta PHSP;
    Enddecay

    Decay phi
    1.000 K+ K- VSS;
    Enddecay

    Decay eta
    1.000 pi+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive MC samples (100k events each) ###
exMC_invisible = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_3097_omega_invisible_eta"
  config.related_dataset = jpsi_data
  config.events = 100000
  config.decay_card = decay_card_invisible
  config.cross_section = :default
end

exMC_omega = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_3097_omega_visible_eta"
  config.related_dataset = jpsi_data
  config.events = 100000
  config.decay_card = decay_card_omega
  config.cross_section = :default
end

exMC_phi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_3097_phi_visible_eta"
  config.related_dataset = jpsi_data
  config.events = 100000
  config.decay_card = decay_card_phi
  config.cross_section = :default
end

### Event selection (BOSS) ###
# All three modes use the same track quality cuts, photon cuts and PID method.

# ------------------------------ Mode 1: invisible (J/psi -> omega(inv) eta, eta -> pi+pi-pi0) --------------
alg_inv = Algorithm.new("OmegainvEta")
alg_inv.set_header(["OmegainvEtaAlg/OmegainvEta.h"])
       .set_constant({ "ECMS" => [:double, 3.097] })
       .set_alias({ "std::vector<double>" => "Vdouble" })

sel_inv = Selection.new
  .select_track {                 # charged track quality + multiplicity
    cos_theta 0.93                # |cos(theta)| < 0.93
    Vz        100.0               # |Vz| < 100 (beam direction)
    Vr        10.0                # Vr < 10 (transverse plane)
    nChrp     "==1"               # exactly one positive track
    nChrn     "==1"               # exactly one negative track
    nNet      "==0"               # net charge zero
  }
  .select_photon {                # photon quality
    tdc_emc_start 0               # EMC timing window start
    tdc_emc_end   14              # EMC timing window end
    angle_to_track 10.0           # angle to nearest charged track > 10 deg
    energyThreshold_b 0.025       # E > 25 MeV in the barrel
    energyThreshold_e 0.050       # E > 50 MeV in the endcap
    nGam ">=2"                    # at least two photons
  }
  .pid(method: :probability) {    # PID: pions against kaons
    prob_cut 0.001
    identify :pion, against: [:kaon]
    npip "==1"
    npim "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # build pi0 from a gamma-gamma pair
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:pip, :pim, :pi0]) {        # 4C fit to pi+pi-pi0 (eta)
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_inv.with_decay_card(decay_card_invisible).apply(sel_inv)
alg_inv.execute_on([jpsi_data, jpsi_incMC, exMC_invisible])

# ------------------------------ Mode 2: visible omega reference (J/psi -> omega eta, omega -> pi+pi-pi0, eta -> pi+pi-pi0) --------------
alg_omega = Algorithm.new("OmegaRefEta")
alg_omega.set_header(["OmegaRefEtaAlg/OmegaRefEta.h"])
         .set_constant({ "ECMS" => [:double, 3.097] })
         .set_alias({ "std::vector<double>" => "Vdouble" })

sel_omega = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        100.0
    Vr        10.0
    nChrp     "==2"               # two pi+ (from omega and eta)
    nChrn     "==2"               # two pi-
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end   14
    angle_to_track 10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=4"                    # at least four photons -> two pi0
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon]
    npip "==2"
    npim "==2"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # build pi0 candidates
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=2"
  }
  .kinematic_fit([:pip, :pip, :pim, :pim, :pi0, :pi0]) {   # 4C fit to pi+pi-pi0 pi+pi-pi0
    nominal
    constrain_four_momentum
    invariant_mass_of(:pip, :pim, :pi0).within(0.65, 0.98)  # M(omega) window
    invariant_mass_of(:pip, :pim, :pi0).within(0.41, 0.65)  # M(eta) window
    chi2_cut 200
  }

alg_omega.with_decay_card(decay_card_omega).apply(sel_omega)
alg_omega.execute_on([jpsi_data, jpsi_incMC, exMC_omega])

# ------------------------------ Mode 3: visible phi reference (J/psi -> phi eta, phi -> K+K-, eta -> pi+pi-pi0) --------------
alg_phi = Algorithm.new("PhiRefEta")
alg_phi.set_header(["PhiRefEtaAlg/PhiRefEta.h"])
       .set_constant({ "ECMS" => [:double, 3.097] })
       .set_alias({ "std::vector<double>" => "Vdouble" })

sel_phi = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        100.0
    Vr        10.0
    nChrp     "==2"               # two positive tracks (K+ and pi+)
    nChrn     "==2"               # two negative tracks (K- and pi-)
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end   14
    angle_to_track 10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=2"                    # at least two photons -> one pi0
  }
  .pid(method: :probability) {    # shared PID: pions against kaons
    prob_cut 0.001
    identify :pion, against: [:kaon]
  }
  .assign({ :chrgp => :kp, :chrgn => :km })   # phi -> K+K- without PID: charged tracks as (anti-)kaon candidates
  .kalman_kinematic_fit([:gamma, :gamma]) {   # build pi0 from a gamma-gamma pair
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  .kinematic_fit([:kp, :km, :pip, :pim, :pi0]) {   # 4C fit to K+K- pi+pi- pi0
    nominal
    constrain_four_momentum
    invariant_mass_of(:kp, :km).within(0.987, 1.10)  # M(K+K-) = M(phi) window
    chi2_cut 200
  }

alg_phi.with_decay_card(decay_card_phi).apply(sel_phi)
alg_phi.execute_on([jpsi_data, jpsi_incMC, exMC_phi])