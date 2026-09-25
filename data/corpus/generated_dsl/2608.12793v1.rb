# frozen_string_literal: true
# BOSS-side DSL for e+e- -> e+e- eta' via single-tag two-photon fusion at sqrt(s) = 3.773 GeV
#   Mode I : eta' -> pi+ pi- gamma
#   Mode II: eta' -> pi+ pi- eta, eta -> gamma gamma

### Dataset description ###
psi3770_data  = DatasetManager.real_data.find("712_3773")      # 20.3 fb^-1 psi(3770) data, sqrt(s) = 3.773 GeV
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")   # corresponding inclusive MC

# Decay card — Mode I : eta' -> pi+ pi- gamma
decay_card_modeI = <<~DECAYCARD
    Decay psi(3770)
    1.000 e+ e- eta'  PHSP;
    Enddecay

    Decay eta'
    1.000 pi+ pi- gamma  PHSP;
    Enddecay

    End
DECAYCARD

# Decay card — Mode II : eta' -> pi+ pi- eta, eta -> gamma gamma
decay_card_modeII = <<~DECAYCARD
    Decay psi(3770)
    1.000 e+ e- eta'  PHSP;
    Enddecay

    Decay eta'
    1.000 pi+ pi- eta  PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma  PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC — 500k events for each eta' decay mode
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "ee_ee_etap_pipigamma"
  config.related_dataset = psi3770_data
  config.events          = 500_000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "ee_ee_etap_pipieta"
  config.related_dataset = psi3770_data
  config.events          = 500_000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# ==================================================================
# Mode I : eta' -> pi+ pi- gamma
# ==================================================================
alg_name_modeI = "EtapToPipiGammaTag"
alg_modeI = Algorithm.new(alg_name_modeI)
alg_modeI.set_header(["#{alg_name_modeI}Alg/#{alg_name_modeI}.h"])
        .set_constant({"ECMS" => [:double, 3.773]})
        .set_alias({"std::vector<double>" => "Vdouble"})
        .note(:charge_conjugate_branch,
              "charge-conjugate configuration (tagged e- with the opposite-sign untagged e+ missing) is \
treated identically; only the e+ tagged branch is written out explicitly")

sel_modeI = Selection.new
  .select_track {                 # exactly three charged tracks -> net charge +-1
    cos_theta 0.93                # |cos(theta)| < 0.93
    Vz        10.0                # |Vz| < 10 cm
    Vr        1.0                 # Vr < 1 cm
    nChrp     ">=1"               # with >=1 positive and >=1 negative track
    nChrn     ">=1"
    nTot      "==3"               # and exactly three tracks, the net charge is +-1
  }
  .select_photon {                # good photon selection (baseline)
    tdc_emc_start     0
    tdc_emc_end       14          # shower time in [0, 700] ns
    energyThreshold_b 0.025       # > 25 MeV (barrel)
    energyThreshold_e 0.050       # > 50 MeV (endcap)
    angle_to_track    10.0        # > 10 deg from nearest charged track
    nGam              ">=1"       # at least one good photon (eta' -> pi+ pi- gamma)
  }
  .pid(method: :probability) {    # PID: high-momentum leptons + pi/K separation
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.5,
                                   treat_as_electron_if_energy_above: 0.8
    identify :pion, against: [:kaon]
    npip "==1"                    # one pi+
    npim "==1"                    # one pi-
  }
  .remove([:pip <= :chrgp, :pim <= :chrgn])   # the leftover track is the tagged e+-
  .assign({:chrgp => :ep, :chrgn => :em})
  .kinematic_fit([:gamma, :pip, :pim, :ep, :em]) {   # 1C kinematic fit
    nominal
    miss_track_of :em             # untagged lepton treated as missing (electron mass)
    constrain_four_momentum
    chi2_cut 200                  # loose BOSS cut; tight chi2_1C < 70 applied later in ROOT
  }

alg_modeI.with_decay_card(decay_card_modeI).apply(sel_modeI)
alg_modeI.execute_on([psi3770_data, psi3770_incMC, exMC_modeI])

# ==================================================================
# Mode II : eta' -> pi+ pi- eta, eta -> gamma gamma
# ==================================================================
alg_name_modeII = "EtapToPipiEtaTag"
alg_modeII = Algorithm.new(alg_name_modeII)
alg_modeII.set_header(["#{alg_name_modeII}Alg/#{alg_name_modeII}.h"])
         .set_constant({"ECMS" => [:double, 3.773]})
         .set_alias({"std::vector<double>" => "Vdouble"})
         .note(:charge_conjugate_branch,
               "charge-conjugate configuration (tagged e- with the opposite-sign untagged e+ missing) is \
treated identically; only the e+ tagged branch is written out explicitly")

sel_modeII = Selection.new
  .select_track {                 # same three-track baseline
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=1"
    nChrn     ">=1"
    nTot      "==3"
  }
  .select_photon {                # same photon baseline
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=2"       # at least two good photons (eta -> gamma gamma)
  }
  .pid(method: :probability) {    # same PID baseline
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.5,
                                   treat_as_electron_if_energy_above: 0.8
    identify :pion, against: [:kaon]
    npip "==1"
    npim "==1"
  }
  .remove([:pip <= :chrgp, :pim <= :chrgn])
  .assign({:chrgp => :ep, :chrgn => :em})
  .kalman_kinematic_fit([:gamma, :gamma]) {   # reconstruct eta from gamma gamma (1C mass constraint)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=1"                    # require at least one eta candidate
  }
  .kinematic_fit([:gamma, :gamma, :pip, :pim, :ep, :em]) {   # 2C kinematic fit
    nominal
    miss_track_of :em             # untagged lepton treated as missing (electron mass)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)   # M(gg) = m_eta
    constrain_four_momentum
    chi2_cut 200                  # loose BOSS cut; tight chi2_2C < 200 applied later in ROOT
  }

alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)
alg_modeII.execute_on([psi3770_data, psi3770_incMC, exMC_modeII])