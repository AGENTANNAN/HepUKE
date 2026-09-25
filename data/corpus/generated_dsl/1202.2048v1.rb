# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi (3.097 GeV) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # J/psi inclusive MC
psip_data  = DatasetManager.real_data.find("709_3686")      # psi(2S) (3.686 GeV) real data
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # psi(2S) inclusive MC

# Decay card: J/psi -> rho0 pi0, rho0 -> pi+ pi-, pi0 -> gamma gamma
decay_card_jpsi = <<~DECAYCARD
  Decay J/psi
  1.000 rho0 pi0 PHSP;
  Enddecay

  Decay rho0
  1.000 pi+ pi- VSS;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Decay card: psi(2S) -> pi+ pi- pi0 (phase space), pi0 -> gamma gamma
decay_card_psip = <<~DECAYCARD
  Decay psi(2S)
  1.000 pi+ pi- pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Exclusive MC: 1M events for each signal mode
exMC_jpsi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_jpsi_rho0pi0"
  config.related_dataset = jpsi_data
  config.events = 1_000_000
  config.decay_card = decay_card_jpsi
  config.cross_section = :default
end

exMC_psip = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_psip_pipipi0"
  config.related_dataset = psip_data
  config.events = 1_000_000
  config.decay_card = decay_card_psip
  config.cross_section = :default
end

### Event selection (BOSS) ###

# ---------- J/psi -> pi+ pi- pi0 ----------
alg_jpsi_name = "JpsiPiPiPi0"
alg_jpsi = Algorithm.new(alg_jpsi_name)
alg_jpsi.set_header(["#{alg_jpsi_name}Alg/#{alg_jpsi_name}.h"])
        .set_constant({"ECMS" => [:double, 3.097]})

sel_jpsi = Selection.new
  .select_track {
    cos_theta 0.93        # |cos(theta)| < 0.93
    Vz 10.0              # |Vz| < 10 cm
    Vr 1.0               # Vr < 1 cm
    nChrp "==1"          # exactly one positive track
    nChrn "==1"          # exactly one negative track
    nNet "==0"           # net charge zero
  }
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=2"           # at least two photons
  }
  # no hadron PID: assign the charged tracks directly as pi+ pi-
  .assign({:chrgp => :pip, :chrgn => :pim})
  # Kalman fit: constrain the photon pair to the nominal pi0 mass (best pair kept)
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 50
    npi0 ">=1"           # at least one pi0
  }
  # Nominal 5C fit: 4C four-momentum constraint + pi0 mass constraint
  .kinematic_fit([:pip, :pim, :pi0]) {
    nominal
    constrain_four_momentum
    chi2_cut 50
  }
  # Competing K+ K- pi0 hypothesis: store chi2 only (no chi2_cut, not nominal);
  # the veto against the kaon hypothesis is applied in ROOT
  .assign({:pip => :kp, :pim => :km})
  .kinematic_fit([:kp, :km, :pi0]) {
    constrain_four_momentum
  }

alg_jpsi.with_decay_card(decay_card_jpsi).apply(sel_jpsi)
root_files_jpsi = alg_jpsi.execute_on([jpsi_data, jpsi_incMC, exMC_jpsi])

# ---------- psi(2S) -> pi+ pi- pi0 ----------
alg_psip_name = "PsipPiPiPi0"
alg_psip = Algorithm.new(alg_psip_name)
alg_psip.set_header(["#{alg_psip_name}Alg/#{alg_psip_name}.h"])
        .set_constant({"ECMS" => [:double, 3.686]})
        # psi(2S)-only extra vetoes (not expressible with dedicated DSL methods)
        .note(:electron_suppression, "for the psi(2S) sample, charged tracks with associated EMC energy above 0.8 GeV are removed to suppress electrons")
        .note(:muon_veto, "for the psi(2S) sample, an analysis-level veto requires MUC penetration depth < 40 cm")
        .note(:pipim_mass_veto, "for the psi(2S) sample, an analysis-level veto requires m(pi+pi-) < 3 GeV/c^2")

sel_psip = Selection.new
  .select_track {
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp "==1"
    nChrn "==1"
    nNet "==0"
  }
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=2"
  }
  # no hadron PID: assign the charged tracks directly as pi+ pi-
  .assign({:chrgp => :pip, :chrgn => :pim})
  # Kalman fit: constrain the photon pair to the nominal pi0 mass (best pair kept)
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 50
    npi0 ">=1"
  }
  # Nominal 5C fit: 4C four-momentum constraint + pi0 mass constraint
  .kinematic_fit([:pip, :pim, :pi0]) {
    nominal
    constrain_four_momentum
    chi2_cut 50
  }
  # Competing K+ K- pi0 hypothesis: store chi2 only (no chi2_cut, not nominal)
  .assign({:pip => :kp, :pim => :km})
  .kinematic_fit([:kp, :km, :pi0]) {
    constrain_four_momentum
  }

alg_psip.with_decay_card(decay_card_psip).apply(sel_psip)
root_files_psip = alg_psip.execute_on([psip_data, psip_incMC, exMC_psip])