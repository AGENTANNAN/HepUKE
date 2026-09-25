# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
# J/psi peak data and inclusive MC (sample name = BOSS_version_CMS_energy)
jpsi_data  = DatasetManager.real_data.find("708_3097")      # Real J/psi data at 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # Corresponding inclusive J/psi MC

# --- Decay card: charged topology J/psi -> gamma eta, eta -> pi+ pi- ---
decay_card_charged = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta    PHSP;
    Enddecay

    Decay eta
    1.0000 pi+ pi-    PHSP;
    Enddecay

    End
DECAYCARD

# --- Decay card: neutral topology J/psi -> gamma eta, eta -> pi0 pi0 -> 4 gamma ---
decay_card_neutral = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma eta    PHSP;
    Enddecay

    Decay eta
    1.0000 pi0 pi0    PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma    PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC: 200k events for charged topology (gamma pi+ pi-)
exMC_charged = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_gamma_eta_charged"
  config.related_dataset = jpsi_data
  config.events         = 200000
  config.decay_card     = decay_card_charged
  config.cross_section  = :default
end

# Exclusive MC: 200k events for neutral topology (gamma pi0 pi0 -> 4 gamma)
exMC_neutral = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_gamma_eta_neutral"
  config.related_dataset = jpsi_data
  config.events         = 200000
  config.decay_card     = decay_card_neutral
  config.cross_section  = :default
end

### Event selection (BOSS) ###
# Two independent final states (gamma pi+ pi- and gamma pi0 pi0) => two Algorithm objects.

# ===== Mode I: charged topology  J/psi -> gamma pi+ pi- =====
alg_name_ch = "JpsiGammaEtaCharged"
alg_ch = Algorithm.new(alg_name_ch)
alg_ch.set_header(["#{alg_name_ch}Alg/#{alg_name_ch}.h"])
      .set_constant({"ECMS" => [:double, 3.097]})   # J/psi peak energy (GeV)

sel_ch = Selection.new
  .select_track {                    # exactly one pi+ and one pi-
      cos_theta 0.93                 # |cos theta| < 0.93
      Vz        10.0                 # |Vz| < 10 cm
      Vr        1.0                  # |Vr| < 1 cm
      nChrp     "==1"
      nChrn     "==1"
      nNet      "==0"
  }
  .select_photon {                   # at least one radiative photon
      energyThreshold_b 0.025        # 25 MeV (barrel)
      energyThreshold_e 0.050        # 50 MeV (endcap)
      tdc_emc_start     0            # EMC timing window 0-14
      tdc_emc_end       14
      angle_to_track    20.0         # >= 20 degrees from any charged track
      nGam              ">=1"
  }
  .pid(method: :probability) {       # probability PID, cut 0.001
      prob_cut 0.001
      identify :pion, against: [:kaon]   # Prob(pi) > Prob(K); pi+ and pi- together
      npip "==1"
      npim "==1"
  }
  # Nominal 4C fit to gamma pi+ pi-
  .kinematic_fit([:gamma, :pip, :pim]) {
      nominal
      constrain_four_momentum
      chi2_cut 200                   # loose cut; tight cut optimized in ROOT
  }
  # Competing-hypothesis fit: reassign pions as kaons, store chi2 for gamma K+ K-
  .assign({:pip => :kp, :pim => :km})
  .kinematic_fit([:gamma, :kp, :km]) {
      constrain_four_momentum         # no chi2_cut / no nominal -> only chi2 stored for ROOT veto
  }

alg_ch.with_decay_card(decay_card_charged).apply(sel_ch)

# ===== Mode II: neutral topology  J/psi -> gamma pi0 pi0 -> 5 gamma =====
alg_name_neu = "JpsiGammaEtaNeutral"
alg_neu = Algorithm.new(alg_name_neu)
alg_neu.set_header(["#{alg_name_neu}Alg/#{alg_name_neu}.h"])
       .set_constant({"ECMS" => [:double, 3.097]})

sel_neu = Selection.new
  .select_track {                    # no charged tracks required
      nChrp "==0"
      nChrn "==0"
      nNet  "==0"
  }
  .select_photon {                   # at least five photons (1 radiative + 4 from 2 pi0)
      energyThreshold_b 0.025        # same energy thresholds
      energyThreshold_e 0.050
      angle_to_track    20.0         # same angular separation from tracks
      nGam              ">=5"        # no EMC timing cut for this mode
  }
  # Reconstruct two pi0 from photon pairs, each constrained to the nominal pi0 mass
  .kalman_kinematic_fit([:gamma, :gamma]) {
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 25
      npi0 ">=2"
  }
  # Nominal 4C fit to gamma pi0 pi0
  .kinematic_fit([:gamma, :pi0, :pi0]) {
      nominal
      constrain_four_momentum
      chi2_cut 200
  }

# The most energetic photon is used as the radiative photon in the fit combination.
alg_neu.note(:radiative_photon_selection,
  "the most energetic photon among the selected showers is taken as the radiative
   photon in the gamma pi0 pi0 final state; the remaining photon pairs are used to
   build the two pi0 candidates")

alg_neu.with_decay_card(decay_card_neutral).apply(sel_neu)

# Execute both algorithms on real data, inclusive MC and their exclusive MC samples.
root_files_ch  = alg_ch.execute_on([jpsi_data, jpsi_incMC, exMC_charged])
root_files_neu = alg_neu.execute_on([jpsi_data, jpsi_incMC, exMC_neutral])