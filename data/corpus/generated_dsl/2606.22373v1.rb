# Core DSL classes and dependencies are loaded automatically at execution

### Dataset description ###
# J/psi data at sqrt(s) = 3.097 GeV, with its inclusive MC sample
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Decay card for J/psi -> gamma eta', eta' -> pi+ pi- eta,
# with BOTH charge-conjugate eta decays eta -> e+ mu- and eta -> e- mu+
decay_card_signal = <<~DECAYCARD
    Alias eta_epmum eta
    Alias eta_emmup eta

    Decay J/psi
    1.0000 gamma eta'                PHSP;
    Enddecay

    Decay eta'
    0.5000 pi+ pi- eta_epmum         PHSP;
    0.5000 pi+ pi- eta_emmup         PHSP;
    Enddecay

    Decay eta_epmum
    1.0000 e+ mu-                    PHSP;
    Enddecay

    Decay eta_emmup
    1.0000 e- mu+                    PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for the full signal chain (500k events)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_gamma_etaprime_emu"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "JpsiGammaEtaPrimeEMu"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.097]})           # sqrt(s) = 3.097 GeV
   .set_alias({"std::vector<double>" => "Vdouble"})
   .note(:additional_lepton_id,
         "additional lepton ID applied to the signal l+ l- : the electron candidate must " \
         "satisfy E/p > 0.8; the muon candidate must have an EMC energy deposit in " \
         "(0.1, 0.3) GeV")

event_selection = Selection.new
event_selection
  .select_track {
    cos_theta 0.93      # |cos(theta)| < 0.93
    Vz        10.0      # |Vz| < 10 cm
    Vr        1.0       # Vr < 1 cm
    nChrp     "==2"     # exactly two positive tracks
    nChrn     "==2"     # exactly two negative tracks
    nNet      "==0"     # net charge zero
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0     # > 10 deg from nearest charged track
    energyThreshold_b 0.025    # barrel E > 25 MeV
    energyThreshold_e 0.050    # endcap E > 50 MeV
    nGam              ">=1"    # at least one photon
  }
  .pid(method: :probability) {
    prob_cut 0.001
    # high-momentum tracks (p > 1.0 GeV) -> leptons; e if EMC energy > 0.6 GeV, else mu
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    identify :pion, against: [:kaon, :proton]
    npip "==1"                 # exactly 1 pi+
    npim "==1"                 # exactly 1 pi-
    nlp  "==1"                 # exactly 1 l+
    nlm  "==1"                 # exactly 1 l-
  }
  # Nominal 4C kinematic fit to gamma pi+ pi- l+ l-
  .kinematic_fit([:gamma, :pip, :pim, :lp, :lm]) {
    nominal
    constrain_four_momentum
    chi2_cut 36
    invariant_mass_of(:pip, :pim, :lp, :lm).within(0.946, 0.970)  # eta' signal window
    invariant_mass_of(:lp, :lm).within(0.538, 0.558)              # eta signal window
  }
  # Competing gamma pi+ pi- e+ e- hypothesis (chi2 stored; veto applied in ROOT)
  .assign({:lp => :ep, :lm => :em})
  .kinematic_fit([:gamma, :pip, :pim, :ep, :em]) {
    constrain_four_momentum
  }
  # Competing gamma pi+ pi- pi+ pi- hypothesis (chi2 stored; veto applied in ROOT)
  .assign({:lp => :pip, :lm => :pim})
  .kinematic_fit([:gamma, :pip, :pim, :pip, :pim]) {
    constrain_four_momentum
  }

# Generate the algorithm for the signal decay card and apply the selection chain
alg.with_decay_card(decay_card_signal).apply(event_selection)
# Execute on real data, inclusive MC and exclusive signal MC
root_files = alg.execute_on([jpsi_data, jpsi_incMC, exMC_signal])