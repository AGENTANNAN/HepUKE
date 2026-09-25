# ============================================================================
# e+e- -> eta psi(2S), psi(2S) -> pi+pi- J/psi, J/psi -> l+l- (l = e or mu),
# eta -> gamma gamma
# BESIII energy scan 4.236 - 4.600 GeV (14 real-data points)
# ============================================================================

### Datasets ###
# The 14 real-data energy points of the 4.236 - 4.600 GeV scan (BOSS 703)
data_points = [
  DatasetManager.real_data.find("703_4237"),   # 4235.7 MeV
  DatasetManager.real_data.find("703_4245"),   # 4241.7 MeV
  DatasetManager.real_data.find("703_4246"),   # 4243.8 MeV
  DatasetManager.real_data.find("703_4260"),   # 4258.0 MeV
  DatasetManager.real_data.find("703_4270"),   # 4266.8 MeV
  DatasetManager.real_data.find("703_4280"),   # 4277.7 MeV
  DatasetManager.real_data.find("703_4310"),   # 4307.9 MeV
  DatasetManager.real_data.find("703_4360"),   # 4358.3 MeV
  DatasetManager.real_data.find("703_4390"),   # 4387.4 MeV
  DatasetManager.real_data.find("703_4420"),   # 4415.6 MeV
  DatasetManager.real_data.find("703_4470"),   # 4467.1 MeV
  DatasetManager.real_data.find("703_4530"),   # 4527.1 MeV
  DatasetManager.real_data.find("703_4575"),   # 4574.5 MeV
  DatasetManager.real_data.find("703_4600")    # 4599.5 MeV
]

# Inclusive MC sample at 4.258 GeV (503 pb^-1 scan point of the same campaign)
incMC = DatasetManager.inclusive_mc.find("703_4260")

### Decay card for the signal process (EvtGen format) ###
# e+e- -> eta psi(2S); psi(2S) -> pi+pi- J/psi; J/psi -> e+e- (PHOTOS, VLL); eta -> gamma gamma
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.000 eta psi(2S)        PHSP;
    Enddecay

    Decay psi(2S)
    1.000 pi+ pi- J/psi      PHSP;
    Enddecay

    Decay J/psi
    1.000 e+ e-              PHOTOS VLL;
    Enddecay

    Decay eta
    1.000 gamma gamma        PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive MC ###
# Same signal process generated at every scan point: 100k events per energy point
exMCs_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_etapsip"        # becomes sig_etapsip_703_4237, ... per point
  config.events        = 100_000              # 100k events at each energy point
  config.decay_card    = decay_card_signal    # shared decay card for all points
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "EtaPsi2S"
eta_psip_alg = Algorithm.new(alg_name)
eta_psip_alg
  .set_header(["#{alg_name}Alg/#{alg_name}.h"])
  .set_constant({"ECMS" => [:double, 4.26]})   # nominal CMS energy (GeV); per-point sqrt(s) is taken from each dataset at execution
  .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection
  .select_track {                          # Charged track selection
      cos_theta   0.93                     # |cos(theta)| < 0.93
      Vz          10.0                     # |Vz| < 10 cm
      Vr          1.0                      # Vr < 1 cm in the transverse plane
      nChrp       "==2"                    # exactly two positively charged tracks
      nChrn       "==2"                    # exactly two negatively charged tracks
      nNet        "==0"                    # net charge zero (four charged tracks in total)
  }
  .select_photon {                         # Photon selection
      tdc_emc_start     0                  # TDC start time
      tdc_emc_end       14                 # TDC end time
      angle_to_track    10.0               # > 10 degrees from the nearest charged track
      energyThreshold_b 0.025              # 25 MeV minimum shower energy in the barrel
      energyThreshold_e 0.050              # 50 MeV minimum shower energy in the endcap
      nGam              ">=2"              # at least two photons (eta -> gamma gamma)
  }
  .pid(method: :probability) {             # Probability-based particle identification
      # High-momentum lepton identification: p > 1.0 GeV/c -> lepton;
      # EMC energy > 1.0 GeV -> electron, otherwise muon
      identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                     treat_as_electron_if_energy_above: 1.0
      identify :pion, against: [:kaon]     # pi+/pi- separation against kaons
      npip "==1"                           # one pi+
      npim "==1"                           # one pi-
      nlp  "==1"                           # one l+
      nlm  "==1"                           # one l-
  }
  # 4C kinematic fit to gamma gamma pi+ pi- l+ l- (best chi2 combination retained)
  .kinematic_fit([:gamma, :gamma, :pip, :pim, :lp, :lm]) {
      nominal                              # nominal fit: its corrected four-momenta are used
      constrain_four_momentum              # 4C energy-momentum constraint
      chi2_cut 40                          # chi2 < 40
  }
  # NOTE: the J/psi (3064.6-3140.8 MeV), eta (507.1-579.1 MeV) and psi(2S)
  # (3680.3-3692.5 MeV) mass windows, the eta' -> pi+pi- eta veto
  # (M(pi+pi-gamma gamma) > 1.0 GeV) and the smallest-chi2_4C multi-photon
  # combination choice are applied downstream in the ROOT analysis.

# Render the selection into the BOSS algorithm and run it
eta_psip_alg.with_decay_card(decay_card_signal).apply(event_selection)

root_files = eta_psip_alg.execute_on(data_points + [incMC] + exMCs_signal)