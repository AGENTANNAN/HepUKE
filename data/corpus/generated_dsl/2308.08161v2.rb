# =============================================================================
# e+e- -> eta phi,  phi -> K+ K-,  eta -> gamma gamma
# Born cross-section scan at 23 centre-of-mass energies (3.773 - 4.600 GeV)
# BOSS side: dataset preparation + event selection up to the nominal 4C fit.
# The phi mass window, the eta gamma-gamma opening-angle cut and the unbinned
# M(gamma gamma) fit are all post-kinematic-fit (ROOT-level) and are therefore
# NOT expressed here.
# =============================================================================

# ---- 23 c.m. energy points: corresponding real data + matching inclusive MC ----
sample_names = [
  "712_3773", "712_3780",                                        # psi(3770) region
  "703_4009",
  "705_4130", "705_4160",
  "703_4180", "703_4190", "703_4200", "703_4210", "703_4220",
  "703_4230", "703_4237", "703_4246", "703_4260", "703_4270",
  "703_4280", "703_4310", "703_4360", "703_4420",
  "703_4470", "703_4530", "703_4575", "703_4600"
]
data_points  = sample_names.map { |n| DatasetManager.real_data.find(n) }
incMC_points = sample_names.map { |n| DatasetManager.inclusive_mc.find(n) }

# ---- ConExc decay card: mode 23 = vpho -> phi eta (phi -> K+K-, eta -> gamma gamma)
#      No "Particle vpho" line: the framework injects the per-energy-particle
#      for the multi-energy scan (DSL auto-detects the ConExc token).
decay_card_signal = <<~DECAYCARD
  Decay vpho
  1.000 ConExc 23;
  Enddecay

  Decay phi
  1.000 K+ K- VSS;
  Enddecay

  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# ---- ConExc exclusive signal MC: 500k events at each energy point ----
exMC_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_etaphi_conexc"
  config.events        = 500_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

# ---- Algorithm ----
alg_name   = "EtaPhi"
etaphi_alg = Algorithm.new(alg_name)
etaphi_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
          .set_constant({"ECMS" => [:double, 4.26]})   # nominal value; ConExc/scan injects per point
          .set_alias({"std::vector<double>" => "Vdouble"})

# ---- Event selection: ends at the nominal 4C kinematic fit ----
event_selection = Selection.new
event_selection
  .select_track {                  # exactly one positive + one negative track
    cos_theta 0.93                 # |cos(theta)| < 0.93
    Vz        10.0                 # |Vz| < 10 cm
    Vr        1.0                  # Vr < 1 cm
    nChrp     "==1"                # one positive track
    nChrn     "==1"                # one negative track
    nNet      "==0"                # net charge zero
  }
  .select_photon {                 # two photons from eta -> gamma gamma
    tdc_emc_start     0            # EMC timing 0-700 ns
    tdc_emc_end       14
    energyThreshold_b 0.025        # 25 MeV (barrel)
    energyThreshold_e 0.050        # 50 MeV (endcap)
    nGam              ">=2"        # at least two photons
  }
  .pid(method: :probability) {      # kaon ID: CL(K) > CL(pi), CL(K) > CL(p)
    prob_cut 0.001                 # probability cut
    identify :kaon, against: [:pion, :proton]
    nkp "==1"                      # one K+
    nkm "==1"                      # one K-
  }
  .kinematic_fit([:kp, :km, :gamma, :gamma]) {   # 4C fit to K+ K- gamma gamma
    nominal
    constrain_four_momentum
    chi2_cut 100
  }

# ---- Render the BOSS algorithm and execute on all datasets ----
etaphi_alg.with_decay_card(decay_card_signal).apply(event_selection)

root_files = etaphi_alg.execute_on(data_points + incMC_points + exMC_signal)