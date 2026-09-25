# === Dataset preparation ===
# 4.178 GeV reference point of the η h_c scan
data_4180  = DatasetManager.real_data.find("703_4180")      # 4.178 GeV real data
incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")   # matching inclusive MC

# Decay card: e+e- -> eta h_c, h_c -> gamma eta_c, eta_c -> p pbar, eta -> gamma gamma
decay_card_signal = <<~DECAYCARD
  Decay psi(4260)
  1.000 eta h_c PHSP;
  Enddecay

  Decay h_c
  1.000 gamma eta_c PHSP;
  Enddecay

  Decay eta_c
  1.000 p+ anti-p- PHSP;
  Enddecay

  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# 22 scan points spanning 4.129 - 4.600 GeV (same signal MC for every energy point)
scan_points = [
  "705_4130", "705_4160", "703_4180", "703_4190", "703_4200",
  "703_4210", "703_4220", "703_4230", "703_4237", "703_4246",
  "703_4260", "703_4270", "703_4280", "705_4290", "705_4315",
  "705_4340", "703_4360", "705_4380", "705_4400", "703_4420",
  "705_4440", "703_4600"
].map { |name| DatasetManager.real_data.find(name) }

# 500k exclusive MC events for eta_c -> p pbar, one sample per energy point
exMC_scan = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "eta_hc_etac_ppbar"
  config.events        = 500_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

# === Event selection (BOSS) ===
alg_name = "EtaHcEtacPPbar"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.178]})   # 4.178 GeV centre-of-mass energy
            .note(:candidate_selection,
                  "when several eta candidates fall in the loose hc recoil-mass " \
                  "window [3.480, 3.600] GeV/c^2, the eta_c (p pbar) candidate " \
                  "closest to the nominal eta_c mass is retained")

event_selection = Selection.new
  .select_track {                 # exactly one positive and one negative track, net charge 0
    cos_theta 0.93                # |cos(theta)| < 0.93
    Vz        10.0                # |Vz| < 10 cm
    Vr        1.0                 # Vr < 1 cm
    nChrp     "==1"               # exactly one positive track
    nChrn     "==1"               # exactly one negative track
    nNet      "==0"               # net charge zero
  }
  .select_photon {                # >= 3 photons (two from eta -> gamma gamma + hc -> gamma eta_c)
    tdc_emc_start     0           # TDC start time
    tdc_emc_end       14          # TDC end time
    angle_to_track    10.0        # at least 10 deg from any charged track
    energyThreshold_b 0.025       # 25 MeV barrel threshold
    energyThreshold_e 0.050       # 50 MeV endcap threshold
    nGam              ">=3"       # at least three photons
  }
  .pid(method: :probability) {    # probability PID, p / pbar separated from K and pi
    prob_cut 0.001                # PID probability > 0.001
    identify :proton, against: [:kaon, :pion]
    nprp     "==1"                # exactly one proton
    nprm     "==1"                # exactly one anti-proton
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # reconstruct eta from gamma gamma (1C)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25                   # chi2 < 25
    neta     "==1"                # exactly one eta candidate
  }
  .kinematic_fit([:eta, :gamma, :prp, :prm]) {   # nominal 4C fit on eta gamma p pbar
    nominal
    constrain_four_momentum
    chi2_cut 25                   # chi2 < 25
  }

my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on the 4.178 GeV data + inclusive MC and the per-energy-point signal MC
root_files = my_algorithm.execute_on([data_4180, incMC_4180] + exMC_scan)