# ============================================================================
# e+e- -> pi0 pi0 psi(3686), psi(3686) -> pi+ pi- J/psi,
# J/psi -> e+e- / mu+mu-  (50% / 50%), pi0 -> gamma gamma
# XYZ scan region, c.m. energy 4.008 - 4.951 GeV
# ============================================================================

### Datasets ###
# Real data at the scanned energy points (BOSS_version_ECMS convention)
data_points = [
  DatasetManager.real_data.find("703_4009"),   # ~4.008 GeV
  DatasetManager.real_data.find("703_4180"),   # ~4.178 GeV
  DatasetManager.real_data.find("703_4230"),   # ~4.226 GeV
  DatasetManager.real_data.find("703_4260"),   # ~4.258 GeV
  DatasetManager.real_data.find("703_4360"),   # ~4.358 GeV
  DatasetManager.real_data.find("703_4420"),   # ~4.416 GeV
  DatasetManager.real_data.find("703_4600"),   # ~4.600 GeV
  DatasetManager.real_data.find("706_4680"),   # ~4.682 GeV
  DatasetManager.real_data.find("707_4946")    # ~4.951 GeV
]

# Inclusive MC at the same energy points
incMC_points = [
  DatasetManager.inclusive_mc.find("703_4009"),
  DatasetManager.inclusive_mc.find("703_4180"),
  DatasetManager.inclusive_mc.find("703_4230"),
  DatasetManager.inclusive_mc.find("703_4260"),
  DatasetManager.inclusive_mc.find("703_4360"),
  DatasetManager.inclusive_mc.find("703_4420"),
  DatasetManager.inclusive_mc.find("703_4600"),
  DatasetManager.inclusive_mc.find("706_4680"),
  DatasetManager.inclusive_mc.find("707_4946")
]

# Decay card (EvtGen): the description gives no explicit top mother, so the
# BESIII KKMC convention top mother psi(4260) is used.
# J/psi leptonic decays are fixed to 50% e+e- and 50% mu+mu-.
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 pi0 pi0 psi(2S) PHSP;
    Enddecay

    Decay psi(2S)
    1.0000 pi+ pi- J/psi PHSP;
    Enddecay

    Decay J/psi
    0.5000 e+ e-    PHOTOS VLL;
    0.5000 mu+ mu-  PHOTOS VLL;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 200k-event exclusive MC at each energy point (same card / cross section)
exMCs = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_xyz_pi0pi0psi2S_ll"
  config.events        = 200_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end
exMCs.each { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }

### Event selection (BOSS) ###
# Both J/psi decay modes (e+e- and mu+mu-) share identical final-state criteria,
# so a single Algorithm instance is used for all processes.
alg_name = "Pi0Pi0Psi2Sll"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({ "ECMS" => [:double, 4.260] })   # representative scan point
            .set_alias({ "std::vector<double>" => "Vdouble" })

event_selection = Selection.new
  .select_track {                 # Charged track selection
    cos_theta 0.93                # |cos(theta)| < 0.93
    Vz        10.0                # |Vz| < 10 cm
    Vr        1.0                 # Vr < 1 cm
    nChrp     "==2"               # exactly two positive tracks
    nChrn     "==2"               # exactly two negative tracks
    nNet      "==0"               # net charge zero
  }
  .select_photon {                # Photon selection
    tdc_emc_start     0           # TDC window 0 - 700 (DSL units)
    tdc_emc_end       14
    angle_to_track    10.0        # >= 10 deg away from any charged track
    energyThreshold_b 0.025       # 25 MeV barrel
    energyThreshold_e 0.050       # 50 MeV endcap
    nGam              ">=4"       # at least four photons
  }
  .pid(method: :probability) {    # PID - probability method
    prob_cut 0.001
    # Tracks with p > 1.1 GeV/c treated as leptons; e/mu separated via EMC energy
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.1,
                                    treat_as_electron_if_energy_above: 0.45
    identify :pion, against: [:kaon, :proton]   # pi+ and pi- vs K / p
    npip "==1"                    # one pi+
    npim "==1"                    # one pi-
    nlp  "==1"                    # one positive lepton
    nlm  "==1"                    # one negative lepton
  }
  # Photons paired into two pi0 candidates via a mass-constrained (Kalman) fit;
  # the pairing minimises (M(gg) - m_pi0)^2.
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=2"                    # two pi0 candidates from the photons
  }
  # Nominal 4C kinematic fit on gamma gamma gamma gamma pi+ pi- l+ l-.
  # When more than four photons are present the best (smallest chi2) combination is kept.
  .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :pip, :pim, :lp, :lm]) {
    nominal
    constrain_four_momentum
    chi2_cut 120
  }
  # 7C refit: 4C + two pi0 mass constraints + J/psi mass constraint,
  # on the same track indices as the nominal fit.
  .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :pip, :pim, :lp, :lm]) {
    use_track_index_from_nominal_kmfit
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)
    chi2_cut 200
  }

my_algorithm
  .note(:pid_correction_method,
        "charged-track PID is non-standard: tracks with p>1.1 GeV/c are treated as
         leptons and tracks with p<0.75 GeV/c are assigned to pions; e/mu separation
         uses EMC E/p>0.7 for electrons and EMC E<0.45 GeV for muons.
         identify_high_momentum_leptons exposes only a momentum threshold and a single
         EMC-energy threshold, so the E/p and soft-pion criteria are implemented in the
         generated BOSS code.")
  .note(:pi0_pairing,
        "the two pi0 are formed from the selected photons by minimising (M(gg)-m_pi0)^2;
         each pair is required to satisfy |M(gg)-m_pi0| < 20 MeV/c^2; when more than four
         photons are present the combination with the smallest 4C chi2 is retained
         (default best-combination selection of the kinematic fit).")
  .note(:energy_scan_ecms,
        "ECMS varies over the 4.008-4.951 GeV scan; the value 4.260 GeV set here is
         representative and the correct per-point c.m. energy must be supplied for the
         4C/7C kinematic fits at each energy point.")

my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_algorithm.execute_on(data_points + incMC_points + exMCs)