# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")      # psi(2S) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")    # corresponding inclusive MC sample

# Decay card for the signal process psi(2S) -> pi0 h_c, h_c -> gamma eta_c, pi0 -> gamma gamma
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi0 h_c PHSP;
    Enddecay

    Decay h_c
    1.0000 gamma eta_c PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 300k exclusive MC events for the signal process
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_pi0_hc_etac"
  config.related_dataset = psip_data       # associated real dataset
  config.events          = 300000          # 300k signal events
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) ###
alg_name = "Pi0HcEtac"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})      # sqrt(s) = 3.686 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {               # Charged track selection (no explicit PID is applied)
                  cos_theta 0.93             # |cos(theta)| < 0.93
                  Vz        10.0             # |Vz| < 10 cm
                  Vr        1.0              # Vr < 1 cm
                  nTot      ">=2"            # at least two charged tracks in total
                }
               .select_photon {              # Photon selection
                  tdc_emc_start     0        # EMC TDC start
                  tdc_emc_end       14       # EMC TDC end
                  angle_to_track    10.0     # angle to the nearest charged track > 10 degrees
                  energyThreshold_b 0.025    # barrel energy threshold 25 MeV
                  energyThreshold_e 0.050    # endcap energy threshold 50 MeV
                  nGam ">=2"                 # at least two photons
                }
               .kalman_kinematic_fit([:gamma, :gamma]) {   # 1C Kalman fit reconstructing pi0 from photon pairs
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                  chi2_cut 200               # 1C fit chi2 < 200
                  npi0 ">=1"                 # require at least one pi0 candidate
                }
               .kinematic_fit([:pi0, :gamma]) {            # final signal-extraction fit (tagged h_c channel: pi0 + E1 gamma)
                  nominal
                  constrain_four_momentum
                  chi2_cut 200
                }

# BOSS-side procedures that have no dedicated DSL construct
my_algorithm
  .note(:e1_photon_selection, "In the tagged h_c channel an additional E1 photon is required with 0.465 < E_gamma < 0.535 GeV that does not combine with any other photon to form a pi0; the per-photon energy window together with the anti-pi0 pairing veto has no dedicated DSL method and is applied as an event-level requirement")
  .note(:pi0_candidate_selection, "If several pi0 candidates fall in the recoil-mass signal window 3.500-3.550 GeV/c^2 (M_recoil(pi0) = M(h_c)), the candidate with the smallest 1C Kalman-fit chi2 is retained")
  .note(:background_veto, "psi(2S) -> J/psi pi+pi- vetoed by requiring the pi+pi- recoil mass to lie outside M_J/psi +/- 0.004 GeV/c^2; psi(2S) -> J/psi pi0pi0 vetoed by requiring the pi0pi0 recoil mass to lie outside [M_J/psi - 0.008, M_J/psi + 0.038] GeV/c^2")
  .note(:emc_total_energy, "Total EMC energy required to satisfy 0.6 < E_EMC < 3.2 GeV to suppress Bhabha and continuum backgrounds")
  .note(:signal_extraction_fit, "The pi0 recoil mass is computed from four-momentum conservation and used for the final signal-extraction fit")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

# Execute on real data, inclusive MC and the signal exclusive MC
root_files = my_algorithm.execute_on([psip_data, psip_incMC, exMC_signal])