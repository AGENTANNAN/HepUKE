# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# Nine energy points of the 4.180 - 4.280 GeV e+e- -> omega chi_c0 scan (BOSS 703 samples)
scan_samples = %w[703_4180 703_4190 703_4200 703_4210 703_4220 703_4230 703_4246 703_4260 703_4280]

data_points  = scan_samples.map { |s| DatasetManager.real_data.find(s) }   # real data, one dataset per point
incMC_points = scan_samples.map { |s| DatasetManager.inclusive_mc.find(s) } # corresponding inclusive MC

# Decay card: e+e- -> omega chi_c0 with chi_c0 -> pi+ pi-
decay_card_pipi = <<~DECAYCARD
    Decay psi(4260)
    1.0000 omega chi_c0 PHSP;
    Enddecay

    Decay chi_c0
    1.0000 pi+ pi- PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Decay card: e+e- -> omega chi_c0 with chi_c0 -> K+ K-
decay_card_KK = <<~DECAYCARD
    Decay psi(4260)
    1.0000 omega chi_c0 PHSP;
    Enddecay

    Decay chi_c0
    1.0000 K+ K- PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# KKMC-generated exclusive MC for both chi_c0 decay modes: one sample per energy point
exMCs_pipi = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_omegachic0_pipi"   # becomes ..._703_4180, ..._703_4190, ...
  config.events        = 100000
  config.decay_card    = decay_card_pipi
  config.cross_section = :default
end

exMCs_KK = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_omegachic0_KK"     # becomes ..._703_4180, ..._703_4190, ...
  config.events        = 100000
  config.decay_card    = decay_card_KK
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "OmegaChiC0"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.260]})
            .note(:measured_ecms, "center-of-mass energy is taken run-by-run from the measured value (MeasuredEcmsSvc) instead of a fixed ECMS constant over the 4.180-4.280 GeV scan; the ECMS constant is only a fallback")

# Single shared selection chain for both chi_c0 decay hypotheses
event_selection = Selection.new
event_selection.select_track {          # charged track selection
                  cos_theta   0.93      # |cos(theta)| < 0.93
                  Vz          10.0      # |Vz| < 10 cm
                  Vr          1.0       # Vr < 1 cm
                  nChrp       "==2"     # two positively charged tracks
                  nChrn       "==2"     # two negatively charged tracks
                  nNet        "==0"     # net charge zero
                }
               .select_photon {         # photon selection
                  tdc_emc_start     0        # EMC timing window 0-700 ns
                  tdc_emc_end       14
                  angle_to_track    10.0     # at least 10 deg separation from any charged track
                  energyThreshold_b 0.025    # 25 MeV in the barrel
                  energyThreshold_e 0.050    # 50 MeV in the endcap
                  nGam              ">=2"    # two photons for pi0 -> gamma gamma
                }
               # gamma gamma Kalman fit: constrain the two photons to the nominal pi0 mass
               .kalman_kinematic_fit([:gamma, :gamma]) {
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                  chi2_cut 25
                  npi0     ">=1"
                }
               # chi_c0 -> pi+ pi- hypothesis: all four charged tracks (chi_c0 and omega
               # daughters) are pions; 5C fit = 4-momentum conservation + pi0 mass constraint
               .assign({:chrgp => :pip, :chrgn => :pim})
               .kinematic_fit([:pip, :pip, :pim, :pim, :pi0]) {
                  nominal                  # nominal fit, its four-momenta are stored
                  constrain_four_momentum  # 4C energy-momentum conservation (+ pi0 mass from the Kalman fit -> 5C)
                  chi2_cut 100
                }
               # chi_c0 -> K+ K- hypothesis: kaons and pions from chi_c0 and omega are clearly
               # separated in momentum, so p > 1 GeV/c tracks are the chi_c0 daughters (K) and
               # p < 1 GeV/c tracks are the omega daughters (pions)
               .assign({:chrgp => :kp, :chrgn => :km})
               .remove(:kp)  { condition "three_momentum_of(:kp)  < 1.0" }  # keep only high-momentum K+
               .remove(:km)  { condition "three_momentum_of(:km)  < 1.0" }  # keep only high-momentum K-
               .remove(:pip) { condition "three_momentum_of(:pip) > 1.0" }  # keep only low-momentum pi+
               .remove(:pim) { condition "three_momentum_of(:pim) > 1.0" }  # keep only low-momentum pi-
               .kinematic_fit([:kp, :km, :pip, :pim, :pi0]) {
                  # competing hypothesis: chi2 is stored (no chi2_cut, not nominal); the choice
                  # of the smaller 5C chi2 between pi+pi- and K+K- is made in the ROOT stage
                  constrain_four_momentum
                }

my_algorithm.with_decay_card(decay_card_pipi).apply(event_selection)

# Execute on real data, inclusive MC and both exclusive MC sets
datasets = data_points + incMC_points + exMCs_pipi + exMCs_KK
root_files = my_algorithm.execute_on(datasets)