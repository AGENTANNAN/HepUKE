# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
data_3686 = DatasetManager.real_data.find("709_3686")        # ψ(2S) real data at 3.686 GeV (BOSS 709)
incMC_3686 = DatasetManager.inclusive_mc.find("709_3686")    # Corresponding inclusive MC sample

# Decay card for the signal process: ψ(2S) → π+π- J/ψ, J/ψ → γ η_c, η_c → γγ
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi+ pi- J/psi   PHSP;
    Enddecay

    Decay J/psi
    1.0000 gamma eta_c     PHSP;
    Enddecay

    Decay eta_c
    1.0000 gamma gamma     PHSP;
    Enddecay

    End
DECAYCARD

# Create the exclusive MC sample for the signal process (200k events)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3686_pipi_jpsi_gamma_etac"
  config.related_dataset = data_3686               # tie the MC to the real data conditions
  config.events          = 200000                  # 200k events requested
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "PiPiJpsiGammaEtac"
my_algorithm = Algorithm.new(alg_name)                                   # Algorithm for ψ(2S) → π+π- J/ψ(→ γ η_c(→ γγ))
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])                # header file for the algorithm
            .set_constant({"ECMS" => [:double, 3.686]})                  # ECMS = 3.686 GeV

# Build the event selection chain
event_selection = Selection.new
event_selection.select_track {          # Charged track selection
                  cos_theta 0.93        # |cosθ| < 0.93
                  Vz        10.0        # |Vz| < 10 cm along the beam axis
                  Vr        1.0         # Vr < 1 cm in the transverse plane
                  nChrp     "==1"       # exactly one positive track
                  nChrn     "==1"       # exactly one negative track
                  nNet      "==0"       # net charge zero
                }
               .select_photon {         # Photon selection
                  tdc_emc_start     0    # TDC start time
                  tdc_emc_end       14   # TDC end time
                  angle_to_track    10.0 # min angle to nearest charged track (degrees)
                  energyThreshold_b 0.025 # EMC barrel energy threshold (25 MeV)
                  energyThreshold_e 0.050 # EMC endcap energy threshold (50 MeV)
                  nGam              ">=3" # at least three photons
                }
               .pid(method: :probability) {   # PID by the probability method
                  prob_cut 0.001              # probability > 0.001
                  identify :pion, against: [:electron]  # π+ and π- identified against electrons
                  npip ">=1"                  # at least one π+
                  npim ">=1"                  # at least one π-
                }
               .kinematic_fit([:pip, :pim, :gamma, :gamma, :gamma]) {  # 4C fit on π+π-γγγ
                  nominal                  # mark as the nominal fit
                  constrain_four_momentum  # constrain total four-momentum to the CMS energy
                  chi2_cut 200             # loose χ² cut (tight χ² applied later in ROOT)
                }

# Generate the complete algorithm for the process defined in the decay card
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
# Execute on the real data, inclusive MC and signal exclusive MC
root_files = my_algorithm.execute_on([data_3686, incMC_3686, exMC_signal])