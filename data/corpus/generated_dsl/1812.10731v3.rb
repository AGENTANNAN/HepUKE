# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
data_4600 = DatasetManager.real_data.find("703_4600")        # 4.600 GeV real data (567 pb^-1)
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")    # matching inclusive MC sample

# Decay card for the signal process (EvtGen format, EvtGen particle names)
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.000 Lambda_c+ anti-Lambda_c- PHSP;
    Enddecay

    Decay Lambda_c+
    1.000 Lambda0 eta pi+ PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.000 anti-p- K+ pi- PHSP;
    Enddecay

    Decay Lambda0
    1.000 p+ pi- PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Create the exclusive MC sample for the signal process (200,000 events)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_4600_LcToLambdaEtaPi"
  config.related_dataset = data_4600
  config.events          = 200000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "LambdacToLambdaEtaPi"
my_alg = Algorithm.new(alg_name)
my_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
       .set_constant({"ECMS" => [:double, 4.600]})
       .note(:lambda_daughter_track_cuts, "the Λ daughter p and π- are allowed a looser charged-track quality than the global selection: |Vz| < 20 cm and no Vr cut (global selection: |cosθ| < 0.93, |Vz| < 10 cm, Vr < 1 cm)")
       .note(:lambda_selection, "Λ candidates built from pπ- required to have invariant mass in [1.111, 1.121] GeV/c2, secondary vertex fit chi2 < 100, and flight distance significance > 2 sigma; these criteria are applied in BOSS outside the DSL secondary_vertex_fit primitive")

event_selection = Selection.new
event_selection.select_track {          # Charged track selection
                  cos_theta   0.93      # |cosθ| < 0.93
                  Vz          10.0      # |Vz| < 10 cm
                  Vr          1.0       # Vr < 1 cm
                  nChrp       ">=2"     # at least two positive tracks
                  nChrn       ">=2"     # at least two negative tracks
                  nNet        "==0"     # net charge zero
                }
               .select_photon {          # Photon selection
                  tdc_emc_start     0     # TDC start time
                  tdc_emc_end       14    # TDC end time
                  angle_to_track    10.0  # angle to any charged track > 10 degrees
                  energyThreshold_b 0.025 # EMC barrel threshold 25 MeV
                  energyThreshold_e 0.050 # EMC endcap threshold 50 MeV
                  nGam              ">=2" # at least two photon candidates
                }
               .pid(method: :probability) {   # PID with the probability method
                  prob_cut 0.001              # PID probability > 0.001
                  identify :proton, against: [:kaon, :pion]  # p+ and anti-p- vs K and π
                  identify :pion,   against: [:kaon]         # π+ and π- vs K
                  nprp ">=1"                                 # at least one proton; anti-proton count unconstrained
                }
               .kalman_kinematic_fit([:gamma, :gamma]) {   # Reconstruct η from two photons (1-C mass constraint)
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
                  chi2_cut 25             # χ² < 25
                  neta ">=1"              # at least one η candidate
                }
               .secondary_vertex_fit([:prp, :pim]) {       # Build Λ from pπ- via secondary vertex fit
                  build_virtual_particle(:Lambda).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
                }
               .kinematic_fit([:pip, :Lambda, :eta]) {     # Final 4C kinematic fit to π+Λη
                  nominal                 # nominal fit — corrected four-momenta are kept
                  constrain_four_momentum # 4C energy-momentum constraint
                  chi2_cut 200            # χ² < 200
                }

# Generate the algorithm for the process in the decay card and execute
my_alg.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_alg.execute_on([data_4600, incMC_4600, exMC_signal])