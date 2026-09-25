### Dataset description ###
psip_data = DatasetManager.real_data.find("709_3686")        # ψ(3686) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")    # corresponding inclusive MC sample

# Decay card for the signal process ψ(3686) → γ χ_cJ, χ_cJ → φφ (phase space), φ → K+K−
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c1 P2GC1;
    Enddecay

    Decay chi_c1
    1.000 phi phi PHSP;
    Enddecay

    Decay phi
    1.000 K+ K- VSS;
    Enddecay

    End
DECAYCARD

# Exclusive MC for the signal process: 1,000,000 events
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3686_gammaChicJ_phiphi"  # signal exclusive MC
  config.related_dataset = psip_data                      # associated real dataset
  config.events          = 1_000_000                      # 1,000,000 events
  config.decay_card      = decay_card_signal              # decay card above
  config.cross_section   = :default                       # default cross section
end
# Save the MC configuration for later use (optional)
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) ###
alg_name = "GammaChicJToPhiPhi"  # ψ(3686) → γ χ_cJ, χ_cJ → φφ, φ → K+K−
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]}) # CMS energy 3.686 GeV

event_selection = Selection.new
event_selection.select_track {              # Charged track selection
                  cos_theta 0.93            # |cosθ| < 0.93
                  Vz        10.0            # |Vz| < 10 cm
                  Vr        1.0             # Vr < 1 cm
                  nChrp     "==2"           # exactly two positive tracks (K+)
                  nChrn     "==2"           # exactly two negative tracks (K−)
                  nNet      "==0"           # net charge zero
                }
               .select_photon {             # Photon selection
                  tdc_emc_start      0      # EMC timing window start (0)
                  tdc_emc_end        14     # EMC timing window end (700 ns)
                  angle_to_track     10.0   # opening angle > 10° from any charged track
                  energyThreshold_b  0.025  # E > 25 MeV in the barrel
                  energyThreshold_e  0.050  # E > 50 MeV in the endcap
                  nGam               ">=1"  # at least one photon
                }
               .pid(method: :probability) { # Kaon identification (probability method)
                  prob_cut 0.001            # PID probability > 0.001
                  identify :kaon, against: [:pion, :proton]  # K+ and K− vs π and p
                  nkp "==2"                 # exactly two K+
                  nkm "==2"                 # exactly two K−
                }
               .kinematic_fit([:gamma, :kp, :km, :kp, :km]) {  # 4C fit to γ K+K−K+K−
                  nominal                   # nominal fit — its four-momenta are used
                  constrain_four_momentum   # energy–momentum conservation
                  chi2_cut 60               # χ² < 60
                }

# Generate the complete algorithm for the process in the decay card
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
# Execute on real data, inclusive MC and signal exclusive MC
root_files = my_algorithm.execute_on([psip_data, psip_incMC, exMC_signal])