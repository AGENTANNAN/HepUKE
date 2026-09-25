### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi (3.097 GeV) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # Corresponding inclusive MC sample

# Decay card for the signal chain: J/psi -> phi eta, phi -> K+K-, eta -> pi0 + invisible
# (the invisible eta decay is simulated as eta -> pi0 nu anti-nu), pi0 -> gamma gamma
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 phi eta PHSP;
    Enddecay

    Decay phi
    1.0000 K+ K- VSS;
    Enddecay

    Decay eta
    1.0000 pi0 nu anti-nu PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 500k-event exclusive MC for the signal chain
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_phi_eta_invisible_mc"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "Jpsi2PhiEtaInv"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})   # CMS energy of J/psi
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {                    # Charged track selection
                  cos_theta 0.93                  # |cos(theta)| < 0.93
                  Vz        10.0                  # |Vz| < 10 cm
                  Vr        1.0                   # Vr < 1 cm
                  nChrp     "==1"                 # Exactly one positive track
                  nChrn     "==1"                 # Exactly one negative track  -> exactly two tracks
                  nNet      "==0"                 # Net charge zero
                }
               .select_photon {                   # Photon selection
                  tdc_emc_start     0             # TDC EMC start
                  tdc_emc_end       14            # TDC EMC end
                  angle_to_track    10.0          # Reject photons within 10 degrees of a track
                  energyThreshold_b 0.025         # 25 MeV in the EMC barrel
                  energyThreshold_e 0.050         # 50 MeV in the EMC endcap
                  nGam              ">=2"         # At least two photons
                }
               .pid(method: :probability) {       # Kaon PID with the probability method
                  prob_cut 0.001                  # Kaon likelihood > 0
                  identify :kaon, against: [:pion, :proton]  # K+ and K- vs pi and p
                  nkp ">=1"                       # At least one K+
                  nkm ">=1"                       # At least one K-
                }
               .kinematic_fit([:kp, :km, :gamma, :gamma]) {  # 4C kinematic fit of K+K-gammagamma
                  nominal                         # Nominal fit
                  constrain_four_momentum         # 4-momentum conservation (4C)
                  invariant_mass_of(:kp, :km).within(1.00, 1.04)      # M(K+K-) in [1.00, 1.04] GeV (phi)
                  invariant_mass_of(:gamma, :gamma).within(0.115, 0.150)  # M(gammagamma) in [0.115, 0.150] GeV
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # gamma gamma constrained to pi0 mass
                  chi2_cut 200                    # Loose chi^2 < 200 (tight cut applied in ROOT)
                }

# Generate the complete algorithm for the process defined in the decay card
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
# Execute on the real data, inclusive MC and exclusive signal MC
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])