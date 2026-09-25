### Dataset description ###
# Full J/psi real data set at sqrt(s) = 3.097 GeV and the matching inclusive MC
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Decay card for the signal process (EvtGen format, EvtGen particle names):
# J/psi -> Xi0 anti-Xi0 ; anti-Xi0 -> anti-Lambda0 pi0 ; Xi0 -> Lambda0 Lambda0
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 Xi0 anti-Xi0 PHSP;
    Enddecay

    Decay anti-Xi0
    1.0000 anti-Lambda0 pi0 PHSP;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    Decay Xi0
    1.0000 Lambda0 Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for the signal: 2.25M events of J/psi -> Xi0 anti-Xi0 -> (Lambda Lambda)(anti-Lambda pi0)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_to_xi0_xibar0_ll_lbarpi0"
  config.related_dataset = jpsi_data
  config.events          = 2_250_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "JpsiToXi0Xibar0"
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {            # Charged track selection
                  cos_theta 0.93          # |cos(theta)| < 0.93
                  Vz        10.0          # |Vz| < 10 cm
                  Vr        1.0           # |Vr| < 1 cm
                  nChrp     "==3"         # exactly three positive tracks
                  nChrn     "==3"         # exactly three negative tracks
                  nNet      "==0"         # net charge zero
                }
               .select_photon {           # Photon selection
                  tdc_emc_start    0      # EMC time window 0-700 ns
                  tdc_emc_end      14
                  angle_to_track   10.0   # angle to nearest charged track > 10 degrees
                  energyThreshold_b 0.025 # E > 25 MeV (barrel)
                  energyThreshold_e 0.050 # E > 50 MeV (endcap)
                  nGam             ">=2"  # at least two photons (for pi0 -> gamma gamma)
                }
               .pid(method: :probability) {   # Probability-based PID
                  prob_cut 0.001              # PID probability > 0.001
                  identify :proton, against: [:kaon, :pion]   # protons (largest likelihood vs K, pi); identifies p and anti-p
                  identify :pion,   against: [:kaon, :proton] # pions  (largest likelihood vs K, p)
                  nprp "==2"                  # 2 protons   (both Lambda -> p pi-)
                  nprm "==1"                  # 1 anti-proton (anti-Lambda -> anti-p pi+)
                  npip "==1"                  # 1 pi+
                  npim "==2"                  # 2 pi-
                }
               .secondary_vertex_fit([:prm, :pip]) {  # Reconstruct the tag anti-Lambda from anti-p pi+
                  build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list        # used tracks removed
                }
               .kalman_kinematic_fit([:gamma, :gamma]) {  # Kalman 1C fit of gamma gamma to the pi0 mass
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                  chi2_cut 25                             # chi2 < 25
                  npi0 ">=1"                              # at least one pi0
                }
               .secondary_vertex_fit([:prp, :pim]) {      # Form Lambda candidates from all remaining p pi- pairs
                  build_virtual_particle(:Lambda).by_minimizing_mass_difference
                }
               .kinematic_fit([:Lambda_bar, :pi0, :Lambda, :Lambda]) {  # Nominal 4C fit: J/psi -> anti-Lambda pi0 Lambda Lambda
                  nominal
                  constrain_four_momentum                 # constrain total four-momentum to CMS energy
                  chi2_cut 200                            # loose cut in BOSS; tight cut applied in ROOT
                }

my_Algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_Algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])