### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")     # psi(2S) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # corresponding inclusive MC

# Decay cards: one per chi_cJ mode (J = 0, 1, 2). All three share the same final
# state p K- anti-Lambda eta, differing only in the intermediate chi_cJ.
decay_card_chic0 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c0 P2GC0;
    Enddecay

    Decay chi_c0
    1.000 p+ K- anti-Lambda0 eta PHSP;
    Enddecay

    Decay anti-Lambda0
    1.000 anti-p- pi+ HypWK;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

decay_card_chic1 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c1 P2GC1;
    Enddecay

    Decay chi_c1
    1.000 p+ K- anti-Lambda0 eta PHSP;
    Enddecay

    Decay anti-Lambda0
    1.000 anti-p- pi+ HypWK;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

decay_card_chic2 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c2 P2GC2;
    Enddecay

    Decay chi_c2
    1.000 p+ K- anti-Lambda0 eta PHSP;
    Enddecay

    Decay anti-Lambda0
    1.000 anti-p- pi+ HypWK;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC: 500k events for each of the three chi_cJ modes
exMC_chic0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_chic0_pKLambdaEta"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_chic0
  config.cross_section   = :default
end

exMC_chic1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_chic1_pKLambdaEta"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_chic1
  config.cross_section   = :default
end

exMC_chic2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_chic2_pKLambdaEta"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_chic2
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# All three chi_cJ modes share identical final state and selection chain -> one Algorithm.
alg_name = "ChicJpKLambdaEta"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})   # CMS energy 3.686 GeV

event_selection = Selection.new
event_selection.select_track {          # Charged track selection
                  cos_theta 0.93        # |cos(theta)| < 0.93
                  Vz        10.0        # |Vz| < 10 cm
                  Vr        1.0         # Vr < 1 cm
                  nChrp     ">=2"       # at least 2 positive tracks
                  nChrn     ">=2"       # at least 2 negative tracks
                  nNet      "==0"       # net charge zero
                }
               .select_photon {         # Photon selection
                  tdc_emc_start     0
                  tdc_emc_end       14      # EMC time window 0-14 (<=700 ns)
                  angle_to_track    10.0    # >=10 degrees from any charged track
                  energyThreshold_b 0.025   # >25 MeV in the barrel
                  energyThreshold_e 0.050   # >50 MeV in the endcap
                  nGam              ">=3"   # at least 3 photon candidates
                }
               .pid(method: :probability) {   # PID by probability method
                  prob_cut 0.001              # prob > 0.001
                  identify :proton, against: [:kaon, :pion]   # p+/pbar vs K, pi
                  identify :kaon,   against: [:proton, :pion] # K+/K- vs p, pi
                  nprp ">=1"                  # at least one proton
                  nkm  ">=1"                  # at least one K-
                }
               .remove([:prp <= :chrgp, :km <= :chrgn])  # drop identified p and K- from the charged lists used for anti-Lambda
               .select_isolated_photon {  # isolated photon, >=20 deg from the (anti-)proton track
                  angle_to_prm_track 20.0
                  nGam ">=3"
                }
               .assign({:chrgp => :pip, :chrgn => :prm})  # remaining positive -> pi+, negative -> anti-p
               .secondary_vertex_fit([:prm, :pip]) {      # reconstruct anti-Lambda -> anti-p pi+
                  build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
                }
               .kinematic_fit([:prp, :km, :Lambda_bar, :gamma, :gamma, :gamma]) {  # nominal 4C fit on p K- anti-Lambda + 3 photons
                  nominal
                  constrain_four_momentum
                  chi2_cut 200            # loose chi2 cut; tight chi2_4C < 35 applied in ROOT
                }

# Any decay card of the three works for header generation (identical final state).
my_algorithm.with_decay_card(decay_card_chic1).apply(event_selection)

# Execute on real data, inclusive MC, and the three signal exclusive MC samples
root_files = my_algorithm.execute_on([psip_data, psip_incMC, exMC_chic0, exMC_chic1, exMC_chic2])