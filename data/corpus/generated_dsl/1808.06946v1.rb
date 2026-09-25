### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")        # J/psi real data at 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")     # Corresponding inclusive MC sample

# Decay card for the signal process (EvtGen format)
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma K_S0 K_S0 PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC: 1,000,000 events of J/psi -> gamma K_S K_S, K_S -> pi+ pi-
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "jpsi_to_gamma_ksks_exclusive_mc"
  config.related_dataset = jpsi_data
  config.events          = 1000000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "JpsiToGammaKsKs"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})       # CMS energy 3.097 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {                                # Charged track selection
                  cos_theta  0.93                             # |cos(theta)| < 0.93
                  Vz         10.0                             # |Vz| < 10 cm
                  Vr         1.0                              # Vr < 1 cm
                  nChrp      "==2"                            # Exactly 2 positive tracks
                  nChrn      "==2"                            # Exactly 2 negative tracks
                  nNet       "==0"                            # Net charge zero (4 tracks total)
                }
               .select_photon {                               # Photon selection
                  tdc_emc_start     0
                  tdc_emc_end       14
                  angle_to_track    10.0
                  energyThreshold_b 0.025                     # > 25 MeV in the barrel
                  energyThreshold_e 0.050                     # > 50 MeV in the endcap
                  nGam              ">=1"                     # At least one good photon
                }
               .assign({:chrgp => :pip, :chrgn => :pim})      # No PID; all tracks assumed pions
               .secondary_vertex_fit([:pip, :pim]) {          # First K_S from a pi+ pi- pair
                  build_virtual_particle(:K_S0).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list     # Do not reuse these tracks
                }
               .secondary_vertex_fit([:pip, :pim]) {          # Second K_S from the remaining pi+ pi- pair
                  build_virtual_particle(:K_S0).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
                }
               # 6C kinematic fit to gamma K_S K_S: 4-momentum conservation + one mass
               # constraint per K_S (4C + 2 x 1C = 6C)
               .kinematic_fit([:gamma, :K_S0, :K_S0]) {
                  nominal                                     # Nominal fit: corrected four-momenta used downstream
                  constrain_four_momentum                     # 4C energy-momentum conservation
                  invariant_mass_of(:K_S0).constrain_to_nominal_mass_of(:K_S0)  # 1C mass constraint for each K_S
                  chi2_cut 60                                 # chi^2 < 60
                }

# K_S flight-length significance cuts (L/sigma_L > 0 for each K_S and combined
# sqrt((L1/sigma1)^2 + (L2/sigma2)^2) > 2.2) have no dedicated DSL method; recorded as a note.
my_algorithm
  .note(:decay_length_significance, "each secondary-vertex-fitted K_S required to satisfy
    L/sigma_L > 0, and the two K_S together required to satisfy
    sqrt((L1/sigma1)^2 + (L2/sigma2)^2) > 2.2; applied to the K_S candidates built by
    the two secondary_vertex_fit steps before the 6C kinematic fit")
  .note(:background_veto, "inclusive MC (~0.5%) and continuum (~0.7%) backgrounds estimated
    to be negligible and are not rejected by an explicit BOSS-level veto; inclusive MC is
    still processed for efficiency/background bookkeeping")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

# Execute on real data, inclusive MC and the exclusive signal MC
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])