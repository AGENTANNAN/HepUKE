### Dataset description ###
# J/psi (3.097 GeV) real data and inclusive MC; sample name convention 708_3097
jpsi_data = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Decay card for the signal: J/psi -> phi K_S0 K_S0, phi -> K+ K-, (tagged) K_S0 -> pi+ pi-.
# The other K_S0 is the invisible signal candidate and is left undecayed (alias keeps it distinct).
decay_card_signal = <<~DECAYCARD
    Alias K_S0_vis K_S0

    Decay J/psi
    1.0 phi K_S0_vis K_S0 PHSP;
    Enddecay

    Decay phi
    1.0 K+ K- VSS;
    Enddecay

    Decay K_S0_vis
    1.0 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for the signal process (100k events)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "jpsi_phi_ksks_invisible"
  config.related_dataset = jpsi_data
  config.events = 100000
  config.decay_card = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "JpsiPhiKsKs"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})

event_selection = Selection.new
  .select_track {                 # Charged track selection
      cos_theta 0.93              # |cos(theta)| < 0.93
      Vz 100.0                    # |Vz| < 100 cm
      Vr 10.0                     # Vr < 10 cm
      nChrp "==2"                 # exactly 2 positive tracks (K+, pi+)
      nChrn "==2"                 # exactly 2 negative tracks (K-, pi-)
      nNet "==0"                  # exactly four charged tracks, net charge 0
  }
  .select_photon {                # Photon selection (no minimum multiplicity)
      tdc_emc_start 0             # EMC TDC start
      tdc_emc_end 14              # EMC TDC end
      angle_to_track 10.0         # angle to nearest charged track > 10 degrees
      energyThreshold_b 0.025     # 25 MeV in the barrel
      energyThreshold_e 0.050     # 50 MeV in the endcap
  }
  .pid(method: :probability) {    # PID by probability method, K/pi separation
      prob_cut 0.001
      identify :kaon, against: [:pion, :proton]   # K+ and K- (charge-conjugation shorthand)
      identify :pion, against: [:kaon, :proton]   # pi+ and pi- (charge-conjugation shorthand)
      nkp "==1"                   # exactly 1 K+
      nkm "==1"                   # exactly 1 K-
      npip "==1"                  # exactly 1 pi+
      npim "==1"                  # exactly 1 pi-
  }
  .secondary_vertex_fit([:pip, :pim]) {           # Reconstruct the tagged K_S0 from pi+ pi-
      build_virtual_particle(:K_S0).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list    # remove the pions used in the fit
  }
  .kinematic_fit([:kp, :km, :K_S0]) {             # Fit K+ K- and the tagged K_S0 ...
      nominal                                     # ... flagged as the nominal fit
      miss_track_of :K_S0                         # the other K_S0 treated as a missing particle
      constrain_four_momentum                     # four-momentum constraint
      chi2_cut 200                                # chi^2 < 200
  }

my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])