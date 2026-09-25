### Dataset description ###
# 26-point J/psi scan between 3000 and 3119.88 MeV (accompanying PRD version;
# same event selection as the Letter). The DSL enumerates the subset of scan
# points resolvable from the current table and captures the full scope via a note.
scan_data_points = [
  DatasetManager.real_data.find("713_Rscan_3080"),
  DatasetManager.real_data.find("713_Rscan_3020"),
  DatasetManager.real_data.find("713_Rscan_3000"),
  DatasetManager.real_data.find("713_Rscan_2981"),
  DatasetManager.real_data.find("713_Rscan_2950")
]
inc_mc_points = scan_data_points.map { |d| DatasetManager.inclusive_mc.find("713_#{d.sample_name}") }

# ------------------------------------------------------------------
# Signal decay card: e+e- -> K_S K+ pi- via ConExc mode index 9.
# 'Particle vpho' omitted so that the DSL will inject the correct per-scan-point sqrt(s).
# ------------------------------------------------------------------
decay_card_ks_k_pi = <<~DECAYCARD
    Decay vpho
    1 ConExc 9;
    Enddecay
    Decay vhdr
    1 K_S0 K+ pi- PHSP;
    Enddecay
    Decay K_S0
    1 pi+ pi- PHSP;
    Enddecay
    End
DECAYCARD

sig_mc = DatasetManager.create_exclusive_mc_for(scan_data_points) do |config|
  config.sample_name   = "ee_KsKpi_ConExc"
  config.events        = 200000
  config.decay_card    = decay_card_ks_k_pi
  config.cross_section = :default
end
sig_mc.each { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }

########################################################################
# Algorithm: e+e- -> K_S0 K+ pi-  (with K_S0 -> pi+ pi-)
########################################################################
alg = Algorithm.new("EEtoKsKPi_PRD")
alg.set_header(["EEtoKsKPi_PRDAlg/EEtoKsKPi_PRD.h"])
   .set_constant({"ECMS" => [:double, 3.080]})   # nominal; per-run value used at runtime

sel = Selection.new
sel.select_track {                       # Non-K_S0 tracks: |Vz| < 10 cm, |Vxy| < 1 cm; |cos(theta)| < 0.93
      cos_theta 0.93
      Vz        10.0
      Vr        1.0
      nChrp     ">=2"                    # K+ from signal and pi+ from K_S0
      nChrn     ">=2"                    # pi- from signal and pi- from K_S0
      nNet      "==0"
    }
   .select_photon {
      tdc_emc_start      0
      tdc_emc_end        14
      angle_to_track     10.0
      energyThreshold_b  0.025
      energyThreshold_e  0.050
    }
   .pid(method: :probability) {          # K/pi ID from combined dE/dx + TOF likelihoods
      prob_cut 0.001
      identify :kaon, against: [:pion, :proton]
      identify :pion, against: [:kaon, :proton]
    }
   # Reconstruct K_S0 -> pi+ pi-: common vertex fit with chi2<100, |M(pi+pi-) - m_KS0| < 15 MeV/c^2,
   # decay length > 2*sigma. Longest-decay-length pick made in ROOT for multi-KS events.
   .secondary_vertex_fit([:pip, :pim]) {
      build_virtual_particle(:K_S0).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
    }
   # 4C kinematic fit under e+e- -> K_S0 K+ pi-
   .kinematic_fit([:K_S0, :kp, :pim]) {
      nominal
      constrain_four_momentum
      chi2_cut 200
    }
   # 5C fit: 4C + M(pi+ pi-) constrained to K_S0 mass, for improved kinematic-variable resolution
   .kinematic_fit([:K_S0, :kp, :pim]) {
      constrain_four_momentum
      invariant_mass_of(:K_S0).constrain_to_nominal_mass_of(:K_S0)
    }

alg.note(:dataset_scan_points,
         "Full analysis uses 26 e+e- energy points from 3000.00 to 3119.88 MeV (440.7 pb^-1); " \
         "the DSL enumerates the subset resolvable in the current 713 Rscan_29xx/30xx table. " \
         "Add the remaining points to scan_data_points once the table is completed.")
   .note(:ks_selection,
         "K_S0 -> pi+ pi- with the two tracks satisfying |Vz| < 20 cm (no PID); common-vertex fit " \
         "chi^2 < 100; |M(pi+pi-) - m_KS0| < 15 MeV/c^2; decay length > 2*sigma from the IP.")
   .note(:ks_multi_candidate,
         "If multiple K_S0 candidates survive in one event, the one with the longest decay length " \
         "is retained (arbitration performed in ROOT after the BOSS ntuple is produced).")
   .note(:gamma_conversion_veto,
         "Reject gamma-conversion Bhabha background by requiring the opening angle between the two " \
         "pions from the K_S0 decay to be greater than 30 degrees (theta_pipi > 30 deg).")
   .note(:pwa_downstream,
         "PWA fit (helicity amplitude, TF-PWA framework) is performed downstream in ROOT; the BOSS " \
         "algorithm only produces the pre-PWA candidate NTuple after the 5C fit.")

alg.with_decay_card(decay_card_ks_k_pi).apply(sel)
alg.execute_on(scan_data_points + inc_mc_points + sig_mc)
