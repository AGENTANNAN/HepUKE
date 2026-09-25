### Dataset description ###
# 26-point J/psi scan between 3000.00 and 3119.88 MeV (total 440.7 pb^-1).
# The full set of scan points is not enumerated in the sample table used here;
# we list the ones that can be resolved and rely on the .note() field to record
# the intended full list. Any additional points would be added the same way.
scan_data_points = [
  DatasetManager.real_data.find("713_Rscan_3080"),
  DatasetManager.real_data.find("713_Rscan_3020"),
  DatasetManager.real_data.find("713_Rscan_3000"),
  DatasetManager.real_data.find("713_Rscan_2981"),
  DatasetManager.real_data.find("713_Rscan_2950")
]
inc_mc_points = scan_data_points.map { |d| DatasetManager.inclusive_mc.find("713_#{d.sample_name}") }

# ------------------------------------------------------------------
# Signal decay cards - use ConExc for the continuum + J/psi lineshape scan.
# Mode index 9 is e+e- -> K_S K+ pi- (ConExc mode table).
# The 'Particle vpho' line is intentionally omitted so that the DSL will
# inject the correct per-energy-point value automatically.
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

sig_mc_ks_k_pi = DatasetManager.create_exclusive_mc_for(scan_data_points) do |config|
  config.sample_name   = "ee_KsKpi_ConExc"
  config.events        = 200000
  config.decay_card    = decay_card_ks_k_pi
  config.cross_section = :default
end
sig_mc_ks_k_pi.each { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }

########################################################################
# Algorithm: e+e- -> K_S0 K+ pi-  (with subsequent K_S0 -> pi+ pi-)
########################################################################
alg = Algorithm.new("EEtoKsKPi")
alg.set_header(["EEtoKsKPiAlg/EEtoKsKPi.h"])
   .set_constant({"ECMS" => [:double, 3.080]})   # nominal; per-run value used at runtime

sel = Selection.new
sel.select_track {                       # For non-K_S tracks: |Vz| < 10 cm, |Vxy| < 1 cm; K_S tracks unconstrained
      cos_theta 0.93
      Vz        10.0
      Vr        1.0
      nChrp     ">=2"                    # K+ and pi+ from K_S0
      nChrn     ">=2"                    # pi- and pi- from K_S0
      nNet      "==0"
    }
   .select_photon {                      # Kept minimal; the analysis uses no radiative photons in the fit
      tdc_emc_start      0
      tdc_emc_end        14
      angle_to_track     10.0
      energyThreshold_b  0.025
      energyThreshold_e  0.050
    }
   .pid(method: :probability) {
      prob_cut 0.001
      identify :kaon, against: [:pion, :proton]
      identify :pion, against: [:kaon, :proton]
    }
   # Reconstruct K_S0 -> pi+ pi- via a secondary vertex fit (longest decay length pick made in ROOT)
   .secondary_vertex_fit([:pip, :pim]) {
      build_virtual_particle(:K_S0).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
    }
   # 4C kinematic fit under e+e- -> K_S0 K+ pi- hypothesis
   .kinematic_fit([:K_S0, :kp, :pim]) {
      nominal
      constrain_four_momentum
      chi2_cut 200
    }
   # 5C fit: 4C plus M(pi+ pi-) constrained to K_S0 nominal mass, for improved resolution
   .kinematic_fit([:K_S0, :kp, :pim]) {
      constrain_four_momentum
      invariant_mass_of(:K_S0).constrain_to_nominal_mass_of(:K_S0)
    }

alg.note(:dataset_scan_points,
         "Full analysis uses 26 e+e- energy points between 3000.00 and 3119.88 MeV. " \
         "The DSL enumerates the subset available in the current dataset table (713 Rscan_29xx/30xx). " \
         "Additional points must be added to scan_data_points when the table is completed.")
   .note(:ks_multi_candidate,
         "If multiple K_S0 candidates survive per event, the one with the longest decay length is retained.")
   .note(:gamma_conversion_veto,
         "Reject gamma-conversion Bhabha events by requiring the opening angle between the two " \
         "pions from the K_S0 decay to be greater than 30 degrees.")
   .note(:selection_ref,
         "Selections for pi+/-, K+/-, K_S0 candidates follow Ref. [43] (accompanying PRD version): " \
         "pi/K distinguished by combined MDC dE/dx and TOF likelihoods.")
   .note(:pwa_reference,
         "PWA is performed downstream in ROOT (TF-PWA framework); this BOSS spec only produces " \
         "the pre-PWA candidate NTuple after the kinematic fits.")

alg.with_decay_card(decay_card_ks_k_pi).apply(sel)
alg.execute_on(scan_data_points + inc_mc_points + sig_mc_ks_k_pi)
