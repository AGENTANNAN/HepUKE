### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")     # psi(3686) real data (2.712 x 10^9 events)
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # psi(3686) inclusive MC

# --- Decay card: psi(3686) -> gamma chi_cJ (using chi_c1 as generic), chi_c1 -> p pbar KS0 KS0 ---
signal_decay_card = <<~DECAYCARD
  Alias K_S0_a K_S0
  Alias K_S0_b K_S0

  Decay psi(3686)
  1.000 gamma chi_c1                       PHSP;
  Enddecay

  Decay chi_c1
  1.000 p+ anti-p- K_S0_a K_S0_b           PHSP;
  Enddecay

  Decay K_S0_a
  1.000 pi+ pi-                            PHSP;
  Enddecay

  Decay K_S0_b
  1.000 pi+ pi-                            PHSP;
  Enddecay

  End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |c|
  c.sample_name    = "chi_cJ_ppKsKs"
  c.related_dataset = psip_data
  c.events         = 1_000_000
  c.decay_card     = signal_decay_card
  c.cross_section  = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection ###
alg_name = "ChicJToPpKsKs"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({ "ECMS" => [:double, 3.686] })
   .set_alias({ "std::vector<double>" => "Vdouble" })

sel = Selection.new
sel.select_track {
      cos_theta 0.93
      Vz 10.0
      Vr 1.0
      nChrp ">=3"
      nChrn ">=3"
      nNet  "==0"
    }
   .select_photon {
      tdc_emc_start 0
      tdc_emc_end 14
      angle_to_track 10.0
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      nGam ">=1"
    }
   .pid(method: :probability) {
      prob_cut 0.001
      identify :proton, against: [:kaon, :pion]
      nprp ">=1"
      nprm ">=1"
    }
   # First KS0 -> pi+ pi-
   .secondary_vertex_fit([:pip, :pim]) {
      build_virtual_particle(:K_S0).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
   }
   # Second KS0 -> pi+ pi-
   .secondary_vertex_fit([:pip, :pim]) {
      build_virtual_particle(:K_S0).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
   }
   # 4C kinematic fit e+e- -> gamma p pbar KS0 KS0
   .kinematic_fit([:gamma, :prp, :prm, :K_S0, :K_S0]) {
      nominal
      constrain_four_momentum
      chi2_cut 200
   }

alg.note(:ks0_mass_window,
         "|M(pi+ pi-) - m_KS0| < 0.012 GeV/c^2 (3 sigma), vertex + secondary vertex fit, " \
         "decay length > 2 sigma away from IP.")
   .note(:lambda_veto,
         "|M(p pi-) - m_Lambda| > 0.012 GeV/c^2 and |M(pbar pi+) - m_Lambda_bar| > 0.012 GeV/c^2 " \
         "to suppress Lambda/anti-Lambda backgrounds.")
   .note(:chi2_4c_cut,
         "Optimized chi^2_4C < 50 in the paper; loose 200 used in BOSS per Rule T3, tight cut in ROOT.")
   .note(:ks0_sideband,
         "1D KS0 sideband (0.454,0.478) U (0.518,0.542) GeV/c^2 used to build 2D KS0 KS0 sideband " \
         "regions SB1 and SB2 for peaking-background subtraction (ROOT stage).")
   .note(:helix_correction,
         "Helix-parameter correction applied to charged tracks before 4C kinematic fit " \
         "(systematic-uncertainty study).")

alg.with_decay_card(signal_decay_card).apply(sel)
alg.execute_on([psip_data, psip_incMC, exMC_signal])
