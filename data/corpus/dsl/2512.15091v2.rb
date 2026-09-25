### Dataset description ###
# XYZ scan data at c.m. energies 4.16 - 4.34 GeV (total 10.9 fb^-1, 15 energy points).
datasets_xyz = [
  DatasetManager.real_data.find("703_4180"),
  DatasetManager.real_data.find("703_4190"),
  DatasetManager.real_data.find("703_4200"),
  DatasetManager.real_data.find("703_4210"),
  DatasetManager.real_data.find("703_4220"),
  DatasetManager.real_data.find("703_4230"),
  DatasetManager.real_data.find("703_4237"),
  DatasetManager.real_data.find("703_4246"),
  DatasetManager.real_data.find("703_4260"),
  DatasetManager.real_data.find("703_4270"),
  DatasetManager.real_data.find("703_4280"),
  DatasetManager.real_data.find("705_4290"),
  DatasetManager.real_data.find("705_4315"),
  DatasetManager.real_data.find("705_4340"),
]

# --- Decay cards ---
decay_card_KsKpi = <<~DECAYCARD
  Decay psi(4260)
  1.000 gamma X_1(3872)                    PHSP;
  Enddecay

  Decay X_1(3872)
  1.000 K_S0 K+ pi-                        PHSP;
  Enddecay

  Decay K_S0
  1.000 pi+ pi-                            PHSP;
  Enddecay

  End
DECAYCARD

decay_card_KstarK = <<~DECAYCARD
  Decay psi(4260)
  1.000 gamma X_1(3872)                    PHSP;
  Enddecay

  Decay X_1(3872)
  0.500 K*+ K-                             VSS;
  0.500 K*0 anti-K0                        VSS;
  Enddecay

  Decay K*+
  1.000 K0 pi+                             VSS;
  Enddecay

  Decay K*0
  1.000 K+ pi-                             VSS;
  Enddecay

  Decay K0
  1.000 K_S0                               PHSP;
  Enddecay

  Decay anti-K0
  1.000 K_S0                               PHSP;
  Enddecay

  Decay K_S0
  1.000 pi+ pi-                            PHSP;
  Enddecay

  End
DECAYCARD

# --- Exclusive MC samples for each energy point ---
exMC_KsKpi_list = DatasetManager.create_exclusive_mc_for(datasets_xyz) do |c|
  c.sample_name    = "X3872_KsKpi"
  c.events         = 100_000
  c.decay_card     = decay_card_KsKpi
  c.cross_section  = :default
end
exMC_KsKpi_list.each { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }

exMC_KstarK_list = DatasetManager.create_exclusive_mc_for(datasets_xyz) do |c|
  c.sample_name    = "X3872_KstarK"
  c.events         = 100_000
  c.decay_card     = decay_card_KstarK
  c.cross_section  = :default
end
exMC_KstarK_list.each { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }

### Event selection ###
alg_name = "X3872ToKsKPi"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({ "ECMS" => [:double, 4.226] })
   .set_alias({ "std::vector<double>" => "Vdouble" })

sel = Selection.new
sel.select_track {
      cos_theta 0.93
      Vz 10.0
      Vr 1.0
      nChrp "==2"
      nChrn "==2"
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
       identify :kaon, against: [:pion, :proton]
       nkp ">=1 || nkm >=1"
    }
   .secondary_vertex_fit([:pip, :pim]) {
      build_virtual_particle(:K_S0).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
   }
   .kinematic_fit([:gamma, :K_S0, :kp, :pim]) {
      nominal
      constrain_four_momentum
      chi2_cut 200
   }
   # Alternative fits for photon-veto/photon-add veto stored to NTuple (Rule T2)
   .kinematic_fit([:K_S0, :kp, :pim]) {
      constrain_four_momentum
   }
   .kinematic_fit([:gamma, :gamma, :K_S0, :kp, :pim]) {
      constrain_four_momentum
   }

alg.note(:ks0_mass_window,
         "|M(pi+ pi-) - m_KS0| < 22 MeV/c^2, decay length > 2 sigma from IP, secondary vertex fit " \
         "smallest chi2 combination kept.")
   .note(:ks0_sideband,
         "KS0 sideband 0.432 < M(pi+pi-) < 0.454 and 0.542 < M(pi+pi-) < 0.564 GeV/c^2 used to " \
         "estimate non-KS0 background (ROOT stage).")
   .note(:vertex_fit_IP,
         "Vertex fit to K, remaining pi and KS0 to common IP vertex.")
   .note(:pi0_veto,
         "M_recoil^2(KS0 K+- pi-+) < 0.012 (GeV/c^2)^2 to suppress pi0 backgrounds (ROOT stage).")
   .note(:kstar_selection,
         "K*(892) candidates: 0.748 < M(K+- pi-+) < 1.036 GeV/c^2 for neutral K* or " \
         "0.743 < M(KS0 pi+-) < 1.05 GeV/c^2 for charged K* (ROOT stage).")
   .note(:helix_correction,
         "helix-parameter correction applied to charged tracks before 4C kinematic fit.")

alg.with_decay_card(decay_card_KsKpi).apply(sel)
alg.execute_on(datasets_xyz + exMC_KsKpi_list + exMC_KstarK_list)
