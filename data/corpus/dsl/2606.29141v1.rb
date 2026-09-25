### dataset description ###
# 15 c.m. energy points from sqrt(s) = 4.4156 to 4.9509 GeV; total 8.5 fb^-1.
# Representative scan points listed here; expand data_scan_points for production.
data_scan_points = [
  DatasetManager.real_data.find("703_4260"),  # 4.260 GeV region (proxy sample name)
  DatasetManager.real_data.find("703_4600"),  # 4.600 GeV region
  DatasetManager.real_data.find("703_4680"),  # 4.680 GeV region
]

incMC_samples = data_scan_points.map do |ds|
  DatasetManager.inclusive_mc.find("#{ds.boss.gsub('.', '')[0,3]}_#{ds.energy.to_i}")
end

# -------------- Decay cards for the four (Ds+,Ds-) sub-modes --------------

# Mode I: Ds+ -> K+ K- pi+ (D_DALITZ), Ds- -> K+ K- pi- (D_DALITZ)
decay_card_KKpi_KKpi = <<~DECAYCARD
    Decay psi(4260)
    1.000  pi+  pi-  D_s+  D_s-           PHSP;
    Enddecay

    Decay D_s+
    1.000  K+  K-  pi+                    D_DALITZ;
    Enddecay

    Decay D_s-
    1.000  K+  K-  pi-                    D_DALITZ;
    Enddecay

    End
DECAYCARD

# Mode II: Ds+ -> K+ K- pi+ (D_DALITZ), Ds- -> K_S0 K- (PHSP)
decay_card_KKpi_KsK = <<~DECAYCARD
    Decay psi(4260)
    1.000  pi+  pi-  D_s+  D_s-           PHSP;
    Enddecay

    Decay D_s+
    1.000  K+  K-  pi+                    D_DALITZ;
    Enddecay

    Decay D_s-
    1.000  K_S0  K-                       PHSP;
    Enddecay

    Decay K_S0
    1.000  pi+  pi-                       PHSP;
    Enddecay

    End
DECAYCARD

# Mode III: Ds+ -> K_S0 K+ (PHSP), Ds- -> K+ K- pi- (D_DALITZ)
decay_card_KsK_KKpi = <<~DECAYCARD
    Decay psi(4260)
    1.000  pi+  pi-  D_s+  D_s-           PHSP;
    Enddecay

    Decay D_s+
    1.000  K_S0  K+                       PHSP;
    Enddecay

    Decay D_s-
    1.000  K+  K-  pi-                    D_DALITZ;
    Enddecay

    Decay K_S0
    1.000  pi+  pi-                       PHSP;
    Enddecay

    End
DECAYCARD

# Mode IV: Ds+ -> K_S0 K+ (PHSP), Ds- -> K_S0 K- (PHSP)
decay_card_KsK_KsK = <<~DECAYCARD
    Alias K_S0_1 K_S0
    Alias K_S0_2 K_S0

    Decay psi(4260)
    1.000  pi+  pi-  D_s+  D_s-           PHSP;
    Enddecay

    Decay D_s+
    1.000  K_S0_1  K+                     PHSP;
    Enddecay

    Decay D_s-
    1.000  K_S0_2  K-                     PHSP;
    Enddecay

    Decay K_S0_1
    1.000  pi+  pi-                       PHSP;
    Enddecay

    Decay K_S0_2
    1.000  pi+  pi-                       PHSP;
    Enddecay

    End
DECAYCARD

# -------------- Exclusive MC per sub-mode across scan points ---------------
exMC_KKpi_KKpi = DatasetManager.create_exclusive_mc_for(data_scan_points) do |c|
  c.sample_name = "exmc_pipiDsDs_KKpi_KKpi"
  c.events = 100_000
  c.decay_card = decay_card_KKpi_KKpi
  c.cross_section = :default
end

exMC_KKpi_KsK = DatasetManager.create_exclusive_mc_for(data_scan_points) do |c|
  c.sample_name = "exmc_pipiDsDs_KKpi_KsK"
  c.events = 100_000
  c.decay_card = decay_card_KKpi_KsK
  c.cross_section = :default
end

exMC_KsK_KKpi = DatasetManager.create_exclusive_mc_for(data_scan_points) do |c|
  c.sample_name = "exmc_pipiDsDs_KsK_KKpi"
  c.events = 100_000
  c.decay_card = decay_card_KsK_KKpi
  c.cross_section = :default
end

exMC_KsK_KsK = DatasetManager.create_exclusive_mc_for(data_scan_points) do |c|
  c.sample_name = "exmc_pipiDsDs_KsK_KsK"
  c.events = 100_000
  c.decay_card = decay_card_KsK_KsK
  c.cross_section = :default
end

# ============ Algorithm I: (Ds+ -> K+K-pi+) + (Ds- -> K+K-pi-) ============
# Final state: pi+ pi- K+ K- K+ K- pi+ pi- => 3(pi+) 3(pi-) 2(K+) 2(K-)
alg_KKpi_KKpi = Algorithm.new("PipiDsDsKKpiKKpi")
alg_KKpi_KKpi.set_header(["PipiDsDsKKpiKKpiAlg/PipiDsDsKKpiKKpi.h"])
             .set_alias({"std::vector<double>" => "Vdouble"})
sel_KKpi_KKpi = Selection.new
sel_KKpi_KKpi.select_track {
                  cos_theta 0.93
                  Vz 10.0
                  Vr 1.0
                  nChrp ">=3"
                  nChrn ">=3"
                  nNet "==0"
                }
               .pid(method: :probability) {
                  prob_cut 0.001
                  identify :pion, against: [:kaon]
                  identify :kaon, against: [:pion]
                  npip ">=3"
                  npim ">=3"
                  nkp  ">=2"
                  nkm  ">=2"
               }
               .kinematic_fit([:pip, :pip, :pip, :pim, :pim, :pim, :kp, :kp, :km, :km]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 200
               }
alg_KKpi_KKpi
  .note(:Ds_mode_pairing, "Ds+ -> K+K-pi+ (D_DALITZ) paired with Ds- -> K+K-pi- (D_DALITZ); optimal (Ds+,Ds-) combination chosen by minimising |avg(M(Ds+)+M(Ds-))/2 - M_PDG(Ds)|")
  .note(:phi_veto_KKpi, "for D_s^+/- -> K+ K- pi+/-, phi (K+K-) mass window: M(K+K-) < 1.2 GeV/c^2; K*bar(892)^0 (K-pi+) window: M(K-pi+) < 1.0 GeV/c^2 (also on charge-conjugate)")
  .note(:Ds_signal_mass_window, "|M(Ds-) - M_PDG(Ds-)| < 15 MeV/c^2 imposed on the recoiling Ds when fitting the recoil mass RM(pi+ pi- Ds+/-)")
  .with_decay_card(decay_card_KKpi_KKpi)
  .apply(sel_KKpi_KKpi)
alg_KKpi_KKpi.execute_on(data_scan_points + incMC_samples + exMC_KKpi_KKpi)

# ============ Algorithm II: (Ds+ -> K+K-pi+) + (Ds- -> K_S0 K-) ============
# Ds-'s K_S0 -> pi+ pi- is reconstructed via a secondary vertex fit.
# Final state: pi+ pi- (K+ K- pi+) (K_S0 K-)
# After Ks reconstruction: pi+ pi- K+ K- pi+ K_S0 K-
alg_KKpi_KsK = Algorithm.new("PipiDsDsKKpiKsK")
alg_KKpi_KsK.set_header(["PipiDsDsKKpiKsKAlg/PipiDsDsKKpiKsK.h"])
            .set_alias({"std::vector<double>" => "Vdouble"})
sel_KKpi_KsK = Selection.new
sel_KKpi_KsK.select_track {
                  cos_theta 0.93
                  Vz 20.0        # loosened to accommodate K_S0 daughter tracks
                  Vr 10.0
                  nChrp ">=3"
                  nChrn ">=3"
                  nNet "==0"
                }
              .pid(method: :probability) {
                  prob_cut 0.001
                  identify :pion, against: [:kaon]
                  identify :kaon, against: [:pion]
                  npip ">=3"
                  npim ">=3"
                  nkp  ">=1"
                  nkm  ">=2"
               }
              .secondary_vertex_fit([:pip, :pim]) {
                  build_virtual_particle(:K_S0).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
               }
              .kinematic_fit([:pip, :pip, :pim, :pim, :kp, :km, :km, :K_S0]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 200
               }
alg_KKpi_KsK
  .note(:KS0_selection, "K_S0 daughter tracks: |Vz|<20 cm, |cos theta|<0.93, no PID; vertex-fitted invariant mass in (0.487, 0.511) GeV/c^2 and decay length > 2*resolution")
  .note(:phi_veto_KKpi, "for D_s+ -> K+ K- pi+, phi and K*(892) mass windows applied as in mode I")
  .note(:Ds_pairing, "optimal (Ds+,Ds-) combination selected by minimising Delta_hat_M = |0.5(M(Ds+)+M(Ds-)) - M_PDG(Ds)|")
  .with_decay_card(decay_card_KKpi_KsK)
  .apply(sel_KKpi_KsK)
alg_KKpi_KsK.execute_on(data_scan_points + incMC_samples + exMC_KKpi_KsK)

# ============ Algorithm III: (Ds+ -> K_S0 K+) + (Ds- -> K+K-pi-) ============
alg_KsK_KKpi = Algorithm.new("PipiDsDsKsKKKpi")
alg_KsK_KKpi.set_header(["PipiDsDsKsKKKpiAlg/PipiDsDsKsKKKpi.h"])
            .set_alias({"std::vector<double>" => "Vdouble"})
sel_KsK_KKpi = Selection.new
sel_KsK_KKpi.select_track {
                  cos_theta 0.93
                  Vz 20.0
                  Vr 10.0
                  nChrp ">=3"
                  nChrn ">=3"
                  nNet "==0"
                }
              .pid(method: :probability) {
                  prob_cut 0.001
                  identify :pion, against: [:kaon]
                  identify :kaon, against: [:pion]
                  npip ">=3"
                  npim ">=3"
                  nkp  ">=2"
                  nkm  ">=1"
               }
              .secondary_vertex_fit([:pip, :pim]) {
                  build_virtual_particle(:K_S0).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
               }
              .kinematic_fit([:pip, :pip, :pim, :pim, :kp, :kp, :km, :K_S0]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 200
               }
alg_KsK_KKpi
  .note(:phi_veto_KKpi, "for D_s- -> K+ K- pi-, phi and K*(892) mass windows applied as in mode I")
  .with_decay_card(decay_card_KsK_KKpi)
  .apply(sel_KsK_KKpi)
alg_KsK_KKpi.execute_on(data_scan_points + incMC_samples + exMC_KsK_KKpi)

# ============ Algorithm IV: (Ds+ -> K_S0 K+) + (Ds- -> K_S0 K-) ============
alg_KsK_KsK = Algorithm.new("PipiDsDsKsKKsK")
alg_KsK_KsK.set_header(["PipiDsDsKsKKsKAlg/PipiDsDsKsKKsK.h"])
           .set_alias({"std::vector<double>" => "Vdouble"})
sel_KsK_KsK = Selection.new
sel_KsK_KsK.select_track {
                  cos_theta 0.93
                  Vz 20.0
                  Vr 10.0
                  nChrp ">=3"
                  nChrn ">=3"
                  nNet "==0"
                }
             .pid(method: :probability) {
                  prob_cut 0.001
                  identify :pion, against: [:kaon]
                  identify :kaon, against: [:pion]
                  npip ">=3"
                  npim ">=3"
                  nkp  ">=1"
                  nkm  ">=1"
              }
             .secondary_vertex_fit([:pip, :pim]) {
                  build_virtual_particle(:K_S0).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
              }
             .secondary_vertex_fit([:pip, :pim]) {
                  build_virtual_particle(:K_S0).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
              }
             .kinematic_fit([:pip, :pim, :kp, :km, :K_S0, :K_S0]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 200
              }
alg_KsK_KsK
  .note(:KS0_two_candidates, "two independent K_S0 -> pi+ pi- candidates required; both reconstructed via secondary vertex fits and both used in the kinematic fit")
  .with_decay_card(decay_card_KsK_KsK)
  .apply(sel_KsK_KsK)
alg_KsK_KsK.execute_on(data_scan_points + incMC_samples + exMC_KsK_KsK)
