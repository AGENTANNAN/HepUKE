# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# XYZ-scan real-data c.m. energy points from 4.180 to 4.340 GeV (no inclusive MC)
data_points = %w[
  703_4180 703_4190 703_4200 703_4210 703_4220 703_4230 703_4237
  703_4246 703_4260 703_4270 703_4280 705_4290 703_4310 705_4315 705_4340
].map { |name| DatasetManager.real_data.find(name) }

# ---------------- Decay cards (EvtGen) ----------------
# Mode I : psi(4260) -> gamma X(3872), X(3872) -> K_S0 K+ pi-, K_S0 -> pi+ pi-
decay_card_mode1 = <<~DECAYCARD
  Decay psi(4260)
  1.000 gamma X(3872) PHSP;
  Enddecay

  Decay X(3872)
  1.000 K_S0 K+ pi- PHSP;
  Enddecay

  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# Mode II: psi(4260) -> gamma X(3872), X(3872) -> K*(892) K
#          (K*+ K- / K*0 anti-K0), K* -> K pi, K0 -> K_S0 -> pi+ pi-
decay_card_mode2 = <<~DECAYCARD
  Decay psi(4260)
  1.000 gamma X(3872) PHSP;
  Enddecay

  Decay X(3872)
  0.500 K*+ K- PHSP;
  0.500 K*0 anti-K0 PHSP;
  Enddecay

  Decay K*+
  1.000 K0 pi+ VSS;
  Enddecay

  Decay K*0
  1.000 K+ pi- VSS;
  Enddecay

  Decay K0
  1.000 K_S0 PHSP;
  Enddecay

  Decay anti-K0
  1.000 K_S0 PHSP;
  Enddecay

  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# ---------------- Exclusive MC: 100k events per energy point, per mode ----------------
exMC_mode1 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_X3872_KsKpi"
  config.events        = 100_000
  config.decay_card    = decay_card_mode1
  config.cross_section = :default
end

exMC_mode2 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_X3872_KstarK"
  config.events        = 100_000
  config.decay_card    = decay_card_mode2
  config.cross_section = :default
end

### Event selection (BOSS) ###

# ================= Mode I: X(3872) -> K_S0 K+ pi- =================
alg_name1 = "X3872KsKpi"
alg1 = Algorithm.new(alg_name1)
alg1.set_header(["#{alg_name1}Alg/#{alg_name1}.h"])
    .set_constant({"ECMS" => [:double, 4.226]})
    .note(:helix_correction,
          "helix-parameter correction applied to all charged tracks before the 4C kinematic fit; the efficiency difference is estimated by re-running the BOSS selection with and without the correction")
    .note(:kaon_multiplicity,
          "at least one identified charged kaon (K+/-) required after PID (the K_S0 daughters are pions)")
    .note(:k_s0_selection,
          "K_S0 candidates built from pi+ pi- keeping the smallest |M(pi+pi-) - m_K_S0|; required |M(pi+pi-) - m_K_S0| < 22 MeV/c^2 and a decay length > 2 sigma from the IP")
    .note(:background_veto,
          "pi0 background suppressed by M_recoil^2(K_S0 K+/- pi-/+ ) < 0.012 (GeV/c^2)^2; K_S0 sidebands 0.432-0.454 and 0.542-0.564 GeV/c^2 used to estimate the non-K_S0 background")

event_selection_mode1 = Selection.new
event_selection_mode1
  .select_track {
      cos_theta 0.93     # |cos(theta)| < 0.93
      Vz        10.0     # |Vz| < 10 cm
      Vr        1.0      # Vr < 1 cm
      nChrp     "==2"    # exactly two positive tracks
      nChrn     "==2"    # exactly two negative tracks
      nNet      "==0"    # net charge zero
  }
  .select_photon {
      energyThreshold_b 0.025   # barrel energy > 25 MeV
      energyThreshold_e 0.050   # endcap energy > 50 MeV
      tdc_emc_start     0       # EMC TDC in [0, 14]
      tdc_emc_end       14
      angle_to_track    10.0    # angle to any track > 10 degrees
      nGam              ">=1"   # at least one photon
  }
  .pid(method: :probability) {
      prob_cut 0.001
      identify :kaon, against: [:pion, :proton]
  }
  .remove([:kp <= :chrgp, :km <= :chrgn])     # keep kaons out of the pion lists
  .assign({:chrgp => :pip, :chrgn => :pim})   # remaining tracks are pions
  .secondary_vertex_fit([:pip, :pim]) {       # K_S0 from pi+ pi-
      build_virtual_particle(:K_S0).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:gamma, :K_S0, :kp, :pim]) {   # nominal 4C fit of gamma K_S0 K+ pi-
      nominal
      constrain_four_momentum
      chi2_cut 200
      vertex_fit([2, 3])                         # K_S0(daughters), kaon and pion to a common IP
  }
  .kinematic_fit([:K_S0, :kp, :pim]) {           # 3C fit without the photon (saved)
      constrain_three_momentum
  }
  .kinematic_fit([:gamma, :gamma, :K_S0, :kp, :pim]) {  # 5C fit with an extra photon (saved)
      constrain_four_momentum
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  }

alg1.with_decay_card(decay_card_mode1).apply(event_selection_mode1)

# ================= Mode II: X(3872) -> K*(892) K =================
alg_name2 = "X3872KstarK"
alg2 = Algorithm.new(alg_name2)
alg2.set_header(["#{alg_name2}Alg/#{alg_name2}.h"])
    .set_constant({"ECMS" => [:double, 4.226]})
    .note(:helix_correction,
          "helix-parameter correction applied to all charged tracks before the 4C kinematic fit")
    .note(:kaon_multiplicity,
          "at least one identified charged kaon (K+/-) required after PID")
    .note(:k_s0_selection,
          "K_S0 candidates built from pi+ pi- keeping the smallest |M(pi+pi-) - m_K_S0|; required |M(pi+pi-) - m_K_S0| < 22 MeV/c^2 and a decay length > 2 sigma from the IP")
    .note(:kstar_candidates,
          "K*(892) candidates required to satisfy 0.748 < M(K+/- pi-/+ ) < 1.036 GeV/c^2 (neutral K*) or 0.743 < M(K_S0 pi+/- ) < 1.05 GeV/c^2 (charged K*)")
    .note(:background_veto,
          "pi0 background suppressed by M_recoil^2(K_S0 K+/- pi-/+ ) < 0.012 (GeV/c^2)^2")

event_selection_mode2 = Selection.new
event_selection_mode2
  .select_track {
      cos_theta 0.93
      Vz        10.0
      Vr        1.0
      nChrp     "==2"
      nChrn     "==2"
      nNet      "==0"
  }
  .select_photon {
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      tdc_emc_start     0
      tdc_emc_end       14
      angle_to_track    10.0
      nGam              ">=1"
  }
  .pid(method: :probability) {
      prob_cut 0.001
      identify :kaon, against: [:pion, :proton]
  }
  .remove([:kp <= :chrgp, :km <= :chrgn])
  .assign({:chrgp => :pip, :chrgn => :pim})
  .secondary_vertex_fit([:pip, :pim]) {
      build_virtual_particle(:K_S0).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:gamma, :K_S0, :kp, :pim]) {
      nominal
      constrain_four_momentum
      chi2_cut 200
      vertex_fit([2, 3])
  }
  .kinematic_fit([:K_S0, :kp, :pim]) {
      constrain_three_momentum
  }
  .kinematic_fit([:gamma, :gamma, :K_S0, :kp, :pim]) {
      constrain_four_momentum
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  }

alg2.with_decay_card(decay_card_mode2).apply(event_selection_mode2)

# ---------------- Execute ----------------
alg1.execute_on(data_points + exMC_mode1)
alg2.execute_on(data_points + exMC_mode2)