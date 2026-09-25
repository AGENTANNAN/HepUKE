# =====================================================================
# e+e- -> pi+ pi- Ds+ Ds-   (15-point c.m. energy scan, 4.4156-4.9509 GeV)
# BOSS / event-selection (HepScript) specification
# =====================================================================

### Dataset description ###
# 15 c.m. energy scan points, total 8.5 fb^-1
scan_points = [
  DatasetManager.real_data.find("703_4420"),   # ~4.4156 GeV
  DatasetManager.real_data.find("703_4470"),
  DatasetManager.real_data.find("703_4530"),
  DatasetManager.real_data.find("703_4575"),
  DatasetManager.real_data.find("703_4600"),
  DatasetManager.real_data.find("706_4620"),
  DatasetManager.real_data.find("706_4640"),
  DatasetManager.real_data.find("706_4660"),
  DatasetManager.real_data.find("706_4680"),
  DatasetManager.real_data.find("706_4700"),
  DatasetManager.real_data.find("707_4740"),
  DatasetManager.real_data.find("707_4750"),
  DatasetManager.real_data.find("707_4780"),
  DatasetManager.real_data.find("707_4840"),
  DatasetManager.real_data.find("707_4946")    # ~4.9509 GeV
]

# corresponding inclusive MC samples (background estimation)
scan_incMC = [
  DatasetManager.inclusive_mc.find("703_4420"),
  DatasetManager.inclusive_mc.find("703_4470"),
  DatasetManager.inclusive_mc.find("703_4530"),
  DatasetManager.inclusive_mc.find("703_4575"),
  DatasetManager.inclusive_mc.find("703_4600"),
  DatasetManager.inclusive_mc.find("706_4620"),
  DatasetManager.inclusive_mc.find("706_4640"),
  DatasetManager.inclusive_mc.find("706_4660"),
  DatasetManager.inclusive_mc.find("706_4680"),
  DatasetManager.inclusive_mc.find("706_4700"),
  DatasetManager.inclusive_mc.find("707_4740"),
  DatasetManager.inclusive_mc.find("707_4750"),
  DatasetManager.inclusive_mc.find("707_4780"),
  DatasetManager.inclusive_mc.find("707_4840"),
  DatasetManager.inclusive_mc.find("707_4946")
]

### Decay cards (EvtGen syntax) — one per Ds-pair sub-mode ###

# Mode I : (K+K-pi+)(K+K-pi-)
decay_card_mode1 = <<~DECAYCARD
  Decay psi(4260)
  1.0 pi+ pi- Ds+ Ds- PHSP;
  Enddecay

  Decay Ds+
  1.0 K+ K- pi+ DS_DALITZ;
  Enddecay

  Decay Ds-
  1.0 K+ K- pi- DS_DALITZ;
  Enddecay
  End
DECAYCARD

# Mode II : (K+K-pi+)(K_S0 K-)
decay_card_mode2 = <<~DECAYCARD
  Decay psi(4260)
  1.0 pi+ pi- Ds+ Ds- PHSP;
  Enddecay

  Decay Ds+
  1.0 K+ K- pi+ DS_DALITZ;
  Enddecay

  Decay Ds-
  1.0 K- K_S0 PHSP;
  Enddecay

  Decay K_S0
  1.0 pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

# Mode III : (K_S0 K+)(K+K-pi-)
decay_card_mode3 = <<~DECAYCARD
  Decay psi(4260)
  1.0 pi+ pi- Ds+ Ds- PHSP;
  Enddecay

  Decay Ds+
  1.0 K+ K_S0 PHSP;
  Enddecay

  Decay Ds-
  1.0 K+ K- pi- DS_DALITZ;
  Enddecay

  Decay K_S0
  1.0 pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

# Mode IV : (K_S0 K+)(K_S0 K-)
decay_card_mode4 = <<~DECAYCARD
  Decay psi(4260)
  1.0 pi+ pi- Ds+ Ds- PHSP;
  Enddecay

  Decay Ds+
  1.0 K+ K_S0 PHSP;
  Enddecay

  Decay Ds-
  1.0 K- K_S0 PHSP;
  Enddecay

  Decay K_S0
  1.0 pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

### Exclusive MC — 100k events for each sub-mode at every energy point ###
exMC_mode1 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "exmc_pipidsdsbar_mode1"
  config.events        = 100_000
  config.decay_card    = decay_card_mode1
  config.cross_section = :default
end

exMC_mode2 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "exmc_pipidsdsbar_mode2"
  config.events        = 100_000
  config.decay_card    = decay_card_mode2
  config.cross_section = :default
end

exMC_mode3 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "exmc_pipidsdsbar_mode3"
  config.events        = 100_000
  config.decay_card    = decay_card_mode3
  config.cross_section = :default
end

exMC_mode4 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "exmc_pipidsdsbar_mode4"
  config.events        = 100_000
  config.decay_card    = decay_card_mode4
  config.cross_section = :default
end

# =====================================================================
### Event selection (BOSS) ###
# =====================================================================

# ---------------------------------------------------------------------
# Mode I : e+e- -> pi+pi- (Ds+ -> K+K-pi+) (Ds- -> K+K-pi-)
# ---------------------------------------------------------------------
alg_name_1 = "DsDsbarMode1"
alg1 = Algorithm.new(alg_name_1)
alg1.set_header(["#{alg_name_1}Alg/#{alg_name_1}.h"])
    .set_constant({"ECMS" => [:double, 4.6]})   # representative c.m. energy; per-point value taken at execute_on
    .note(:energy_point_handling,
          "the analysis runs over 15 c.m. energy points (4.4156-4.9509 GeV); the ECMS constant is
           resolved per data / MC sample from the dataset c.m. energy when the job is generated")
    .note(:ds_pairing_criterion,
          "the (Ds+, Ds-) candidate pairing is chosen by minimising |0.5*(M(Ds+)+M(Ds-)) - M_PDG(Ds)|")
    .note(:recoil_mass_window,
          "|M(Ds) - M_PDG(Ds)| < 15 MeV/c^2 applied on the recoiling Ds in the RM(pi+pi-Ds) fit")

sel1 = Selection.new
sel1.select_track {
        cos_theta 0.93
        Vz        10.0
        Vr        1.0
        nChrp     ">=3"
        nChrn     ">=3"
        nNet      "==0"
    }
    .pid(method: :probability) {
        prob_cut 0.001
        identify :pion, against: [:kaon]   # pi+ and pi-
        identify :kaon, against: [:pion]   # K+ and K-
        npip ">=3"
        npim ">=3"
        nkp  ">=2"
        nkm  ">=2"
    }
    # nominal 4C fit on the full final state (pi+ pi- K+ K- pi+ K+ K- pi-)
    .kinematic_fit([:pip, :pip, :pim, :pim, :kp, :kp, :km, :km]) {
        nominal
        constrain_four_momentum
        invariant_mass_of(:kp, :km).larger_than(1.2)    # phi   veto  M(K+K-)  < 1.2 GeV/c^2
        invariant_mass_of(:km, :pip).larger_than(1.0)   # K*(892) veto M(K-pi+) < 1.0 GeV/c^2
        chi2_cut 200
    }

alg1.with_decay_card(decay_card_mode1).apply(sel1)

# ---------------------------------------------------------------------
# Mode II : e+e- -> pi+pi- (Ds+ -> K+K-pi+) (Ds- -> K_S0 K-)
# ---------------------------------------------------------------------
alg_name_2 = "DsDsbarMode2"
alg2 = Algorithm.new(alg_name_2)
alg2.set_header(["#{alg_name_2}Alg/#{alg_name_2}.h"])
    .set_constant({"ECMS" => [:double, 4.6]})
    .note(:energy_point_handling,
          "15 c.m. energy points (4.4156-4.9509 GeV); ECMS resolved per dataset at execute_on")
    .note(:ks_mass_window,
          "K_S0 candidates restricted to 0.487 < M(pi+pi-) < 0.511 GeV/c^2")
    .note(:ks_decay_length,
          "K_S0 decay length required to be greater than twice its resolution")
    .note(:ds_pairing_criterion,
          "the (Ds+, Ds-) candidate pairing is chosen by minimising |0.5*(M(Ds+)+M(Ds-)) - M_PDG(Ds)|")
    .note(:recoil_mass_window,
          "|M(Ds) - M_PDG(Ds)| < 15 MeV/c^2 applied on the recoiling Ds in the RM(pi+pi-Ds) fit")

sel2 = Selection.new
sel2.select_track {
        cos_theta 0.93
        Vz        20.0     # loosened for the K_S0 modes
        Vr        10.0     # loosened for the K_S0 modes
        nChrp     ">=3"
        nChrn     ">=3"
        nNet      "==0"
    }
    .pid(method: :probability) {
        prob_cut 0.001
        identify :pion, against: [:kaon]
        identify :kaon, against: [:pion]
        nkp ">=1"
        nkm ">=2"
    }
    # build K_S0 from pi+pi- via a secondary vertex fit; remove used daughters
    .secondary_vertex_fit([:pip, :pim]) {
        build_virtual_particle(:K_S0).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    # nominal 4C fit on the full final state (pi+ pi- K+ K- pi+ K- K_S0)
    .kinematic_fit([:pip, :pip, :pim, :kp, :km, :km, :K_S0]) {
        nominal
        constrain_four_momentum
        invariant_mass_of(:kp, :km).larger_than(1.2)    # phi   veto (Ds+ -> K+K-pi+)
        invariant_mass_of(:km, :pip).larger_than(1.0)   # K*(892) veto
        chi2_cut 200
    }

alg2.with_decay_card(decay_card_mode2).apply(sel2)

# ---------------------------------------------------------------------
# Mode III : e+e- -> pi+pi- (Ds+ -> K_S0 K+) (Ds- -> K+K-pi-)
# ---------------------------------------------------------------------
alg_name_3 = "DsDsbarMode3"
alg3 = Algorithm.new(alg_name_3)
alg3.set_header(["#{alg_name_3}Alg/#{alg_name_3}.h"])
    .set_constant({"ECMS" => [:double, 4.6]})
    .note(:energy_point_handling,
          "15 c.m. energy points (4.4156-4.9509 GeV); ECMS resolved per dataset at execute_on")
    .note(:ks_mass_window,
          "K_S0 candidates restricted to 0.487 < M(pi+pi-) < 0.511 GeV/c^2")
    .note(:ks_decay_length,
          "K_S0 decay length required to be greater than twice its resolution")
    .note(:ds_pairing_criterion,
          "the (Ds+, Ds-) candidate pairing is chosen by minimising |0.5*(M(Ds+)+M(Ds-)) - M_PDG(Ds)|")
    .note(:recoil_mass_window,
          "|M(Ds) - M_PDG(Ds)| < 15 MeV/c^2 applied on the recoiling Ds in the RM(pi+pi-Ds) fit")

sel3 = Selection.new
sel3.select_track {
        cos_theta 0.93
        Vz        20.0
        Vr        10.0
        nChrp     ">=3"
        nChrn     ">=3"
        nNet      "==0"
    }
    .pid(method: :probability) {
        prob_cut 0.001
        identify :pion, against: [:kaon]
        identify :kaon, against: [:pion]
        nkp ">=2"
        nkm ">=1"
    }
    .secondary_vertex_fit([:pip, :pim]) {
        build_virtual_particle(:K_S0).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    # nominal 4C fit on the full final state (pi+ pi- K+ K+ K- pi- K_S0)
    .kinematic_fit([:pip, :pim, :pim, :kp, :kp, :km, :K_S0]) {
        nominal
        constrain_four_momentum
        invariant_mass_of(:kp, :km).larger_than(1.2)    # phi   veto (Ds- -> K+K-pi-)
        invariant_mass_of(:kp, :pim).larger_than(1.0)   # K*(892) veto
        chi2_cut 200
    }

alg3.with_decay_card(decay_card_mode3).apply(sel3)

# ---------------------------------------------------------------------
# Mode IV : e+e- -> pi+pi- (Ds+ -> K_S0 K+) (Ds- -> K_S0 K-)
# ---------------------------------------------------------------------
alg_name_4 = "DsDsbarMode4"
alg4 = Algorithm.new(alg_name_4)
alg4.set_header(["#{alg_name_4}Alg/#{alg_name_4}.h"])
    .set_constant({"ECMS" => [:double, 4.6]})
    .note(:energy_point_handling,
          "15 c.m. energy points (4.4156-4.9509 GeV); ECMS resolved per dataset at execute_on")
    .note(:ks_mass_window,
          "each K_S0 candidate restricted to 0.487 < M(pi+pi-) < 0.511 GeV/c^2")
    .note(:ks_decay_length,
          "each K_S0 decay length required to be greater than twice its resolution")
    .note(:ds_pairing_criterion,
          "the (Ds+, Ds-) candidate pairing is chosen by minimising |0.5*(M(Ds+)+M(Ds-)) - M_PDG(Ds)|")
    .note(:recoil_mass_window,
          "|M(Ds) - M_PDG(Ds)| < 15 MeV/c^2 applied on the recoiling Ds in the RM(pi+pi-Ds) fit")

sel4 = Selection.new
sel4.select_track {
        cos_theta 0.93
        Vz        20.0
        Vr        10.0
        nChrp     ">=3"
        nChrn     ">=3"
        nNet      "==0"
    }
    .pid(method: :probability) {
        prob_cut 0.001
        identify :pion, against: [:kaon]
        identify :kaon, against: [:pion]
        nkp ">=1"
        nkm ">=1"
    }
    # two independent K_S0 -> pi+pi- candidates
    .secondary_vertex_fit([:pip, :pim]) {
        build_virtual_particle(:K_S0).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    .secondary_vertex_fit([:pip, :pim]) {
        build_virtual_particle(:K_S0).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    # nominal 4C fit on the full final state (pi+ pi- K+ K- K_S0 K_S0)
    .kinematic_fit([:pip, :pim, :kp, :km, :K_S0, :K_S0]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
    }

alg4.with_decay_card(decay_card_mode4).apply(sel4)

# =====================================================================
### Run on real data + inclusive MC + exclusive MC ###
# =====================================================================
root_files_1 = alg1.execute_on(scan_points + scan_incMC + exMC_mode1)
root_files_2 = alg2.execute_on(scan_points + scan_incMC + exMC_mode2)
root_files_3 = alg3.execute_on(scan_points + scan_incMC + exMC_mode3)
root_files_4 = alg4.execute_on(scan_points + scan_incMC + exMC_mode4)