### Dataset description ###
# Five BESIII energy points (BOSS 706): 4.628, 4.641, 4.661, 4.681, 4.698 GeV
data_points = [
  DatasetManager.real_data.find("706_4620"),   # 4.628 GeV
  DatasetManager.real_data.find("706_4640"),   # 4.641 GeV
  DatasetManager.real_data.find("706_4660"),   # 4.661 GeV
  DatasetManager.real_data.find("706_4680"),   # 4.681 GeV
  DatasetManager.real_data.find("706_4700"),   # 4.698 GeV
]
incMC_points = [
  DatasetManager.inclusive_mc.find("706_4620"),
  DatasetManager.inclusive_mc.find("706_4640"),
  DatasetManager.inclusive_mc.find("706_4660"),
  DatasetManager.inclusive_mc.find("706_4680"),
  DatasetManager.inclusive_mc.find("706_4700"),
]

# Decay card for e+e- -> K+ Ds- anti-D*0, anti-D*0 -> anti-D0 pi0, anti-D0 -> K+ pi-, pi0 -> gamma gamma
# Mode I: Ds- -> K+ K- pi-
decay_card_modeI = <<~DECAYCARD
    Decay psi(4260)
    1.0000 K+ D_s- anti-D*0 PHSP;
    Enddecay

    Decay anti-D*0
    1.0000 anti-D0 pi0 PHSP;
    Enddecay

    Decay anti-D0
    1.0000 K+ pi- PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay D_s-
    1.0000 K+ K- pi- PHSP;
    Enddecay

    End
DECAYCARD

# Mode II: Ds- -> K_S0 K-, K_S0 -> pi+ pi-
decay_card_modeII = <<~DECAYCARD
    Decay psi(4260)
    1.0000 K+ D_s- anti-D*0 PHSP;
    Enddecay

    Decay anti-D*0
    1.0000 anti-D0 pi0 PHSP;
    Enddecay

    Decay anti-D0
    1.0000 K+ pi- PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay D_s-
    1.0000 K_S0 K- PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# 500k-event exclusive MC at each of the five energy points, one sample per Ds- decay mode
exMCs_modeI = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_KpDsmissDstar0_modeI"
  config.events        = 500000
  config.decay_card    = decay_card_modeI
  config.cross_section = :default
end

exMCs_modeII = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_KpDsmissDstar0_modeII"
  config.events        = 500000
  config.decay_card    = decay_card_modeII
  config.cross_section = :default
end

### Event selection (BOSS) ###
# Mode I : Ds- -> K+ K- pi-
alg_modeI = Algorithm.new("KpDsMissDstar0ModeI")
alg_modeI.set_header(["KpDsMissDstar0ModeIAlg/KpDsMissDstar0ModeI.h"])
         .set_constant({"ECMS" => [:double, 4.662]})
         .note(:ecms_multi_energy, "Five data points (4.628/4.641/4.661/4.681/4.698 GeV) share one algorithm; ECMS must be taken per dataset at execute_on time")
         .note(:ds_mass_window, "Mode I: Ds- candidates kept in mass window 1.955-1.980 GeV/c^2 after partial reconstruction (no DSL primitive)")
         .note(:dalitz_region, "Mode I: Ds- candidates restricted to Dalitz regions M(K+K-) < 1.05 GeV/c^2 (phi) or 0.850 < M(K+pi-) < 0.930 GeV/c^2 (K*(892))")

sel_modeI = Selection.new
sel_modeI.select_track {
    cos_theta 0.93   # |cos(theta)| < 0.93
    Vz        10.0   # |Vz| < 10 cm
    Vr        1.0    # Vr < 1 cm
    nChrp     ">=2"  # at least two positive tracks
    nChrn     "==2"  # exactly two negative tracks
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]   # K+ and K- (charge-conjugation shorthand)
    nkp ">=2"   # at least two K+
    nkm "==1"   # exactly one K-
  }
  .remove([:kp <= :chrgp, :km <= :chrgn])      # drop identified kaons from the charged lists
  .assign({:chrgp => :pip, :chrgn => :pim})    # remaining tracks assigned as pions
  # Partial reconstruction: tag bachelor K+ and Ds-, leave anti-D*0 missing
  .partial_miss([3]) {                         # recID 3 = anti-D*0 (its whole subtree is missing)
    best_combination_by_mass :D_s, 1.968       # pick Ds- combination closest to nominal 1.968 GeV/c^2
    require_recoil_mass 1.990, 2.027           # RM(K+ Ds-) in [1.990, 2.027] GeV/c^2
  }

alg_modeI.with_decay_card(decay_card_modeI).apply(sel_modeI)
root_files_modeI = alg_modeI.execute_on(data_points + incMC_points + exMCs_modeI)

# Mode II : Ds- -> K_S0 K-, K_S0 -> pi+ pi-
alg_modeII = Algorithm.new("KpDsMissDstar0ModeII")
alg_modeII.set_header(["KpDsMissDstar0ModeIIAlg/KpDsMissDstar0ModeII.h"])
          .set_constant({"ECMS" => [:double, 4.662]})
          .note(:ecms_multi_energy, "Five data points (4.628/4.641/4.661/4.681/4.698 GeV) share one algorithm; ECMS must be taken per dataset at execute_on time")
          .note(:ks_mass_window, "Mode II: K_S0 built by the secondary-vertex fit and required in mass window 0.485-0.511 GeV/c^2 (no DSL primitive)")
          .note(:ds_mass_window, "Mode II: Ds- candidates kept in mass window 1.955-1.985 GeV/c^2 after partial reconstruction (no DSL primitive)")

sel_modeII = Selection.new
sel_modeII.select_track {
    cos_theta 0.93   # |cos(theta)| < 0.93
    Vz        10.0   # |Vz| < 10 cm
    Vr        1.0    # Vr < 1 cm
    nChrp     "==2"  # exactly two positive tracks
    nChrn     "==2"  # exactly two negative tracks
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]   # K+ and K-
    nkp "==1"   # exactly one K+
    nkm "==1"   # exactly one K-
  }
  .remove([:kp <= :chrgp, :km <= :chrgn])      # drop identified kaons from the charged lists
  .assign({:chrgp => :pip, :chrgn => :pim})    # remaining tracks assigned as pions
  .secondary_vertex_fit([:pip, :pim]) {        # combine pi+ pi- into a K_S0
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # Partial reconstruction: tag bachelor K+ and Ds-, leave anti-D*0 missing
  .partial_miss([3]) {                         # recID 3 = anti-D*0
    best_combination_by_mass :D_s, 1.968       # pick Ds- combination closest to nominal 1.968 GeV/c^2
    require_recoil_mass 1.990, 2.027           # RM(K+ Ds-) in [1.990, 2.027] GeV/c^2
  }

alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)
root_files_modeII = alg_modeII.execute_on(data_points + incMC_points + exMCs_modeII)