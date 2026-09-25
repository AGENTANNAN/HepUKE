# =====================================================================
# psi(3686) -> Xi(1530)0 anti-Xi(1530)0  and  psi(3686) -> Xi(1530)0 anti-Xi0
# Single-baryon tagging: tag = Xi(1530)0 -> Xi- pi+ (or charge conjugate);
# the opposite side is the signal, extracted from the recoil mass M_recoil(Xi- pi+).
# =====================================================================

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")      # psi(3686) 3.686 GeV real data
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # matching inclusive MC

# ---------- Mode I decay card: psi(3686) -> Xi(1530)0 anti-Xi(1530)0 ----------
decay_card_modeI = <<~DECAYCARD
    Decay psi(2S)
    1.0000 Xi(1530)0 anti-Xi(1530)0 PHSP;
    Enddecay

    Decay Xi(1530)0
    1.0000 Xi- pi+ PHSP;
    Enddecay

    Decay anti-Xi(1530)0
    1.0000 anti-Xi+ pi- PHSP;
    Enddecay

    Decay Xi-
    1.0000 Lambda0 pi- HypWK;
    Enddecay

    Decay anti-Xi+
    1.0000 anti-Lambda0 pi+ HypWK;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

# ---------- Mode II decay card: psi(3686) -> Xi(1530)0 anti-Xi0 ----------
decay_card_modeII = <<~DECAYCARD
    Decay psi(2S)
    1.0000 Xi(1530)0 anti-Xi0 PHSP;
    Enddecay

    Decay Xi(1530)0
    1.0000 Xi- pi+ PHSP;
    Enddecay

    Decay Xi-
    1.0000 Lambda0 pi- HypWK;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- HypWK;
    Enddecay

    Decay anti-Xi0
    1.0000 anti-Lambda0 pi0 PHSP;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# ---------- Exclusive MC for both modes ----------
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psi2S_Xi1530_Xi1530bar"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psi2S_Xi1530_Xibar0"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

### Event selection (BOSS) ###

# ============================== Mode I ==============================
# psi(3686) -> Xi(1530)0 anti-Xi(1530)0 ; tag Xi(1530)0 -> Xi- pi+ , signal from recoil
alg_modeI = Algorithm.new("Xi1530TagModeI")
alg_modeI.set_header(["Xi1530TagModeIAlg/Xi1530TagModeI.h"])
         .set_constant({"ECMS" => [:double, 3.686]})

sel_modeI = Selection.new
  .select_track {                 # charged track selection
    cos_theta 0.93                # |cos(theta)| < 0.93
    Vz        100.0               # |Vz| < 100 cm
    Vr        10.0                # Vr < 10 cm
    nChrp     ">=3"               # at least 3 positive tracks
    nChrn     ">=3"               # at least 3 negative tracks
  }
  .pid(method: :probability) {    # 1st PID pass: protons / anti-protons
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp ">=1"
    nprm ">=1"
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])   # remove proton tracks from charged lists
  .pid(method: :probability) {    # 2nd PID pass: pions
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
  }
  .secondary_vertex_fit([:prp, :pim]) {       # Lambda -> p pi-
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {       # anti-Lambda -> anti-p pi+
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # Partial reconstruction of the tag Xi(1530)0 (recID 1) -> Xi- pi+ -> Lambda pi- pi+ ;
  # the opposite baryon is inferred from the recoil four-momentum (no 4C fit is performed).
  .partial_rec([1]) {
    best_combination_by_mass :Xi_1530_0, 1.5318   # pick combination closest to M(Xi(1530)0)
  }

alg_modeI
  .note(:decay_length, "Xi- cascade required to have decay length > 0 (vertex displaced from the IP) in the Xi(1530)0 tag reconstruction")
  .note(:mass_window, "Lambda mass window |M(p pi-) - M(Lambda)| < 5 MeV and Xi- mass window |M(Lambda pi-) - M(Xi-)| < 8 MeV applied after the secondary-vertex / cascade reconstruction")
  .note(:background_veto, "psi(3686) -> pi+ pi- J/psi background suppressed by requiring |M_recoil(pi+ pi-) - M(J/psi)| > 4 MeV; recoil-mass veto applied event-by-event")
  .with_decay_card(decay_card_modeI)
  .apply(sel_modeI)

# ============================== Mode II ==============================
# psi(3686) -> Xi(1530)0 anti-Xi0 ; tag Xi(1530)0 -> Xi- pi+ , signal from recoil,
# with an additional pi0 (from anti-Xi0 -> anti-Lambda pi0) reconstructed from two photons
alg_modeII = Algorithm.new("Xi1530TagModeII")
alg_modeII.set_header(["Xi1530TagModeIIAlg/Xi1530TagModeII.h"])
          .set_constant({"ECMS" => [:double, 3.686]})

sel_modeII = Selection.new
  .select_track {                 # charged track selection
    cos_theta 0.93
    Vz        100.0
    Vr        10.0
    nChrp     ">=3"
    nChrn     ">=3"
  }
  .select_photon {                # photons for pi0 -> gamma gamma
    tdc_emc_start 0
    tdc_emc_end   14
    angle_to_track 10.0
    energyThreshold_b 0.025       # > 25 MeV in the barrel
    energyThreshold_e 0.050       # > 50 MeV in the endcap
    nGam ">=2"                    # at least two photons
  }
  .pid(method: :probability) {    # 1st PID pass: protons / anti-protons
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp ">=1"
    nprm ">=1"
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  .pid(method: :probability) {    # 2nd PID pass: pions
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
  }
  .secondary_vertex_fit([:prp, :pim]) {       # Lambda -> p pi-
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:prm, :pip]) {       # anti-Lambda -> anti-p pi+
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # 1C fit: gamma gamma -> pi0
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"                                # at least one pi0 candidate
  }
  # Partial reconstruction of the tag Xi(1530)0 (recID 1); signal side from recoil.
  .partial_rec([1]) {
    best_combination_by_mass :Xi_1530_0, 1.5318
  }

alg_modeII
  .note(:decay_length, "Xi- cascade required to have decay length > 0 (vertex displaced from the IP) in the Xi(1530)0 tag reconstruction")
  .note(:mass_window, "Lambda mass window |M(p pi-) - M(Lambda)| < 5 MeV and Xi- mass window |M(Lambda pi-) - M(Xi-)| < 8 MeV applied after the secondary-vertex / cascade reconstruction")
  .note(:background_veto, "psi(3686) -> pi+ pi- J/psi background suppressed by requiring |M_recoil(pi+ pi-) - M(J/psi)| > 4 MeV; recoil-mass veto applied event-by-event")
  .with_decay_card(decay_card_modeII)
  .apply(sel_modeII)

# ============================== Execution ==============================
root_files_modeI  = alg_modeI.execute_on([psip_data, psip_incMC, exMC_modeI])
root_files_modeII = alg_modeII.execute_on([psip_data, psip_incMC, exMC_modeII])