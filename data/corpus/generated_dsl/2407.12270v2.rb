# ===========================================================================
# Dataset preparation — single-tag Λc+ at threshold, 4.600 – 4.840 GeV
# ===========================================================================
# 11 energy points of the Λc threshold scan (4.600 – 4.840 GeV, 6.1 fb−1)
scan_points = %w[
  703_4600 706_4610 706_4620 706_4640 706_4660 706_4680 706_4700
  707_4740 707_4750 707_4780 707_4840
]
real_data = scan_points.map { |name| DatasetManager.real_data.find(name) }   # real data at the 11 points
incMC     = scan_points.map { |name| DatasetManager.inclusive_mc.find(name) } # corresponding inclusive MC

# Decay card — Mode I: Λc+ → Λ π+ η, Λ → p π−, η → γγ
decay_card_modeI = <<~DECAYCARD
  Decay psi(4260)
  1.0000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay

  Decay Lambda_c+
  1.0000 Lambda0 pi+ eta PHSP;
  Enddecay

  Decay anti-Lambda_c-
  1.0000 anti-p- K+ pi- PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi- HypWK;
  Enddecay

  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Decay card — Mode II: Λc+ → Λ π+ η, Λ → p π−, η → π+ π− π0, π0 → γγ
decay_card_modeII = <<~DECAYCARD
  Decay psi(4260)
  1.0000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay

  Decay Lambda_c+
  1.0000 Lambda0 pi+ eta PHSP;
  Enddecay

  Decay anti-Lambda_c-
  1.0000 anti-p- K+ pi- PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi- HypWK;
  Enddecay

  Decay eta
  1.0000 pi+ pi- pi0 PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# 1M-event exclusive MC for each of the two modes, one sample per energy point
exMCs_modeI = DatasetManager.create_exclusive_mc_for(real_data) do |config|
  config.sample_name   = "LcToLambdaPiEta_gg"          # Λc+ → Λ π+ η, η → γγ
  config.events        = 1_000_000
  config.decay_card    = decay_card_modeI
  config.cross_section = :default
end

exMCs_modeII = DatasetManager.create_exclusive_mc_for(real_data) do |config|
  config.sample_name   = "LcToLambdaPiEta_pipimpi0"    # Λc+ → Λ π+ η, η → π+ π− π0
  config.events        = 1_000_000
  config.decay_card    = decay_card_modeII
  config.cross_section = :default
end

all_datasets = real_data + incMC + exMCs_modeI + exMCs_modeII

# ===========================================================================
# Mode I — Λc+ → Λ π+ η, η → γγ
# ===========================================================================
alg_name_I = "LcToLambdaPiEtaGG"
alg_modeI  = Algorithm.new(alg_name_I)
alg_modeI.set_header(["#{alg_name_I}Alg/#{alg_name_I}.h"])
         .set_constant({"ECMS" => [:double, 4.600]})   # threshold point (see note on the scan)
         .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeI = Selection.new
  .select_track {                 # charged track quality + multiplicity (at least 2 positive, 1 negative)
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=2"
    nChrn     ">=1"
  }
  .select_photon {                # photon selection (≥2 photons)
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=2"
  }
  .pid(method: :probability) {    # probability PID: p vs K/π and π vs K
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    identify :pion,   against: [:kaon]
    nprp ">=1"
    npip ">=1"
  }
  .remove([:prp <= :chrgp])       # remove identified protons from the charged-pion candidates
  .secondary_vertex_fit([:prp, :pim]) {   # build Λ from p π− (minimise mass difference)
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # reconstruct η from γγ (1-C mass constraint)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 20
    neta ">=1"
  }
  .kinematic_fit([:Lambda, :pip, :pip, :eta]) {   # 4C fit — list carries the duplicated π+ as written
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_modeI
  .note(:best_candidate_selection,
    "keep the Λc+ candidate with the minimum |ΔE| per event, within -0.1 < ΔE < 0.1 GeV; " \
    "ΔE is formed from the fitted Λc+ four-momentum after the 4C fit, so the window is applied at the ROOT level")
  .note(:ecms_energy_scan,
    "ECMS is listed at the 4.600 GeV threshold point only; the same selection is run over the 11 scan points " \
    "4.600-4.840 GeV, so the 4C constraint energy must follow the per-run measured beam energy")

alg_modeI.with_decay_card(decay_card_modeI).apply(sel_modeI)
alg_modeI.execute_on(all_datasets)

# ===========================================================================
# Mode II — Λc+ → Λ π+ η, η → π+ π− π0, π0 → γγ
# ===========================================================================
alg_name_II = "LcToLambdaPiEtaPipPimPi0"
alg_modeII  = Algorithm.new(alg_name_II)
alg_modeII.set_header(["#{alg_name_II}Alg/#{alg_name_II}.h"])
          .set_constant({"ECMS" => [:double, 4.600]})
          .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeII = Selection.new
  .select_track {                 # same track quality cuts, ≥2 positive and ≥2 negative tracks
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=2"
    nChrn     ">=2"
  }
  .select_photon {                # same photon cuts (≥2 photons)
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=2"
  }
  .pid(method: :probability) {    # p vs K/π, π vs K, plus the extra ≥1 π− requirement
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    identify :pion,   against: [:kaon]
    nprp ">=1"
    npip ">=1"
    npim ">=1"
  }
  .remove([:prp <= :chrgp])       # remove identified protons from the charged-pion candidates
  .secondary_vertex_fit([:prp, :pim]) {   # build Λ from p π− (minimise mass difference)
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # reconstruct π0 from γγ (1-C mass constraint)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200
    npi0 ">=1"
  }
  .kinematic_fit([:Lambda, :pip, :pip, :pim, :pi0]) {   # 4C fit on Λ π+ π+ π− π0 (Λ π+ η)
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_modeII
  .note(:best_candidate_selection,
    "keep the Λc+ candidate with the minimum |ΔE| per event, within -0.1 < ΔE < 0.1 GeV; " \
    "ΔE is formed from the fitted Λc+ four-momentum after the 4C fit, so the window is applied at the ROOT level")
  .note(:ecms_energy_scan,
    "ECMS is listed at the 4.600 GeV threshold point only; the same selection is run over the 11 scan points " \
    "4.600-4.840 GeV, so the 4C constraint energy must follow the per-run measured beam energy")

alg_modeII.with_decay_card(decay_card_modeII).apply(sel_modeII)
alg_modeII.execute_on(all_datasets)