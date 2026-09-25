# =====================================================================
#  e+e- -> Omega- Omega+   (single-baryon tagging)
#  BOSS part: dataset preparation + event selection
# =====================================================================

### ----------------------------- Dataset preparation ----------------------------- ###
# Real-data scan points covering 3.81 - 4.78 GeV
scan_data_names = %w[
  703_3810 703_3872 703_3900 703_4009 703_4090
  703_4180 703_4190 703_4200 703_4210 703_4220 703_4230
  703_4237 703_4245 703_4246 703_4260 703_4270 703_4280
  703_4310 703_4360 703_4390 703_4420 703_4470 703_4530
  703_4575 703_4600
  705_4130 705_4160 705_4290 705_4315 705_4340 705_4380 705_4400 705_4440
  706_4610 706_4620 706_4640 706_4660 706_4680 706_4700
  707_4740 707_4750 707_4780
]
scan_data = scan_data_names.map { |name| DatasetManager.real_data.find(name) }

# Matching inclusive MC samples at the same energy points
scan_incmc_names = %w[
  703_4009 703_4180 703_4190 703_4200 703_4210 703_4220 703_4230
  703_4237 703_4246 703_4260 703_4270 703_4280 703_4360 703_4420 703_4600
  705_4130 705_4160
  706_4610 706_4620 706_4640 706_4660 706_4680 706_4700
  707_4740 707_4750 707_4780
]
scan_incMC = scan_incmc_names.map { |name| DatasetManager.inclusive_mc.find(name) }

# ---------------------------------------------------------------------
# Decay cards (EvtGen).  psi(4260) -> Omega- anti-Omega+, with the
# corresponding Omega decay written out explicitly for each tag mode.
# ---------------------------------------------------------------------
# Omega- tag mode:  Omega- -> K- Lambda,  Lambda -> p+ pi-
decay_card_omega_minus = <<~DECAYCARD
  Decay psi(4260)
  1.0000 Omega- anti-Omega+ PHSP;
  Enddecay

  Decay Omega-
  1.0000 K- Lambda0 PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi- HypWK;
  Enddecay

  Decay anti-Omega+
  1.0000 K+ anti-Lambda0 PHSP;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+ HypWK;
  Enddecay

  End
DECAYCARD

# Omega+ tag mode: charge conjugate,  Omega+ -> K+ anti-Lambda,
# anti-Lambda -> anti-p- pi+
decay_card_omega_plus = <<~DECAYCARD
  Decay psi(4260)
  1.0000 Omega- anti-Omega+ PHSP;
  Enddecay

  Decay Omega-
  1.0000 K- Lambda0 PHSP;
  Enddecay

  Decay Lambda0
  1.0000 p+ pi- HypWK;
  Enddecay

  Decay anti-Omega+
  1.0000 K+ anti-Lambda0 PHSP;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+ HypWK;
  Enddecay

  End
DECAYCARD

# ---------------------------------------------------------------------
# Exclusive MC, 100k events per tag mode
# ---------------------------------------------------------------------
data_ref = DatasetManager.real_data.find("703_4260")   # 4.260 GeV reference point

exmc_omega_minus = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_omega_minus_tag"
  config.related_dataset = data_ref
  config.events          = 100000
  config.decay_card      = decay_card_omega_minus
  config.cross_section   = :default
end

exmc_omega_plus = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_omega_plus_tag"
  config.related_dataset = data_ref
  config.events          = 100000
  config.decay_card      = decay_card_omega_plus
  config.cross_section   = :default
end

### ------------------------- Event selection (BOSS) ------------------------- ###
# ---------------------------------------------------------------------
# Mode I : Omega- tag   Omega- -> K- Lambda, Lambda -> p+ pi-
# ---------------------------------------------------------------------
alg_name_minus = "OmegaMinusTag"
alg_omega_minus = Algorithm.new(alg_name_minus)
alg_omega_minus.set_header(["#{alg_name_minus}Alg/#{alg_name_minus}.h"])
               .set_constant({ "ECMS" => [:double, 4.260] })
               .set_alias({ "std::vector<double>" => "Vdouble" })

sel_omega_minus = Selection.new
sel_omega_minus
  .select_track {                       # 1 positive + 2 negative charged tracks
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==1"
    nChrn     "==2"
    nNet      "==-1"
  }
  .pid(method: :probability) {          # probability PID, prob > 0.001
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]   # p+
    identify :kaon,   against: [:pion, :proton] # K-
    identify :pion,   against: [:kaon, :proton] # pi-
    nprp "==1"
    nkm  "==1"
    npim "==1"
  }
  .secondary_vertex_fit([:prp, :pim]) { # Lambda -> p+ pi-
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .partial_miss([2]) {                  # miss the anti-Omega+ (recID 2); everything else tagged
    best_combination_by_mass :"Omega-", 1.67245
    require_recoil_mass 1.64, 1.70
  }

alg_omega_minus.with_decay_card(decay_card_omega_minus).apply(sel_omega_minus)

# ---------------------------------------------------------------------
# Mode II : Omega+ tag   Omega+ -> K+ anti-Lambda, anti-Lambda -> anti-p- pi+
# ---------------------------------------------------------------------
alg_name_plus = "OmegaPlusTag"
alg_omega_plus = Algorithm.new(alg_name_plus)
alg_omega_plus.set_header(["#{alg_name_plus}Alg/#{alg_name_plus}.h"])
              .set_constant({ "ECMS" => [:double, 4.260] })
              .set_alias({ "std::vector<double>" => "Vdouble" })

sel_omega_plus = Selection.new
sel_omega_plus
  .select_track {                       # 2 positive + 1 negative charged tracks
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"
    nChrn     "==1"
    nNet      "==1"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]   # anti-p-
    identify :kaon,   against: [:pion, :proton] # K+
    identify :pion,   against: [:kaon, :proton] # pi+
    nprm "==1"
    nkp  "==1"
    npip "==1"
  }
  .secondary_vertex_fit([:prm, :pip]) { # anti-Lambda -> anti-p- pi+
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .partial_miss([1]) {                  # miss the Omega- (recID 1); everything else tagged
    best_combination_by_mass :"anti-Omega+", 1.67245
    require_recoil_mass 1.64, 1.70
  }

alg_omega_plus.with_decay_card(decay_card_omega_plus).apply(sel_omega_plus)

### ------------------------------- Execution ------------------------------- ###
root_files_omega_minus = alg_omega_minus.execute_on(scan_data + scan_incMC + [exmc_omega_minus])
root_files_omega_plus  = alg_omega_plus.execute_on(scan_data + scan_incMC + [exmc_omega_plus])