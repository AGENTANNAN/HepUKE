# ==================== Dataset preparation ====================
# ψ(2S) resonance (3.686 GeV) real data and its inclusive MC
psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Off-resonance points used for continuum / QED background studies
data_3773  = DatasetManager.real_data.find("712_3773")   # 3.773 GeV
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")
data_3650  = DatasetManager.real_data.find("709_3650")   # 3.65 GeV
incMC_3650 = DatasetManager.inclusive_mc.find("709_3650")

# ==================== Decay cards (EvtGen format) ====================
# ---- Signal modes (all share the γγ K+K− final state) ----
# ψ' → K+K−π0,  π0 → γγ
decay_card_KKpi0 = <<~DECAYCARD
  Decay psi(2S)
  1.0000 K+ K- pi0 PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# ψ' → K+K−η,  η → γγ
decay_card_KKeta = <<~DECAYCARD
  Decay psi(2S)
  1.0000 K+ K- eta PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# ψ' → ηφ,  η → γγ,  φ → K+K−
decay_card_etaphi = <<~DECAYCARD
  Decay psi(2S)
  1.0000 eta phi PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay

  Decay phi
  1.0000 K+ K- VSS;
  Enddecay

  End
DECAYCARD

# ψ' → π0φ,  π0 → γγ,  φ → K+K−
decay_card_pi0phi = <<~DECAYCARD
  Decay psi(2S)
  1.0000 pi0 phi PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  Decay phi
  1.0000 K+ K- VSS;
  Enddecay

  End
DECAYCARD

# ---- Continuum background modes (top mother psi(4260), KKMC convention) ----
# γ* → K+K−π0
decay_card_cont_KKpi0 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 K+ K- pi0 PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# γ* → ηφ
decay_card_cont_etaphi = <<~DECAYCARD
  Decay psi(4260)
  1.0000 eta phi PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay

  Decay phi
  1.0000 K+ K- VSS;
  Enddecay

  End
DECAYCARD

# ---- Background modes at the ψ(2S) ----
# ψ' → γχc2,  χc2 → K+K−π0
decay_card_bkg_chic2_KKpi0 = <<~DECAYCARD
  Decay psi(2S)
  1.0000 gamma chi_c2 PHSP;
  Enddecay

  Decay chi_c2
  1.0000 K+ K- pi0 PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# ψ' → γχc2,  χc2 → K+K−η
decay_card_bkg_chic2_KKeta = <<~DECAYCARD
  Decay psi(2S)
  1.0000 gamma chi_c2 PHSP;
  Enddecay

  Decay chi_c2
  1.0000 K+ K- eta PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# γγ_FSR K+K− background (two FSR photons + K+K−)
decay_card_bkg_ggFSR = <<~DECAYCARD
  Decay psi(2S)
  1.0000 K+ K- gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# ==================== Exclusive MC samples ====================
# 100k-event signal MC, one per signal decay mode
exMC_KKpi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_KKpi0"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_KKpi0
  config.cross_section   = :default
end

exMC_KKeta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_KKeta"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_KKeta
  config.cross_section   = :default
end

exMC_etaphi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_etaphi"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_etaphi
  config.cross_section   = :default
end

exMC_pi0phi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_pi0phi"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_pi0phi
  config.cross_section   = :default
end

# Continuum background: same card run over the two off-resonance points
exMC_cont_KKpi0 = DatasetManager.create_exclusive_mc_for([data_3773, data_3650]) do |config|
  config.sample_name   = "continuum_to_KKpi0"
  config.events        = 100_000
  config.decay_card    = decay_card_cont_KKpi0
  config.cross_section = :default
end

exMC_cont_etaphi = DatasetManager.create_exclusive_mc_for([data_3773, data_3650]) do |config|
  config.sample_name   = "continuum_to_etaphi"
  config.events        = 100_000
  config.decay_card    = decay_card_cont_etaphi
  config.cross_section = :default
end

# Background modes at the ψ(2S)
exMC_bkg_chic2_KKpi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "bkg_chic2_to_KKpi0"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_bkg_chic2_KKpi0
  config.cross_section   = :default
end

exMC_bkg_chic2_KKeta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "bkg_chic2_to_KKeta"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_bkg_chic2_KKeta
  config.cross_section   = :default
end

exMC_bkg_ggFSR = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "bkg_ggFSR_KK"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_bkg_ggFSR
  config.cross_section   = :default
end

# ==================== Event selection (one common chain for all modes) ====================
common_selection = Selection.new
  .select_track {
      cos_theta 0.93        # |cosθ| < 0.93
      Vz        10.0        # |Vz| < 10 cm
      Vr        1.0         # Vr < 1 cm
      nChrp     "==1"       # exactly one positively charged track
      nChrn     "==1"       # exactly one negatively charged track
      nNet      "==0"       # net charge zero
  }
  .select_photon {
      tdc_emc_start     0
      tdc_emc_end       14
      angle_to_track    10.0     # ≥ 10° away from any charged track
      energyThreshold_b 0.025    # E > 25 MeV, barrel (|cosθ| < 0.80)
      energyThreshold_e 0.050    # E > 50 MeV, endcaps (0.86 < |cosθ| < 0.92)
      nGam              ">=2"    # at least two photons (2–10 applied later in ROOT)
  }
  .pid(method: :probability) {
      prob_cut 0.001
      identify :kaon, against: [:pion, :proton]   # both tracks as K+ / K−
      nkp      "==1"
      nkm      "==1"
  }
  .kinematic_fit([:gamma, :gamma, :kp, :km]) {
      nominal
      constrain_four_momentum   # 4C fit to the initial e+e− four-momentum (γγ K+K− hypothesis)
      chi2_cut 200              # loose cut in BOSS; χ² ≤ 20 applied at the ROOT level
  }

# ==================== Algorithms (one per signal decay mode) ====================
alg_name_KKpi0 = "PsipKKpi0"
alg_KKpi0 = Algorithm.new(alg_name_KKpi0)
alg_KKpi0.set_header(["#{alg_name_KKpi0}Alg/#{alg_name_KKpi0}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})
         .set_alias({"std::vector<double>" => "Vdouble"})
alg_KKpi0.note(:background_veto,
  "γγ_FSR K+K− QED background generated as a direct 4-body psi(2S) -> K+ K- gamma gamma
   phase-space card; the exact FSR photon energy/angular spectrum is not modelled by the generator")
alg_KKpi0.with_decay_card(decay_card_KKpi0).apply(common_selection.dup)

alg_name_KKeta = "PsipKKeta"
alg_KKeta = Algorithm.new(alg_name_KKeta)
alg_KKeta.set_header(["#{alg_name_KKeta}Alg/#{alg_name_KKeta}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})
         .set_alias({"std::vector<double>" => "Vdouble"})
alg_KKeta.with_decay_card(decay_card_KKeta).apply(common_selection.dup)

alg_name_etaphi = "PsipEtaPhi"
alg_etaphi = Algorithm.new(alg_name_etaphi)
alg_etaphi.set_header(["#{alg_name_etaphi}Alg/#{alg_name_etaphi}.h"])
          .set_constant({"ECMS" => [:double, 3.686]})
          .set_alias({"std::vector<double>" => "Vdouble"})
alg_etaphi.with_decay_card(decay_card_etaphi).apply(common_selection.dup)

alg_name_pi0phi = "PsipPi0Phi"
alg_pi0phi = Algorithm.new(alg_name_pi0phi)
alg_pi0phi.set_header(["#{alg_name_pi0phi}Alg/#{alg_name_pi0phi}.h"])
          .set_constant({"ECMS" => [:double, 3.686]})
          .set_alias({"std::vector<double>" => "Vdouble"})
alg_pi0phi.with_decay_card(decay_card_pi0phi).apply(common_selection.dup)

# ==================== Execution ====================
# Signal modes: ψ(2S) data + inclusive MC + signal MC
alg_KKpi0.execute_on([psip_data, psip_incMC, exMC_KKpi0])
alg_KKeta.execute_on([psip_data, psip_incMC, exMC_KKeta])
alg_etaphi.execute_on([psip_data, psip_incMC, exMC_etaphi])
alg_pi0phi.execute_on([psip_data, psip_incMC, exMC_pi0phi])

# Radiative χc2 and γγ_FSR background modes at the ψ(2S) (same chain)
alg_KKpi0.execute_on([psip_data, psip_incMC,
                      exMC_bkg_chic2_KKpi0, exMC_bkg_chic2_KKeta, exMC_bkg_ggFSR])

# Continuum / QED background at the two off-resonance points
alg_KKpi0.execute_on([data_3773, incMC_3773, data_3650, incMC_3650,
                      *exMC_cont_KKpi0, *exMC_cont_etaphi])