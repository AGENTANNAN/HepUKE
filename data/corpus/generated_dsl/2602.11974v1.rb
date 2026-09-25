# =====================================================================
# Λc+ Λc− threshold scan (4.600–4.700 GeV, seven points)
#   signal        : Λc+ → p η'   (η' → π+π−γ)
#   normalization : Λc+ → p ω    (ω → π+π−π0, π0 → γγ)
# =====================================================================

### Dataset description ###
# Seven scan points from 4.600 to 4.700 GeV (real data + matching inclusive MC)
data_points = [
  DatasetManager.real_data.find("703_4600"),   # 4.600 GeV
  DatasetManager.real_data.find("706_4610"),   # 4.611 GeV
  DatasetManager.real_data.find("706_4620"),   # 4.628 GeV
  DatasetManager.real_data.find("706_4640"),   # 4.641 GeV
  DatasetManager.real_data.find("706_4660"),   # 4.661 GeV
  DatasetManager.real_data.find("706_4680"),   # 4.682 GeV
  DatasetManager.real_data.find("706_4700")    # 4.699 GeV
]

incMC_points = [
  DatasetManager.inclusive_mc.find("703_4600"),
  DatasetManager.inclusive_mc.find("706_4610"),
  DatasetManager.inclusive_mc.find("706_4620"),
  DatasetManager.inclusive_mc.find("706_4640"),
  DatasetManager.inclusive_mc.find("706_4660"),
  DatasetManager.inclusive_mc.find("706_4680"),
  DatasetManager.inclusive_mc.find("706_4700")
]

# ---- decay cards (EvtGen) ----
# Signal: e+e- → Λc+ Λc−, tag side Λc+ → p η' → p π+π−γ
decay_card_signal = <<~DECAYCARD
  Decay psi(4260)
  1.0000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay

  Decay Lambda_c+
  1.0000 p+ eta' PHSP;
  Enddecay

  Decay eta'
  1.0000 pi+ pi- gamma PHSP;
  Enddecay

  Decay anti-Lambda_c-
  1.0000 anti-p- K+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# Normalization: Λc+ → p ω → p π+π−π0, π0 → γγ
decay_card_norm = <<~DECAYCARD
  Decay psi(4260)
  1.0000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay

  Decay Lambda_c+
  1.0000 p+ omega PHSP;
  Enddecay

  Decay omega
  1.0000 pi+ pi- pi0 OMEGA_DALITZ;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  Decay anti-Lambda_c-
  1.0000 anti-p- K+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# ---- exclusive MC: 200k events per mode, one sample per energy point ----
exMC_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_lctop_etap"
  config.events        = 200000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

exMC_norm = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_lctop_omega"
  config.events        = 200000
  config.decay_card    = decay_card_norm
  config.cross_section = :default
end

### Event selection (BOSS) ###

# ================= Signal channel: Λc+ → p η' =================
alg_signal = Algorithm.new("LcToPEtap")
alg_signal.set_header(["LcToPEtapAlg/LcToPEtap.h"])
          .set_constant({"ECMS" => [:double, 4.600]})   # per-point ECMS supplied at job level
          .set_alias({"std::vector<double>" => "Vdouble"})
          .note(:ecms_scan, "the Λc+Λc− threshold scan runs over seven energy points (4.600-4.700 GeV); ECMS is set per scan point when the algorithm is executed on each dataset, the value here is the first point")

sel_signal = Selection.new
sel_signal
  .select_track {
    cos_theta 0.93   # |cosθ| < 0.93
    Vz        10.0   # |Vz| < 10 cm
    Vr        1.0    # Vr < 1 cm
    nChrp     ">=2"  # >= 2 positively charged tracks
    nChrn     ">=1"  # >= 1 negatively charged track
    nTot      ">=3"  # >= 3 charged tracks in total
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=1"   # >= 1 photon
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]  # protons vs K and π
    identify :pion,   against: [:kaon]         # pions vs K
    nprp ">=1"
    npip ">=1"
    npim ">=1"
  }
  # η' formed from π+π−γ by a mass-constrained (1-C) fit
  .kalman_kinematic_fit([:pip, :pim, :gamma]) {
    invariant_mass_of(:pip, :pim, :gamma).constrain_to_nominal_mass_of(:etap)
    invariant_mass_of(:pip, :pim, :gamma).between(0.90, 1.00)  # M(π+π−γ) window
    chi2_cut 200
    netap ">=1"
  }
  # nominal Λc+ mass constraint on p η'  (p π+π−γ), with Λ / K_S0 vetoes
  .kinematic_fit([:prp, :pip, :pim, :gamma]) {
    nominal
    invariant_mass_of(:prp, :pip, :pim, :gamma).constrain_to_nominal_mass_of(:"Lambda_c+")
    invariant_mass_of(:prp, :pim).out_of(1.10, 1.15)  # Λ veto  (M(pπ−))
    invariant_mass_of(:pip, :pim).out_of(0.48, 0.51)  # K_S0 veto (M(π+π−))
    chi2_cut 200
  }

alg_signal.with_decay_card(decay_card_signal).apply(sel_signal)
root_files_signal = alg_signal.execute_on(data_points + incMC_points + exMC_signal)

# ================= Normalization channel: Λc+ → p ω =================
alg_norm = Algorithm.new("LcToPOmega")
alg_norm.set_header(["LcToPOmegaAlg/LcToPOmega.h"])
        .set_constant({"ECMS" => [:double, 4.600]})
        .set_alias({"std::vector<double>" => "Vdouble"})
        .note(:ecms_scan, "ECMS is set per scan point (4.600-4.700 GeV) at job level, same as the signal channel")

sel_norm = Selection.new
sel_norm
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=2"
    nChrn     ">=1"
    nTot      ">=3"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=2"   # >= 2 photons
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    identify :pion,   against: [:kaon]
    nprp ">=1"
    npip ">=1"
    npim ">=1"
  }
  # π0 reconstructed from γγ via a mass-constrained (1-C) fit
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    invariant_mass_of(:gamma, :gamma).between(0.115, 0.150)  # M(γγ) window
    chi2_cut 200
    npi0 ">=1"
  }
  # nominal fit: M(π+π−π0) → ω, M(p π+π−π0) → Λc+, with Λ / K_S0 / Σ+ vetoes
  .kinematic_fit([:prp, :pip, :pim, :pi0]) {
    nominal
    invariant_mass_of(:pip, :pim, :pi0).constrain_to_nominal_mass_of(:omega)
    invariant_mass_of(:prp, :pip, :pim, :pi0).constrain_to_nominal_mass_of(:"Lambda_c+")
    invariant_mass_of(:pip, :pim, :pi0).between(0.73, 0.83)  # ω mass window
    invariant_mass_of(:prp, :pim).out_of(1.10, 1.15)  # Λ veto
    invariant_mass_of(:pip, :pim).out_of(0.48, 0.51)  # K_S0 veto
    invariant_mass_of(:prp, :pi0).out_of(1.17, 1.20)  # Σ+ veto (M(pπ0))
    chi2_cut 200
  }

alg_norm.with_decay_card(decay_card_norm).apply(sel_norm)
root_files_norm = alg_norm.execute_on(data_points + incMC_points + exMC_norm)