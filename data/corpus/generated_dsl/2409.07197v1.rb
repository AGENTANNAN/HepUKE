# =====================================================================
# ψ(3770) → D0 anti-D0 : double-tag measurement of the CP-even fractions
# F_+ of D0 → π+π−π0 and D0 → K+K−π0
# (BOSS part: dataset preparation + event selection up to the final 4C fit)
# Tag side = recoil anti-D0 reconstructed from the pre-stored DTag candidates.
# =====================================================================

# ---------------- datasets ----------------
data_3773  = DatasetManager.real_data.find("712_3773")     # 3.773 GeV ψ(3770) real data
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")  # inclusive MC at 3.773 GeV

# ---------------- signal decay cards (EvtGen) ----------------
# ψ(3770) → D0 anti-D0 ; D0 → π+π−π0
decay_card_pipip0 = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 pi+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# ψ(3770) → D0 anti-D0 ; D0 → K+K−π0
decay_card_KKp0 = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 K+ K- pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# ---------------- exclusive MC for the two signal modes ----------------
exMC_pipip0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0_pipip0"
  config.related_dataset = data_3773
  config.events          = 200_000
  config.decay_card      = decay_card_pipip0
  config.cross_section   = :default
end

exMC_KKp0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0_KKp0"
  config.related_dataset = data_3773
  config.events          = 100_000
  config.decay_card      = decay_card_KKp0
  config.cross_section   = :default
end

# ---------------- tag modes: recoil anti-D0 in ~12 CP / CP-mixed modes ----------------
tag_modes = [
  :D0toKK,          # K+K-
  :D0toPiPi,        # π+π-
  :D0toKsPiPi,      # K_S ππ
  :D0toKPi,         # Kπ
  :D0toKsPi0,       # K_S π0
  :D0toKsPi0Pi0,    # K_S π0π0
  :D0toKsOmega,     # K_S ω
  :D0toKsEta,       # K_S η
  :D0toKsEtaPrime,  # K_S η'
  :D0toKlPi0,       # K_L π0
  :D0toKlPi0Pi0,    # K_L π0π0
  :D0toKlOmega,     # K_L ω
  :D0toKPiPi0,      # CP-mixed
  :D0toKPiPiPi      # CP-mixed
]

# =====================================================================
# Mode 1 : D0 → π+π−π0
# =====================================================================
alg_pipip0 = TagAnalysis.new("FplusPipip0")
alg_pipip0.set_header(["FplusPipip0Alg/FplusPipip0.h"])
          .set_constant({"ECMS" => [:double, 3.773]})
          .with_decay_card(decay_card_pipip0)

# Tag side: recoil anti-D0 in the ~12 CP / CP-mixed modes
alg_pipip0.tag_side(:D0) do |t|
  t.modes(*tag_modes)   # declared tag modes (trim-mode-lists default keeps only these)
  t.charm -1            # pin the tagged side to the anti-D0
end

# Signal side: everything the tag did not use → D0 → π+π−π0
alg_pipip0.signal_side do |s|
  s.photons 2                # π0 → γγ
  s.min_photon_angle 10.0    # γ angle to nearest charged track > 10°
  s.min_photon_energy 0.025  # E_γ > 25 MeV (barrel floor)
  s.charged(pip: 1, pim: 1)  # signal D0 → π+π−π0
  s.require_charge 0         # net charge zero
end

# Kinematic fit over the derived participants: 4C + π0 mass constraint
alg_pipip0.fit do |f|
  f.constrain_four_momentum                                              # 4C (nominal)
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0) # π0 from γγ
  f.chi2_cut 200
end

alg_pipip0
  .note(:detector_level_selection,
        "Charged tracks: |cosθ|<0.93, |Vz|<10 cm, Vr<1 cm, ≥1 positive and ≥1 negative track, net charge 0; " \
        "photons: TDC 0–700 ns, angle to any track > 10°, E>25 MeV (barrel) / 50 MeV (endcap). " \
        "Handled by the DTag track/shower selection; not expressible in the tag DSL.")
  .note(:pi0_reconstruction,
        "π0 from γγ: 1C Kalman fit with χ²<25; keep ≥1 π0 with M(γγ) in [0.115, 0.150] GeV, both photons in the barrel EMC.")
  .note(:background_veto,
        "K_S0 veto for D0→π+π−π0: reject any π+π− pair with M(π+π−) in [0.481, 0.514] GeV.")
  .note(:tag_side_selection,
        "Tag-side D0 candidates selected with MBC and ΔE (±3σ); MBC 1.86–1.87 GeV for K_L0 tag modes. " \
        "Stored unconditionally (no BOSS-level window) — DT yield extracted from a 2D MBC(tag)–MBC(sig) fit at ROOT level.")

alg_pipip0.apply
alg_pipip0.execute_on([data_3773, incMC_3773, exMC_pipip0])

# =====================================================================
# Mode 2 : D0 → K+K−π0
# =====================================================================
alg_KKp0 = TagAnalysis.new("FplusKKp0")
alg_KKp0.set_header(["FplusKKp0Alg/FplusKKp0.h"])
        .set_constant({"ECMS" => [:double, 3.773]})
        .with_decay_card(decay_card_KKp0)

# Tag side: recoil anti-D0 in the ~12 CP / CP-mixed modes
alg_KKp0.tag_side(:D0) do |t|
  t.modes(*tag_modes)
  t.charm -1
end

# Signal side: everything the tag did not use → D0 → K+K−π0
alg_KKp0.signal_side do |s|
  s.photons 2
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
  s.charged(kp: 1, km: 1)    # signal D0 → K+K−π0
  s.require_charge 0
end

# Kinematic fit over the derived participants: 4C + π0 mass constraint
alg_KKp0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_KKp0
  .note(:detector_level_selection,
        "Charged tracks: |cosθ|<0.93, |Vz|<10 cm, Vr<1 cm, ≥1 positive and ≥1 negative track, net charge 0; " \
        "photons: TDC 0–700 ns, angle to any track > 10°, E>25 MeV (barrel) / 50 MeV (endcap). " \
        "Handled by the DTag track/shower selection.")
  .note(:pi0_reconstruction,
        "π0 from γγ: 1C Kalman fit with χ²<25; keep ≥1 π0 with M(γγ) in [0.115, 0.150] GeV, both photons in the barrel EMC.")
  .note(:tag_side_selection,
        "Tag-side D0 candidates selected with MBC and ΔE (±3σ); MBC 1.86–1.87 GeV for K_L0 tag modes. " \
        "Stored unconditionally — DT yield extracted from a 2D MBC(tag)–MBC(sig) fit at ROOT level.")

alg_KKp0.apply
alg_KKp0.execute_on([data_3773, incMC_3773, exMC_KKp0])