# =====================================================================
# BESIII: Search for dark baryon in Ξ- → π- + invisible
#   [arXiv:2505.22140v1]
#
# Using (1.0087 ± 0.0044) × 10^10 J/ψ events. ST-DT hyperon tagging
# (NOT D-tag — Ξ hyperon reconstruction with secondary vertex fits).
# ST: Ξbar+ → π+ Λbar (Λbar → pbar π+)
# DT: signal Ξ- → π- + invisible dark baryon χ.
# 5 m_χ hypotheses: 1.07, 1.10, m_Λ, 1.13, 1.16 GeV/c².
# Ordinary analysis (Algorithm + Selection).
# =====================================================================

jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Single decay card covering the full chain: J/ψ → Ξ- Ξbar+,
# Ξbar+ → π+ Λbar → π+ pbar π+, Ξ- → π- + χ (invisible).
# Since χ is invisible/dark, it is produced as a "stable" invisible particle.
decay_card_signal = <<~DECAYCARD
  Decay J/psi
  1.0  Xi- anti-Xi-  PHSP;
  Enddecay

  Decay anti-Xi-
  1.0  anti-Lambda pi+  PHSP;
  Enddecay

  Decay anti-Lambda
  1.0  anti-p- pi+  PHSP;
  Enddecay

  Decay Xi-
  1.0  pi- chi0  PHSP;
  Enddecay

  End
DECAYCARD

# For each m_χ hypothesis we need a separate MC (mass varies)
# but create_exclusive_mc only supports a single mass card per call.
# The 5 mass hypotheses are handled with the mass as a note for ROOT.
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "XiMinus2PiChibInv"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

alg_name = "XiMinus2PiInvisible"
algorithm = Algorithm.new(alg_name)
algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 3.097]})

event_selection = Selection.new
event_selection.select_track {
                  cos_theta 0.93
                  Vz        10.0
                  Vr        1.0
                  nChrp     ">=2"     # π+(Λbar) + π+(ST tag)
                  nChrn     ">=3"     # pbar + π-(Λbar) + π-(DT signal)
                  nNet      ">=0"
                }
               .select_photon {
                  tdc_emc_start     0
                  tdc_emc_end       14
                  energyThreshold_b 0.025
                  energyThreshold_e 0.050
                }
               .pid(method: :probability) {
                  prob_cut 0.0
                  identify :prm, against: [:pim, :km]     # anti-proton for Λbar
                  identify :pim, against: [:prm, :km]     # π- for signal and Λbar
               }
               .secondary_vertex_fit(:anti_Lambda, daughters: [:prm, :pip]) {
                  chi2_cut 200
               }
               .secondary_vertex_fit(:"anti-Xi-", daughters: [:anti_Lambda, :pip]) {
                  chi2_cut 200
               }

# The DT side: exactly one additional π- + missing χ.
# 1C kinematic fit constraining missing mass to m_χ.
# T2: competing hypothesis veto against Ξ- → π- Λ.
event_selection.partial_rec do |pr|
  pr.missing :chi0
  pr.charged(pim: 1)
end

event_selection.kinematic_fit([:anti_Lambda, :pip, :pim]) do |fit|
  fit.constrain_four_momentum
  fit.invariant_mass_of(:prm, :pip).constrain_to_nominal_mass_of(:anti_Lambda)
  fit.invariant_mass_of(:anti_Lambda, :pip).constrain_to_nominal_mass_of(:"anti-Xi-")
  fit.chi2_cut 20
  fit.nominal
end

algorithm
  .note(:st_selection,
        "ST Ξbar+ reconstructed via vertex fits on pbar π+ (anti-Lambda) and anti-Lambda π+ (Ξbar+). " \
        "Best candidate by minimum sum of |M(pbar π+) - M_Λ| + |M(anti-Lambda π+) - M_Ξ+|. " \
        "|M(pbar π+) - M_Λ| < 0.004 GeV/c². M_recoil(Ξbar+) ∈ (1.290, 1.345) GeV/c² for Ξ- signal.")
  .note(:dt_selection,
        "DT signal: exactly one extra π- in remaining tracks. " \
        "1C kinematic fit: J/ψ → pbar π+ π+ π- + invisible, constraining M(invisible) = m_χ. " \
        "χ²_π-χ < 20. Competing hypothesis veto: χ²_π-Λ > χ²_π-χ (except m_χ = m_Λ case). " \
        "P_π- momentum windows optimized by Punzi significance, range (0.070-0.192) GeV/c " \
        "depending on m_χ hypothesis.")
  .note(:emc_selection,
        "E_EMC used to separate signal from Ξ- → π- Λ(nπ0) background. " \
        "Showers: barrel |cosθ| < 0.80 or endcap 0.86 < |cosθ| < 0.92. " \
        "Isolation angle: >10° from charged tracks, >20° from anti-protons. " \
        "EMC timing [0, 700] ns. Data-driven E_EMC correction for anti-proton interactions.")
  .note(:mass_hypotheses,
        "Five dark baryon mass hypotheses: m_χ = 1.07, 1.10, m_Λ (≈1.116), 1.13, 1.16 GeV/c². " \
        "No significant signal observed. 90% CL ULs: (4.2-65) × 10⁻⁵. " \
        "The kinematic fit miss constraint mass m_χ is varied per hypothesis.")
  .note(:control_sample,
        "Control sample J/ψ → Ξ- Ξbar+ → Λ(pπ-) π- (6C fit) used for E_EMC data-driven correction. " \
        "E_EMC decomposed into E_EMC^π0 (from π0 of Λ→nπ0) and E_EMC^other (from charged/neutron noise).")
  .note(:emi_emc_datadriven,
        "E_EMC^other corrected using J/ψ → Ξ-(π-Λ(nπ0)) Ξbar+(π+ Λbar(pbar π+)) control sample. " \
        "4C kinematic fit with neutron as missing track. " \
        "Sampled per (anti-proton momentum, θ) bin. This correction is performed externally.")

algorithm.with_decay_card(decay_card_signal).apply(event_selection)
algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])