# =============================================================================
# BOSS DSL — ψ(3770) (√s = 3.773 GeV) → D D̄  —  tag + signal (D → η + hadrons, η → γγ)
#
# Two tag analyses are run:
#   (1) D+-tag side: six hadronic tag modes  → D̄0 signal
#   (2) D0-tag side: three hadronic tag modes → D+  signal
#
# Only the modes actually available in the DTag machinery are declared
# (the remaining candidate modes are simply dropped).
# =============================================================================

### Dataset preparation ###
data_3773  = DatasetManager.real_data.find("712_3773")      # ψ(3770) real data at √s = 3.773 GeV
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")   # matching inclusive MC at 3.773 GeV

# Decay card for analysis 1: ψ(3770) → D+ D-, signal D → η + hadrons, η → γγ
decay_card_dp_tag = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- PHSP;
    Enddecay

    Decay D-
    1.0000 eta K- pi+ pi0 PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for analysis 2: ψ(3770) → D0 D̄0, signal D̄0 → η + hadrons, η → γγ
decay_card_d0_tag = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay anti-D0
    1.0000 eta K+ pi- pi0 PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Signal exclusive MC for the two tag analyses
exMC_dp_tag = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_DpTag_eta_hadrons"
  config.related_dataset = data_3773              # associated real dataset
  config.events          = 100000
  config.decay_card      = decay_card_dp_tag
  config.cross_section   = :default
end

exMC_d0_tag = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0Tag_eta_hadrons"
  config.related_dataset = data_3773
  config.events          = 100000
  config.decay_card      = decay_card_d0_tag
  config.cross_section   = :default
end

# =============================================================================
# Tag analysis 1 — D+ tag side (single tag) / D̄0 signal side
# =============================================================================
alg_dp = TagAnalysis.new("DpTagEtaHadrons")
alg_dp.set_header(["DpTagEtaHadronsAlg/DpTagEtaHadrons.h"])
      .set_constant({ "ECMS" => [:double, 3.773] })   # beam energy fixed to 1.8865 GeV per beam
      .with_decay_card(decay_card_dp_tag)

# Tag side: one tag_side call  ->  single tag (ST) of the D+
alg_dp.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi,        # D+ → K- π+ π+
          :DptoKsPi,         # D+ → K_S0 π+
          :DptoKPiPiPi0,     # D+ → K- π+ π+ π0
          :DptoKsPiPi0,      # D+ → K_S0 π+ π0
          :DptoKsPiPiPi,     # D+ → K_S0 π+ π+ π-
          :DptoKKPi          # D+ → K+ K- π+
  t.charm 1                  # pinned to the D+ tag side
end

# Signal side: what the tag did not use (D̄0 → η + hadrons, η → γγ)
alg_dp.signal_side do |s|
  s.photons 2                # exactly two signal-side photons
  s.min_photon_angle 10.0    # photon opening angle > 10°
  s.charged(at_least: true)  # at least one charged track on the signal side
end

# Kinematic fit: 4-momentum conservation + η mass constraint over the derived participants
alg_dp.fit do |f|
  f.constrain_four_momentum                                              # 4-momentum conservation (beam energy fixed)
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta) # 1C: M(γγ) = m(η)
  f.invariant_mass_of(:gamma, :gamma).between(0.520, 0.575)              # |M(γγ) - m(η)| < 27.5 MeV
  f.chi2_cut 200                                                        # χ² < 200
end

alg_dp
  .note(:ks0_reconstruction,
        "K_S0 → π+π- required with L/σ_L > 2 and M(π+π-) ∈ (0.486, 0.510) GeV; secondary-vertex " \
        "selection not expressible in the tag DSL")
  .note(:background_veto,
        "mass-window vetos against K_S0 (0.468-0.528 and 0.438-0.538 GeV), η (0.498-0.578 GeV), " \
        "ω (0.732-0.832 GeV), η' (0.908-1.008 GeV) and φ (0.990-1.390 GeV)")
  .note(:delta_e_signal,
        "per-mode ΔE_sig requirements applied on the signal side")
  .note(:peaking_background_veto,
        "peaking-background veto applied for the Kππ0η signal modes")
  .note(:tag_candidate_ranking,
        "best tag candidate is the one with the smallest |ΔE_tag|")
  .note(:signal_tag_angle,
        "opening angle between the signal D and the tag D required to exceed 160°")

alg_dp.apply                                   # tag spec takes NO Selection argument
root_files_dp = alg_dp.execute_on([data_3773, incMC_3773, exMC_dp_tag])

# =============================================================================
# Tag analysis 2 — D0 tag side (single tag) / D+ signal side
# =============================================================================
alg_d0 = TagAnalysis.new("D0TagEtaHadrons")
alg_d0.set_header(["D0TagEtaHadronsAlg/D0TagEtaHadrons.h"])
      .set_constant({ "ECMS" => [:double, 3.773] })   # beam energy fixed to 1.8865 GeV per beam
      .with_decay_card(decay_card_d0_tag)

# Tag side: one tag_side call  ->  single tag (ST) of the D0
alg_d0.tag_side(:D0) do |t|
  t.modes :D0toKPi,          # D0 → K- π+
          :D0toKPiPi0,       # D0 → K- π+ π0
          :D0toKPiPiPi       # D0 → K- π+ π+ π-
end

# Signal side: D+ → η + hadrons, η → γγ
alg_d0.signal_side do |s|
  s.photons 2                # exactly two signal-side photons
  s.min_photon_angle 10.0    # photon opening angle > 10°
  s.charged(at_least: true)  # at least one charged track on the signal side
end

# Kinematic fit — same hypothesis as analysis 1
alg_d0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.invariant_mass_of(:gamma, :gamma).between(0.520, 0.575)
  f.chi2_cut 200
end

alg_d0
  .note(:ks0_reconstruction,
        "K_S0 → π+π- required with L/σ_L > 2 and M(π+π-) ∈ (0.486, 0.510) GeV; secondary-vertex " \
        "selection not expressible in the tag DSL")
  .note(:background_veto,
        "mass-window vetos against K_S0 (0.468-0.528 and 0.438-0.538 GeV), η (0.498-0.578 GeV), " \
        "ω (0.732-0.832 GeV), η' (0.908-1.008 GeV) and φ (0.990-1.390 GeV)")
  .note(:delta_e_signal,
        "per-mode ΔE_sig requirements applied on the signal side")
  .note(:peaking_background_veto,
        "peaking-background veto applied for the Kππ0η signal modes")
  .note(:tag_candidate_ranking,
        "best tag candidate is the one with the smallest |ΔE_tag|")
  .note(:signal_tag_angle,
        "opening angle between the signal D and the tag D required to exceed 160°")

alg_d0.apply                                   # tag spec takes NO Selection argument
root_files_d0 = alg_d0.execute_on([data_3773, incMC_3773, exMC_d0_tag])