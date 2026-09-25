# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# ψ(3770) at √s = 3.773 GeV — 2.92 fb^-1 of real data and the matching inclusive MC
data_3773 = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# ---------------------------------------------------------------------------
# Decay cards — ψ(3770) → D0 D0bar, one D in a CP-eigenstate tag mode and the
# opposite side in the flavour mode (K+ π-). One card per CP tag mode.
# ---------------------------------------------------------------------------
decay_card_KK = <<~DECAYCARD
    Decay psi(3770)
    1.0 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0 K+ K- PHSP;
    Enddecay

    Decay anti-D0
    1.0 K+ pi- PHSP;
    Enddecay

    End
DECAYCARD

decay_card_PiPi = <<~DECAYCARD
    Decay psi(3770)
    1.0 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0 pi+ pi- PHSP;
    Enddecay

    Decay anti-D0
    1.0 K+ pi- PHSP;
    Enddecay

    End
DECAYCARD

decay_card_KsPi0Pi0 = <<~DECAYCARD
    Decay psi(3770)
    1.0 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0 K_S0 pi0 pi0 PHSP;
    Enddecay

    Decay anti-D0
    1.0 K+ pi- PHSP;
    Enddecay

    Decay K_S0
    1.0 pi+ pi- PHSP;
    Enddecay

    Decay pi0
    1.0 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

decay_card_Pi0Pi0 = <<~DECAYCARD
    Decay psi(3770)
    1.0 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0 pi0 pi0 PHSP;
    Enddecay

    Decay anti-D0
    1.0 K+ pi- PHSP;
    Enddecay

    Decay pi0
    1.0 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

decay_card_Rho0Pi0 = <<~DECAYCARD
    Decay psi(3770)
    1.0 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0 rho0 pi0 PHSP;
    Enddecay

    Decay anti-D0
    1.0 K+ pi- PHSP;
    Enddecay

    Decay rho0
    1.0 pi+ pi- PHSP;
    Enddecay

    Decay pi0
    1.0 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

decay_card_KsPi0 = <<~DECAYCARD
    Decay psi(3770)
    1.0 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0 K_S0 pi0 PHSP;
    Enddecay

    Decay anti-D0
    1.0 K+ pi- PHSP;
    Enddecay

    Decay K_S0
    1.0 pi+ pi- PHSP;
    Enddecay

    Decay pi0
    1.0 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

decay_card_KsEta = <<~DECAYCARD
    Decay psi(3770)
    1.0 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0 K_S0 eta PHSP;
    Enddecay

    Decay anti-D0
    1.0 K+ pi- PHSP;
    Enddecay

    Decay K_S0
    1.0 pi+ pi- PHSP;
    Enddecay

    Decay eta
    1.0 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

decay_card_KsOmega = <<~DECAYCARD
    Decay psi(3770)
    1.0 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0 K_S0 omega PHSP;
    Enddecay

    Decay anti-D0
    1.0 K+ pi- PHSP;
    Enddecay

    Decay K_S0
    1.0 pi+ pi- PHSP;
    Enddecay

    Decay omega
    1.0 pi+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.0 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Kπ–Kπ sample for right-sign / wrong-sign studies
decay_card_KPiKPi = <<~DECAYCARD
    Decay psi(3770)
    1.0 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0 K- pi+ PHSP;
    Enddecay

    Decay anti-D0
    1.0 K+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# ---------------------------------------------------------------------------
# Exclusive MC samples — 200k events per CP-tag mode plus the Kπ–Kπ sample
# ---------------------------------------------------------------------------
exMC_KK = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3773_CP_KK"
  config.related_dataset = data_3773
  config.events          = 200000
  config.decay_card      = decay_card_KK
  config.cross_section   = :default
end

exMC_PiPi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3773_CP_PiPi"
  config.related_dataset = data_3773
  config.events          = 200000
  config.decay_card      = decay_card_PiPi
  config.cross_section   = :default
end

exMC_KsPi0Pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3773_CP_KsPi0Pi0"
  config.related_dataset = data_3773
  config.events          = 200000
  config.decay_card      = decay_card_KsPi0Pi0
  config.cross_section   = :default
end

exMC_Pi0Pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3773_CP_Pi0Pi0"
  config.related_dataset = data_3773
  config.events          = 200000
  config.decay_card      = decay_card_Pi0Pi0
  config.cross_section   = :default
end

exMC_Rho0Pi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3773_CP_Rho0Pi0"
  config.related_dataset = data_3773
  config.events          = 200000
  config.decay_card      = decay_card_Rho0Pi0
  config.cross_section   = :default
end

exMC_KsPi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3773_CP_KsPi0"
  config.related_dataset = data_3773
  config.events          = 200000
  config.decay_card      = decay_card_KsPi0
  config.cross_section   = :default
end

exMC_KsEta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3773_CP_KsEta"
  config.related_dataset = data_3773
  config.events          = 200000
  config.decay_card      = decay_card_KsEta
  config.cross_section   = :default
end

exMC_KsOmega = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3773_CP_KsOmega"
  config.related_dataset = data_3773
  config.events          = 200000
  config.decay_card      = decay_card_KsOmega
  config.cross_section   = :default
end

exMC_KPiKPi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3773_KPiKPi"
  config.related_dataset = data_3773
  config.events          = 200000
  config.decay_card      = decay_card_KPiKPi
  config.cross_section   = :default
end

### Tag-based event selection (BOSS) ###
# Double tag: one side is a CP-eigenstate single tag, the opposite side is the
# flavour tag D0 → K∓π±. mBC / ΔE / tag-mass observables are stored
# unconditionally; the mode-dependent ΔE windows are applied later in ROOT.
alg = TagAnalysis.new("D0CPTagKPi")
alg.set_header(["D0CPTagKPiAlg/D0CPTagKPi.h"])
   .set_constant({"ECMS" => [:double, 3.773]})

# Tag side 1 — the CP-eigenstate single tag (charm of the tagged side not pinned:
# D0 and D0bar both enter, as the CP eigenstates are charge-blind)
alg.tag_side(:D0) do |t|
  t.modes :D0toKK, :D0toPiPi, :D0toKsPi0Pi0, :D0toPi0Pi0,
          :D0toRho0Pi0, :D0toKsPi0, :D0toKsEta, :D0toKsOmega
end

# Tag side 2 — the flavour tag D0 → Kπ on the opposite side
alg.tag_side(:D0) do |t|
  t.modes :D0toKPi
end

# Kinematic fit over the two tag sides: total four-momentum conservation plus a
# constraint of both tag masses to the D0 nominal mass
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

# BOSS-side procedures that have no dedicated DSL construct
alg.note(:delta_e_selection,
         "mode-dependent ΔE windows applied per tag mode (K+K- ±0.025 GeV; " \
         "π+π- and K-π+ ±0.030 GeV; K_S0π0π0 [-0.080, +0.045] GeV); only the " \
         "candidate combination with the least |ΔE| is retained per mode, and " \
         "the windows are applied on the stored ΔE in the ROOT stage")
   .note(:background_veto,
         "in the K+K- and π+π- tag modes, cosmic-ray and Bhabha backgrounds are " \
         "suppressed by requiring a TOF time difference < 5 ns, inconsistency " \
         "with the μ+μ-/e+e- hypothesis, and rejection of events carrying an " \
         "extra EMC shower > 50 MeV or an additional MDC track")
   .note(:intermediate_reconstruction,
         "π0/η candidates built from photon pairs requiring at least one barrel " \
         "photon, mass windows 0.115-0.150 (0.505-0.570) GeV/c², followed by a 1C " \
         "mass constraint; K_S0 candidates built from π+π- via a vertex fit " \
         "(χ² < 100) on loose IP tracks (|Vz| < 20 cm, no transverse requirement) " \
         "then a second fit constraining the momentum to the IP with flight " \
         "significance L/σ_L > 2 and 0.487 < m(π+π-) < 0.511 GeV/c²")
   .note(:rho_omega_mass_windows,
         "the ρ0π0 and K_S0ω tag modes additionally require " \
         "0.60 < m(π+π-) < 0.95 GeV/c² and 0.72 < m(π+π-π0) < 0.84 GeV/c²")

alg.with_decay_card(decay_card_KK).apply

# Execute on the real data, inclusive MC and all signal exclusive MC samples
root_files = alg.execute_on([
  data_3773, incMC_3773,
  exMC_KK, exMC_PiPi, exMC_KsPi0Pi0, exMC_Pi0Pi0,
  exMC_Rho0Pi0, exMC_KsPi0, exMC_KsEta, exMC_KsOmega,
  exMC_KPiKPi
])