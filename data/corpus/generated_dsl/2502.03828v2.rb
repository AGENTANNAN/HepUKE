# =====================================================================
# Semileptonic D -> K1(1270) mu nu at psi(3770)   (BOSS part)
# Single-tag analysis (TagAnalysis surface): the opposite D is tagged
# (D- against D+ signal, anti-D0 against D0 signal); the semileptonic D
# is the signal side, carrying a missing nu_mu.
# =====================================================================

### ------------------------- Dataset preparation ------------------------- ###
psi3770_data  = DatasetManager.real_data.find("712_3773")   # psi(3770) real data (@ 3.773 GeV)
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773") # corresponding inclusive MC

# Decay card for the D+ signal: psi(3770) -> D+ D-,
#   D+ -> K- pi+ pi0 mu+ nu_mu (signal), D- -> K+ pi- pi- (tag), pi0 -> gamma gamma
decay_card_Dplus = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.0000 K- pi+ pi0 mu+ nu_mu PHSP;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi- PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for the D0 signal: psi(3770) -> D0 anti-D0,
#   D0 -> K- pi+ pi- mu+ nu_mu (signal), anti-D0 -> K+ pi- (tag)
decay_card_D0 = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 K- pi+ pi- mu+ nu_mu PHSP;
    Enddecay

    Decay anti-D0
    1.0000 K+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# 1,000,000-event exclusive signal MC for each signal mode
exMC_Dplus = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_Dplus_K1mu"
  config.related_dataset = psi3770_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_Dplus
  config.cross_section   = :default
end

exMC_D0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0_K1mu"
  config.related_dataset = psi3770_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_D0
  config.cross_section   = :default
end

### ----------- D+ -> K1(1270)^0 mu+ nu_mu  (ST against D-) ----------- ###
alg_Dplus = TagAnalysis.new("DplusK1munu")
alg_Dplus.set_header(["DplusK1munuAlg/DplusK1munu.h"])
         .set_constant({"ECMS" => [:double, 3.773]})
         .with_decay_card(decay_card_Dplus)

# Tag side: the opposite D- reconstructed in six charged-D tag modes
alg_Dplus.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1                                      # pin the tagged side to D-
  t.window :mBC, min: 1.863, max: 1.877           # explicit m_BC tag window
end

# Signal side: K- pi+ pi0 mu+ nu_mu  (K1(1270)^0 -> K- pi+ pi0)
alg_Dplus.signal_side do |s|
  s.charged(km: 1, pip: 1, mup: 1)                # one K-, one pi+, one mu+
  s.photons 2                                     # two photons for the pi0
  s.missing :nu_mu                                # massless missing neutrino
  s.require_charge 1                              # net signal-side charge +1
  s.min_photon_angle 10.0                         # photon isolation from tracks
  s.min_photon_energy 0.025                       # standard barrel energy threshold
end

# Kinematic fit: 4C + pi0 mass constraint, chi2 < 200
alg_Dplus.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_Dplus
  .note(:tag_deltaE_windows, "per-mode DeltaE windows on the tag side (mode-dependent central values and widths) are applied on the stored tag DeltaE; the DSL stores tag DeltaE unconditionally and cannot express per-mode tag windows")
  .note(:pid_correction_method, "muon identification uses the probability method with the fixed CL_mu>0.001, CL_mu>CL_K, CL_mu>CL_pi requirements plus an additional EMC deposit < 0.28 GeV cut on the muon candidate")
  .note(:background_veto, "extra-shower vetoes suppress pi0/gamma backgrounds: no additional pi0, E_extra_gamma < 0.20 GeV, and cos_theta between the missing momentum and the extra photon < 0.69")
  .note(:signal_mass_windows, "signal mass window M(K pi pi0 mu) < 1.72 GeV/c^2 and K1(1270) mass window 1.163-1.343 GeV/c^2; peaking backgrounds suppressed by mass and PID cuts")
  .note(:generator, "signal MC should use a dedicated generator based on the D -> K pi pi e nu amplitude analysis rather than PHSP")

alg_Dplus.apply                                  # no Selection argument for a TagAnalysis
alg_Dplus.execute_on([psi3770_data, psi3770_incMC, exMC_Dplus])

### ----------- D0 -> K1(1270)^- mu+ nu_mu  (ST against anti-D0) ----------- ###
alg_D0 = TagAnalysis.new("D0K1munu")
alg_D0.set_header(["D0K1munuAlg/D0K1munu.h"])
      .set_constant({"ECMS" => [:double, 3.773]})
      .with_decay_card(decay_card_D0)

# Tag side: the opposite anti-D0 reconstructed in three neutral-D tag modes
alg_D0.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1                                      # pin the tagged side to anti-D0
  t.window :mBC, min: 1.859, max: 1.873           # explicit m_BC tag window
end

# Signal side: K- pi+ pi- mu+ nu_mu  (K1(1270)^- -> K- pi+ pi-)
alg_D0.signal_side do |s|
  s.charged(km: 1, pip: 1, pim: 1, mup: 1)        # one K-, one pi+, one pi-, one mu+
  s.missing :nu_mu                                # massless missing neutrino
  s.require_charge 0                              # net signal-side charge 0
  s.min_photon_angle 10.0
end

# Kinematic fit: 4C, chi2 < 200
alg_D0.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_D0
  .note(:tag_deltaE_windows, "per-mode DeltaE windows on the tag side (mode-dependent central values and widths) are applied on the stored tag DeltaE; the DSL stores tag DeltaE unconditionally and cannot express per-mode tag windows")
  .note(:pid_correction_method, "muon identification uses the probability method with the fixed CL_mu>0.001, CL_mu>CL_K, CL_mu>CL_pi requirements plus an additional EMC deposit < 0.28 GeV cut on the muon candidate")
  .note(:background_veto, "extra-shower vetoes suppress pi0/gamma backgrounds: no additional pi0, E_extra_gamma < 0.20 GeV, and cos_theta between the missing momentum and the extra photon < 0.54")
  .note(:signal_mass_windows, "signal mass window M(K pi pi mu) < 1.65 GeV/c^2 and K1(1270) mass window 1.163-1.343 GeV/c^2; peaking backgrounds suppressed by mass and PID cuts")
  .note(:generator, "signal MC should use a dedicated generator based on the D -> K pi pi e nu amplitude analysis rather than PHSP")

alg_D0.apply                                    # no Selection argument for a TagAnalysis
alg_D0.execute_on([psi3770_data, psi3770_incMC, exMC_D0])