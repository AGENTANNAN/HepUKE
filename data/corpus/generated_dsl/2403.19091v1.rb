# frozen_string_literal: true
#
# ψ(3770) semileptonic double-tag analysis (BOSS part):
#   D0 → K_S0 π− π0 e+ νe   (tag side: hadronic anti-D0)
#   D+ → K_S0 π+ π− e+ νe   (tag side: hadronic D−)
# Tag-based analysis → TagAnalysis (no Selection / select_track / pid).

### Dataset description ###
data_3770  = DatasetManager.real_data.find("712_3773")       # ψ(3770) real data @ 3.773 GeV (2.93 fb⁻¹)
incMC_3770 = DatasetManager.inclusive_mc.find("712_3773")    # matching inclusive MC sample

# ---------------------------------------------------------------------------
# Decay cards (signal MC, EvtGen syntax)
# ---------------------------------------------------------------------------
# ψ(3770) → D0 anti-D0, signal D0 → K_S0 π− π0 e+ νe (ISGW2 form factor)
decay_card_D0 = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 K_S0 pi- pi0 e+ nu_e ISGW2;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# ψ(3770) → D+ D−, signal D+ → K_S0 π+ π− e+ νe (ISGW2 form factor)
decay_card_Dplus = <<~DECAYCARD
    Decay psi(3770)
    1.000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.000 K_S0 pi+ pi- e+ nu_e ISGW2;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# ---------------------------------------------------------------------------
# Exclusive signal MC (200k events per signal mode)
# ---------------------------------------------------------------------------
exMC_D0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0_KSpipi0e_nu"
  config.related_dataset = data_3770
  config.events          = 200_000
  config.decay_card      = decay_card_D0
  config.cross_section   = :default
end
exMC_D0.save_to_config(format: :yaml, file_path: 'temp_for_test')

exMC_Dplus = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_Dplus_KSpipie_nu"
  config.related_dataset = data_3770
  config.events          = 200_000
  config.decay_card      = decay_card_Dplus
  config.cross_section   = :default
end
exMC_Dplus.save_to_config(format: :yaml, file_path: 'temp_for_test')

# ===========================================================================
# Channel 1 : tag anti-D0 (hadronic) + signal D0 → K_S0 π− π0 e+ νe
# ===========================================================================
alg_D0 = TagAnalysis.new("D0TagKSpipi0Enu")
alg_D0.set_header(["D0TagKSpipi0EnuAlg/D0TagKSpipi0Enu.h"])
      .set_constant({"ECMS" => [:double, 3.773]})
      .with_decay_card(decay_card_D0)

# Tag side: anti-D0 reconstructed from the hadronic tag modes (charm = -1 pins anti-D0)
alg_D0.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi, :D0toKPiPiPiPi0
  t.charm(-1)
end

# Signal side: what the tag did not use — K_S0(→π+π−), π−, π0(→γγ) and e+ with a missing νe
alg_D0.signal_side do |s|
  s.photons 2                        # two photons from π0 → γγ
  s.charged(pip: 1, pim: 2, ep: 1)   # π+π− from K_S0, plus the π− and the e+ of the D0
  s.require_charge 0                 # signal-side net charge = 0
  s.missing :nu_e                    # massless missing neutrino (semileptonic form)
  s.min_photon_energy 0.025          # Eγ > 25 (50) MeV in barrel (endcap)
end

# 4C kinematic fit with the D mass constraint (tag D0) and the π0 mass constraint
alg_D0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # π0 mass-constrained fit
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)            # tagged D mass constraint
  f.chi2_cut 200
end

# BOSS-side procedures without a formal DSL construct
alg_D0
  .note(:ks0_reconstruction,
        "K_S0 reconstructed from π+π− via a secondary-vertex fit, requiring a flight "
        "distance significance > 2σ and a π+π− invariant mass in (0.486, 0.510) GeV/c²")
  .note(:pi0_reconstruction,
        "π0 reconstructed from γγ with photon energies > 25 (50) MeV in the barrel "
        "(endcap) and a γγ invariant mass in (0.115, 0.150) GeV/c², followed by a π0 "
        "mass-constrained fit")
  .note(:pid_correction_method,
        "positron PID: L'(e) > 0.001, L'(e)/(L'(e)+L'(π)+L'(K)) > 0.8 and "
        "E/p > 0.18·χ²(dE/dx) + 0.32")
  .note(:bremsstrahlung_recovery,
        "photons within 5° of the e+ with E > 50 MeV are recovered into the e+ four-momentum")
  .note(:background_veto,
        "D0 signal-side suppression: E(π0) > 0.22 GeV, |cosθ(π0)| < 0.83, and "
        "M(K_S0 π− π0 π+) with the e+ assigned the pion mass < 1.78 GeV/c²")

alg_D0.apply
alg_D0.execute_on([data_3770, incMC_3770, exMC_D0])

# ===========================================================================
# Channel 2 : tag D− (hadronic) + signal D+ → K_S0 π+ π− e+ νe
# ===========================================================================
alg_Dplus = TagAnalysis.new("DpTagKSpipieNu")
alg_Dplus.set_header(["DpTagKSpipieNuAlg/DpTagKSpipieNu.h"])
        .set_constant({"ECMS" => [:double, 3.773]})
        .with_decay_card(decay_card_Dplus)

# Tag side: D− reconstructed from the hadronic tag modes (charm = -1 pins D−)
alg_Dplus.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm(-1)
end

# Signal side: what the tag did not use — K_S0(→π+π−), two π, e+, and a missing νe
alg_Dplus.signal_side do |s|
  s.charged(pip: 2, pim: 2, ep: 1)   # π+π− from K_S0, plus π+π− and the e+ of the D+
  s.missing :nu_e                    # massless missing neutrino (semileptonic form)
end

# 4C kinematic fit with the (tag) D mass constraint
alg_Dplus.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:Dplus)   # tagged D mass constraint
  f.chi2_cut 200
end

# BOSS-side procedures without a formal DSL construct
alg_Dplus
  .note(:ks0_reconstruction,
        "K_S0 reconstructed from π+π− via a secondary-vertex fit, requiring a flight "
        "distance significance > 2σ and a π+π− invariant mass in (0.486, 0.510) GeV/c²")
  .note(:pid_correction_method,
        "positron PID: L'(e) > 0.001, L'(e)/(L'(e)+L'(π)+L'(K)) > 0.8 and "
        "E/p > 0.18·χ²(dE/dx) + 0.32")
  .note(:bremsstrahlung_recovery,
        "photons within 5° of the e+ with E > 50 MeV are recovered into the e+ four-momentum")
  .note(:background_veto,
        "D+ signal-side suppression: M(K_S0 π+ π− π+) with the e+ assigned the pion mass "
        "< 1.83 GeV/c², cosθ(e+, π−) < 0.95, M(K_S0 π+ π− π+ π0) with the e+ assigned the "
        "pion mass < 1.4 GeV/c², and cosθ(missing momentum, most energetic unused shower) < 0.88")

alg_Dplus.apply
alg_Dplus.execute_on([data_3770, incMC_3770, exMC_Dplus])