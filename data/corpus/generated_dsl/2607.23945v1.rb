# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
data_3773  = DatasetManager.real_data.find("712_3773")     # ψ(3770) real data (~20.3 fb⁻¹)
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")  # Corresponding inclusive MC sample

# ----- Decay cards (EvtGen format) -----
# Signal side is the semileptonic D (PHOTOS/ISGW2); tag side is the hadronic anti-D.
decay_card_d0_epi = <<~DECAYCARD
    Decay psi(3770)
    1.0 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0 pi- e+ nu_e PHOTOS ISGW2;
    Enddecay

    Decay anti-D0
    1.0 K+ pi- PHSP;
    Enddecay

    End
DECAYCARD

decay_card_d0_mupi = <<~DECAYCARD
    Decay psi(3770)
    1.0 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0 pi- mu+ nu_mu PHOTOS ISGW2;
    Enddecay

    Decay anti-D0
    1.0 K+ pi- PHSP;
    Enddecay

    End
DECAYCARD

decay_card_dp_epi = <<~DECAYCARD
    Decay psi(3770)
    1.0 D+ D- PHSP;
    Enddecay

    Decay D+
    1.0 pi0 e+ nu_e PHOTOS ISGW2;
    Enddecay

    Decay pi0
    1.0 gamma gamma PHSP;
    Enddecay

    Decay D-
    1.0 K+ pi- pi- PHSP;
    Enddecay

    End
DECAYCARD

decay_card_dp_mupi = <<~DECAYCARD
    Decay psi(3770)
    1.0 D+ D- PHSP;
    Enddecay

    Decay D+
    1.0 pi0 mu+ nu_mu PHOTOS ISGW2;
    Enddecay

    Decay pi0
    1.0 gamma gamma PHSP;
    Enddecay

    Decay D-
    1.0 K+ pi- pi- PHSP;
    Enddecay

    End
DECAYCARD

# ----- Exclusive MC samples (500k events per mode) -----
exMC_d0_epi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_d0_pienu"
  config.related_dataset = data_3773
  config.events          = 500_000
  config.decay_card      = decay_card_d0_epi
  config.cross_section   = :default
end

exMC_d0_mupi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_d0_pimunu"
  config.related_dataset = data_3773
  config.events          = 500_000
  config.decay_card      = decay_card_d0_mupi
  config.cross_section   = :default
end

exMC_dp_epi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_dp_pi0enu"
  config.related_dataset = data_3773
  config.events          = 500_000
  config.decay_card      = decay_card_dp_epi
  config.cross_section   = :default
end

exMC_dp_mupi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_dp_pi0munu"
  config.related_dataset = data_3773
  config.events          = 500_000
  config.decay_card      = decay_card_dp_mupi
  config.cross_section   = :default
end

### Event selection (BOSS, tag-based) ###
# One TagAnalysis per signal mode; one tag_side (hadronic anti-D) + signal_side with a missing neutrino.

# ---------- Mode 1: D0 -> pi- e+ nu_e ----------
alg_d0_epi = TagAnalysis.new("D0ToPiENu")
alg_d0_epi.set_header(["D0ToPiENuAlg/D0ToPiENu.h"])
          .set_constant({ "ECMS" => [:double, 3.773] })
          .with_decay_card(decay_card_d0_epi)

# Tag side: charm = -1 anti-D0, reconstructed through the listed hadronic modes
alg_d0_epi.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi, :D0toKPiPi0Pi0, :D0toKPiPiPiPi0
  t.charm -1
end

# Signal side: one pi- + one e+ (net charge 0) recoiling against a missing neutrino
alg_d0_epi.signal_side do |s|
  s.charged(pim: 1, ep: 1)
  s.require_charge 0
  s.missing :nu_e
end

# Kinematic fit: 4-momentum constraint, chi2 < 200
alg_d0_epi.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_d0_epi.apply
alg_d0_epi.execute_on([data_3773, incMC_3773, exMC_d0_epi])

# ---------- Mode 2: D0 -> pi- mu+ nu_mu ----------
alg_d0_mupi = TagAnalysis.new("D0ToPiMuNu")
alg_d0_mupi.set_header(["D0ToPiMuNuAlg/D0ToPiMuNu.h"])
           .set_constant({ "ECMS" => [:double, 3.773] })
           .with_decay_card(decay_card_d0_mupi)

alg_d0_mupi.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi, :D0toKPiPi0Pi0, :D0toKPiPiPiPi0
  t.charm -1
end

alg_d0_mupi.signal_side do |s|
  s.charged(pim: 1, mup: 1)
  s.require_charge 0
  s.missing :nu_mu
end

alg_d0_mupi.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_d0_mupi.apply
alg_d0_mupi.execute_on([data_3773, incMC_3773, exMC_d0_mupi])

# ---------- Mode 3: D+ -> pi0 e+ nu_e (pi0 -> gamma gamma) ----------
alg_dp_epi = TagAnalysis.new("DpToPi0ENu")
alg_dp_epi.set_header(["DpToPi0ENuAlg/DpToPi0ENu.h"])
          .set_constant({ "ECMS" => [:double, 3.773] })
          .with_decay_card(decay_card_dp_epi)

# Tag side: charm = -1 D-, reconstructed through the listed hadronic modes
alg_dp_epi.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
end

# Signal side: one e+ + two photons (from pi0), net charge +1, missing neutrino
alg_dp_epi.signal_side do |s|
  s.photons 2
  s.charged(ep: 1)
  s.require_charge 1
  s.min_photon_angle 10.0
  s.missing :nu_e
end

# Kinematic fit: 4-momentum + gamma-gamma invariant mass to nominal pi0, chi2 < 200
alg_dp_epi.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_dp_epi.apply
alg_dp_epi.execute_on([data_3773, incMC_3773, exMC_dp_epi])

# ---------- Mode 4: D+ -> pi0 mu+ nu_mu (pi0 -> gamma gamma) ----------
alg_dp_mupi = TagAnalysis.new("DpToPi0MuNu")
alg_dp_mupi.set_header(["DpToPi0MuNuAlg/DpToPi0MuNu.h"])
           .set_constant({ "ECMS" => [:double, 3.773] })
           .with_decay_card(decay_card_dp_mupi)

alg_dp_mupi.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
end

alg_dp_mupi.signal_side do |s|
  s.photons 2
  s.charged(mup: 1)
  s.require_charge 1
  s.min_photon_angle 10.0
  s.missing :nu_mu
end

alg_dp_mupi.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_dp_mupi.apply
alg_dp_mupi.execute_on([data_3773, incMC_3773, exMC_dp_mupi])