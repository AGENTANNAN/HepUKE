# Dataset preparation and event selection for ψ(4178) → Ds+ Ds−
# Double-tag method: Ds− hadronic tag on one side, semileptonic Ds+ → η(′) e+ νe on the other.

### Dataset description ###
data_4178 = DatasetManager.real_data.find("703_4180")       # 3.19 fb⁻¹ real data at 4.178 GeV
incmc_4178 = DatasetManager.inclusive_mc.find("703_4180")   # corresponding inclusive MC sample

# Decay card for the signal Ds+ → η e+ νe (η → γγ)
decay_card_eta = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Ds+ Ds- PHSP;
    Enddecay

    Decay Ds+
    1.0000 eta e+ nu_e PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for the signal Ds+ → η′ e+ νe (η′ → η π+ π−, η → γγ)
decay_card_etap = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Ds+ Ds- PHSP;
    Enddecay

    Decay Ds+
    1.0000 eta' e+ nu_e PHSP;
    Enddecay

    Decay eta'
    1.0000 eta pi+ pi- PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 500k-event exclusive MC for each signal decay mode
exMC_eta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4178_DsToEtaENu"
  config.related_dataset = data_4178
  config.events          = 500_000
  config.decay_card      = decay_card_eta
  config.cross_section   = :default
end

exMC_etap = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4178_DsToEtaPrimeENu"
  config.related_dataset = data_4178
  config.events          = 500_000
  config.decay_card      = decay_card_etap
  config.cross_section   = :default
end

### Event selection (BOSS) — Ds+ → η e+ νe, η → γγ ###
alg_eta = TagAnalysis.new("DsToEtaENu")
alg_eta.set_header(["DsToEtaENuAlg/DsToEtaENu.h"])
       .set_constant({"ECMS" => [:double, 4.178]})
       .with_decay_card(decay_card_eta)

# Tag side: hadronic Ds− from all hadronic tag modes, pinned to the negative-charge side
alg_eta.tag_side(:Ds) do |t|
  t.mode_group(:hadronic)                     # the 14 hadronic Ds− tag modes
  t.charm -1                                  # tag the Ds− (negative-charge species)
  t.window :mBC, min: 2.010, max: 2.073       # M_BC in (2.010, 2.073) GeV/c²
  t.window :deltaE, abs: 0.04                 # |ΔE| < 0.04 GeV
end

# Signal side: the tag's unused tracks/showers
alg_eta.signal_side do |s|
  s.photons 2                                 # two photons from η → γγ
  s.charged(ep: 1)                            # exactly one e+ and no extra charged tracks
  s.missing :nu_e                             # missing neutrino (semileptonic)
  s.min_photon_angle 10.0                     # photon at least 10° from any track
  s.min_photon_energy 0.025                   # photon energy > 25 MeV
end

# Kinematic fit: 4C + γγ mass constraint to the nominal η mass
alg_eta.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.invariant_mass_of(:gamma, :gamma).between(0.50, 0.57)   # M(γγ) in [0.50, 0.57] GeV/c²
  f.chi2_cut 200
end

# BOSS-side procedures not expressible in the tag DSL
alg_eta
  .note(:tag_candidate_ranking, "when several Ds− tag candidates pass, the candidate whose recoil mass M_rec is closest to the Ds*+ mass is chosen; the Ds*+ transition photon (or π0) is then selected by minimizing |ΔE| within (-0.04, 0.04) GeV")
  .note(:track_quality_cuts, "signal-side charged tracks must satisfy the usual quality cuts |cosθ| < 0.93, |Vz| < 10 cm and Vr < 1 cm; these are applied to the tracks returned by the tag before the e+ is used in the fit")
  .note(:extra_energy_cut, "E_extra < 0.3 GeV, computed from all neutral clusters not assigned to the signal/tag photons")
  .note(:bremsstrahlung_recovery, "Bremsstrahlung energy is partly recovered by adding EMC showers within 10° of the e+ to the electron four-momentum before the kinematic fit")

alg_eta.apply                                # TagAnalysis#apply takes no Selection argument
alg_eta.execute_on([data_4178, incmc_4178, exMC_eta])

### Event selection (BOSS) — Ds+ → η′ e+ νe, η′ → η π+ π−, η → γγ ###
alg_etap = TagAnalysis.new("DsToEtaPrimeENu")
alg_etap.set_header(["DsToEtaPrimeENuAlg/DsToEtaPrimeENu.h"])
        .set_constant({"ECMS" => [:double, 4.178]})
        .with_decay_card(decay_card_etap)

# Same tag side as the η mode
alg_etap.tag_side(:Ds) do |t|
  t.mode_group(:hadronic)
  t.charm -1
  t.window :mBC, min: 2.010, max: 2.073
  t.window :deltaE, abs: 0.04
end

# Signal side: additionally one π+ and one π− from η′ → η π+ π−
alg_etap.signal_side do |s|
  s.photons 2                                 # two photons from η → γγ
  s.charged(ep: 1, pip: 1, pim: 1)            # e+, π+, π− and no extra charged tracks
  s.missing :nu_e
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

# Kinematic fit: 4C + γγ → η + π+π−η → η′ constraints
alg_etap.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.invariant_mass_of(:pip, :pim, :gamma, :gamma).constrain_to_nominal_mass_of(:etap)
  f.invariant_mass_of(:ep, :pip, :pim, :gamma, :gamma).between(0.0, 1.9)  # M(η′e+) < 1.9 GeV/c²
  f.chi2_cut 200
end

alg_etap
  .note(:tag_candidate_ranking, "when several Ds− tag candidates pass, the candidate whose recoil mass M_rec is closest to the Ds*+ mass is chosen; the Ds*+ transition photon (or π0) is then selected by minimizing |ΔE| within (-0.04, 0.04) GeV")
  .note(:track_quality_cuts, "signal-side charged tracks must satisfy the usual quality cuts |cosθ| < 0.93, |Vz| < 10 cm and Vr < 1 cm; these are applied to the tracks returned by the tag before the e+/π± are used in the fit")
  .note(:extra_energy_cut, "E_extra < 0.3 GeV, computed from all neutral clusters not assigned to the signal/tag photons")
  .note(:bremsstrahlung_recovery, "Bremsstrahlung energy is partly recovered by adding EMC showers within 10° of the e+ to the electron four-momentum before the kinematic fit")
  .note(:helicity_cut, "for the η′ → γ ρ0 sub-mode a helicity cut |cosθ_hel| < 0.85 is applied")

alg_etap.apply
alg_etap.execute_on([data_4178, incmc_4178, exMC_etap])