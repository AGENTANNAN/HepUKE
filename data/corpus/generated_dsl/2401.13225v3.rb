# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset preparation ###
data_3773  = DatasetManager.real_data.find("712_3773")      # ψ(3770) real data at √s = 3.773 GeV
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")   # corresponding inclusive MC sample

# Decay card for the muon channel: ψ(3770) → D⁻D⁺, D⁻ → K⁺π⁻π⁻ (tag), D⁺ → π⁺π⁻μ⁺ν_μ (signal, PHOTOS/ISGW2)
decay_card_signal_mu = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.0000 pi+ pi- mu+ nu_mu PHOTOS ISGW2;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi- PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for the electron channel: ψ(3770) → D⁻D⁺, D⁻ → K⁺π⁻π⁻ (tag), D⁺ → π⁺π⁻e⁺ν_e (signal, PHOTOS/ISGW2)
decay_card_signal_e = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.0000 pi+ pi- e+ nu_e PHOTOS ISGW2;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi- PHSP;
    Enddecay

    End
DECAYCARD

# 500k-event exclusive MC for the muon channel
exMC_signal_mu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3773_DpToPipPimMuNu"
  config.related_dataset = data_3773
  config.events          = 500000
  config.decay_card      = decay_card_signal_mu
  config.cross_section   = :default
end

# 500k-event exclusive MC for the electron channel
exMC_signal_e = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3773_DpToPipPimENu"
  config.related_dataset = data_3773
  config.events          = 500000
  config.decay_card      = decay_card_signal_e
  config.cross_section   = :default
end

### Event selection (BOSS) — single-tag semileptonic analysis ###
# The event is tagged by a fully reconstructed D⁻ (six hadronic tag modes); the recoil D⁺
# is searched for in π⁺π⁻ℓ⁺ν (ℓ = μ or e) with the neutrino left missing.
# Two independent signal channels → two TagAnalysis objects (Rule T1).

# ------------------------------------------------------------------ #
# Channel I: D⁺ → π⁺π⁻μ⁺ν_μ                                           #
# ------------------------------------------------------------------ #
alg_mu_name = "DpToPipPimMuNuST"
alg_mu = TagAnalysis.new(alg_mu_name)
alg_mu.set_header(["#{alg_mu_name}Alg/#{alg_mu_name}.h"])
      .set_constant({"ECMS" => [:double, 3.773]})   # √s = 3.773 GeV
      .with_decay_card(decay_card_signal_mu)

# Tag side: D⁻ fully reconstructed from the six hadronic modes (single tag)
alg_mu.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi,      # D⁻ → K⁺π⁻π⁻
          :DptoKPiPiPi0,   # D⁻ → K⁺π⁻π⁻π⁰
          :DptoKsPi,       # D⁻ → K_S π⁻
          :DptoKsPiPi0,    # D⁻ → K_S π⁻π⁰
          :DptoKsPiPiPi,   # D⁻ → K_S π⁻π⁺π⁻
          :DptoKKPi        # D⁻ → K⁺K⁻π⁻
  t.charm -1               # pin the tagged side to the D⁻
end

# Signal side: exactly one π⁺, one π⁻ and one same-flavour μ⁺; net event charge 0 is
# guaranteed by charm(-1) tag (D⁻) + charge(+1) signal multiset; neutrino is missing.
alg_mu.signal_side do |s|
  s.charged(pip: 1, pim: 1, mup: 1)
  s.missing :nu_mu          # massless neutrino → 2-arg AddMissTrack p4 treatment
end

# Four-momentum-constrained kinematic fit with the massless neutrino
alg_mu.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200            # retain events with χ² < 200
end

# Non-standard muon identification criterion on the signal side (not expressible with the
# fixed v1 lepton-PID recipe of tag analyses); efficiency handled outside the kinematic fit
alg_mu.note(:pid_correction_method, "signal-side muon identified by CL_mu > 0.001 AND CL_mu > CL_e, CL_pi, CL_K AND an EMC energy deposit within the 0.09-0.31 GeV MIP window; the fixed muon-PID recipe of the tag-analysis signal side is retuned to these criteria in the generated code")

alg_mu.apply                       # takes no Selection argument
root_files_mu = alg_mu.execute_on([data_3773, incMC_3773, exMC_signal_mu])

# ------------------------------------------------------------------ #
# Channel II: D⁺ → π⁺π⁻e⁺ν_e                                          #
# ------------------------------------------------------------------ #
alg_e_name = "DpToPipPimENuST"
alg_e = TagAnalysis.new(alg_e_name)
alg_e.set_header(["#{alg_e_name}Alg/#{alg_e_name}.h"])
     .set_constant({"ECMS" => [:double, 3.773]})
     .with_decay_card(decay_card_signal_e)

# Same tag side as the muon channel
alg_e.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi,      # D⁻ → K⁺π⁻π⁻
          :DptoKPiPiPi0,   # D⁻ → K⁺π⁻π⁻π⁰
          :DptoKsPi,       # D⁻ → K_S π⁻
          :DptoKsPiPi0,    # D⁻ → K_S π⁻π⁰
          :DptoKsPiPiPi,   # D⁻ → K_S π⁻π⁺π⁻
          :DptoKKPi        # D⁻ → K⁺K⁻π⁻
  t.charm -1               # tag the D⁻
end

# Signal side: exactly one π⁺, one π⁻ and one e⁺; neutrino missing
alg_e.signal_side do |s|
  s.charged(pip: 1, pim: 1, ep: 1)
  s.missing :nu_e          # massless neutrino
end

# Four-momentum-constrained kinematic fit with the massless neutrino
alg_e.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end
# No equivalent PID note for the electron channel.

alg_e.apply
root_files_e = alg_e.execute_on([data_3773, incMC_3773, exMC_signal_e])