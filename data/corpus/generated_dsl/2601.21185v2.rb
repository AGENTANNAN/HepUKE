### Dataset preparation ###
dd_data  = DatasetManager.real_data.find("712_3773")      # psi(3770) real data, 20.3 fb^-1 (BOSS 7.1.2)
dd_incMC = DatasetManager.inclusive_mc.find("712_3773")   # corresponding inclusive MC (tag-reformed samples)

# --- Decay cards: one per signal mode (semileptonic D decay + hadronic opposite D) ---
decay_card_d0_e = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D0 anti-D0          PHSP;
  Enddecay

  Decay D0
  1.0000 K- e+ nu_e          PHSP;
  Enddecay

  Decay anti-D0
  1.0000 K+ pi- pi0          PHSP;
  Enddecay

  End
DECAYCARD

decay_card_d0_mu = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D0 anti-D0          PHSP;
  Enddecay

  Decay D0
  1.0000 K- mu+ nu_mu        PHSP;
  Enddecay

  Decay anti-D0
  1.0000 K+ pi- pi0          PHSP;
  Enddecay

  End
DECAYCARD

decay_card_dp_e = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D+ D-               PHSP;
  Enddecay

  Decay D+
  1.0000 anti-K0 e+ nu_e     PHSP;
  Enddecay

  Decay K_S0
  1.0000 pi+ pi-             PHSP;
  Enddecay

  Decay D-
  1.0000 K+ pi- pi-          PHSP;
  Enddecay

  End
DECAYCARD

decay_card_dp_mu = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D+ D-               PHSP;
  Enddecay

  Decay D+
  1.0000 anti-K0 mu+ nu_mu   PHSP;
  Enddecay

  Decay K_S0
  1.0000 pi+ pi-             PHSP;
  Enddecay

  Decay D-
  1.0000 K+ pi- pi-          PHSP;
  Enddecay

  End
DECAYCARD

# --- Exclusive signal MC: 1,000,000 events for each of the four signal modes ---
exMC_d0_e = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_d0_ke_nu"
  config.related_dataset = dd_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_d0_e
  config.cross_section   = :default
end

exMC_d0_mu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_d0_kmu_nu"
  config.related_dataset = dd_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_d0_mu
  config.cross_section   = :default
end

exMC_dp_e = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_dp_kbar0_e_nu"
  config.related_dataset = dd_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_dp_e
  config.cross_section   = :default
end

exMC_dp_mu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_dp_kbar0_mu_nu"
  config.related_dataset = dd_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_dp_mu
  config.cross_section   = :default
end

### Event selection (BOSS) — tag-based (single tag, opposite D hadronic; signal side semileptonic) ###
# Six hadronic tag modes for each tag species
d0_tag_modes = [:D0toKPi, :D0toKPiPi0, :D0toKPiPiPi, :D0toKsPiPi, :D0toKPiPi0Pi0, :D0toKsPiPiPi0]
dp_tag_modes = [:DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKPiPiPiPi0]

# ===== Signal mode 1: D0 -> K- e+ nu_e  (tag: opposite D0) =====
alg_d0_e = TagAnalysis.new("D0ToKEnuTag")
alg_d0_e.set_header(["D0ToKEnuTagAlg/D0ToKEnuTag.h"])
        .set_constant({ "ECMS" => [:double, 3.773] })
        .with_decay_card(decay_card_d0_e)

alg_d0_e.tag_side(:D0) do |t|
  t.modes(*d0_tag_modes)                       # six hadronic D0 tag modes
  t.window :mBC,    min: 1.859, max: 1.873     # tag mBC window (D0)
  t.window :deltaE, abs: 0.02                  # representative symmetric DeltaE window
end

alg_d0_e.signal_side do |s|
  s.charged(km: 1, ep: 1)                      # one K- and one e+ on the signal side
  s.require_charge 0                           # net charge zero
  s.missing :nu_e                              # missing neutrino (massless)
end

alg_d0_e.fit do |f|
  f.constrain_four_momentum                    # 4C kinematic fit
  f.chi2_cut 200                               # chi2 < 200
end

alg_d0_e.note(:tag_deltaE_window,
              "DeltaE tag windows are mode-dependent (+-3.5 sigma); a single representative symmetric window is declared here")
       .note(:tag_candidate_selection,
              "the tag candidate combination with the smallest |DeltaE| is retained")

alg_d0_e.apply
alg_d0_e.execute_on([dd_data, dd_incMC, exMC_d0_e])

# ===== Signal mode 2: D0 -> K- mu+ nu_mu  (tag: opposite D0) =====
alg_d0_mu = TagAnalysis.new("D0ToKMuNuTag")
alg_d0_mu.set_header(["D0ToKMuNuTagAlg/D0ToKMuNuTag.h"])
         .set_constant({ "ECMS" => [:double, 3.773] })
         .with_decay_card(decay_card_d0_mu)

alg_d0_mu.tag_side(:D0) do |t|
  t.modes(*d0_tag_modes)
  t.window :mBC,    min: 1.859, max: 1.873
  t.window :deltaE, abs: 0.02
end

alg_d0_mu.signal_side do |s|
  s.charged(km: 1, mup: 1)                     # one K- and one mu+
  s.require_charge 0
  s.missing :nu_mu                             # missing neutrino (massless)
end

alg_d0_mu.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
  f.invariant_mass_of(:km, :mup).between(0.0, 1.56)   # muon-mode suppression M(K- mu+) < 1.56 GeV
end

alg_d0_mu.note(:tag_deltaE_window,
               "DeltaE tag windows are mode-dependent (+-3.5 sigma); a single representative symmetric window is declared here")
        .note(:tag_candidate_selection,
               "the tag candidate combination with the smallest |DeltaE| is retained")

alg_d0_mu.apply
alg_d0_mu.execute_on([dd_data, dd_incMC, exMC_d0_mu])

# ===== Signal mode 3: D+ -> Kbar0 e+ nu_e (Kbar0 -> K_S0 -> pi+ pi-); tag: opposite D+ =====
alg_dp_e = TagAnalysis.new("DpToKbar0EnuTag")
alg_dp_e.set_header(["DpToKbar0EnuTagAlg/DpToKbar0EnuTag.h"])
        .set_constant({ "ECMS" => [:double, 3.773] })
        .with_decay_card(decay_card_dp_e)

alg_dp_e.tag_side(:Dplus) do |t|
  t.modes(*dp_tag_modes)                       # six hadronic D+ tag modes
  t.window :mBC,    min: 1.863, max: 1.877     # tag mBC window (D+)
  t.window :deltaE, abs: 0.02
end

alg_dp_e.signal_side do |s|
  s.charged(pip: 1, pim: 1, ep: 1)             # pi+ pi- (from K_S0) and one e+
  s.require_charge 1                           # net charge +1
  s.missing :nu_e
end

alg_dp_e.fit do |f|
  f.constrain_four_momentum                    # 4C kinematic fit
  f.chi2_cut 200
  f.invariant_mass_of(:pip, :pim).constrain_to_nominal_mass_of(:K_S0)  # K_S0 mass constraint
end

alg_dp_e.note(:tag_deltaE_window,
              "DeltaE tag windows are mode-dependent (+-3.5 sigma); a single representative symmetric window is declared here")
       .note(:tag_candidate_selection,
              "the tag candidate combination with the smallest |DeltaE| is retained")

alg_dp_e.apply
alg_dp_e.execute_on([dd_data, dd_incMC, exMC_dp_e])

# ===== Signal mode 4: D+ -> Kbar0 mu+ nu_mu (Kbar0 -> K_S0 -> pi+ pi-); tag: opposite D+ =====
alg_dp_mu = TagAnalysis.new("DpToKbar0MuNuTag")
alg_dp_mu.set_header(["DpToKbar0MuNuTagAlg/DpToKbar0MuNuTag.h"])
         .set_constant({ "ECMS" => [:double, 3.773] })
         .with_decay_card(decay_card_dp_mu)

alg_dp_mu.tag_side(:Dplus) do |t|
  t.modes(*dp_tag_modes)
  t.window :mBC,    min: 1.863, max: 1.877
  t.window :deltaE, abs: 0.02
end

alg_dp_mu.signal_side do |s|
  s.charged(pip: 1, pim: 1, mup: 1)            # pi+ pi- (from K_S0) and one mu+
  s.require_charge 1                           # net charge +1
  s.missing :nu_mu
end

alg_dp_mu.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
  f.invariant_mass_of(:pip, :pim).constrain_to_nominal_mass_of(:K_S0)   # K_S0 mass constraint
  f.invariant_mass_of(:pip, :pim, :mup).between(0.0, 1.59)              # M(K_S0 mu+) < 1.59 GeV
end

alg_dp_mu.note(:tag_deltaE_window,
               "DeltaE tag windows are mode-dependent (+-3.5 sigma); a single representative symmetric window is declared here")
        .note(:tag_candidate_selection,
               "the tag candidate combination with the smallest |DeltaE| is retained")

alg_dp_mu.apply
alg_dp_mu.execute_on([dd_data, dd_incMC, exMC_dp_mu])