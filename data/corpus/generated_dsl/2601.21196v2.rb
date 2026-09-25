# ==============================================================
# Dataset preparation : psi(3770) real data + matching inclusive MC
# ==============================================================
psipp_data  = DatasetManager.real_data.find("712_3773")      # 20.3 fb^-1 psi(3770) real data (BOSS 712), sqrt(s)=3.773 GeV
psipp_incMC = DatasetManager.inclusive_mc.find("712_3773")   # matching inclusive MC at 3.773 GeV

# --------------------------------------------------------------
# Decay cards for the four semileptonic signal modes
# (e+e- -> D Dbar at 3.773 GeV; signal D -> semileptonic, opposite D -> generic hadronic tag)
# --------------------------------------------------------------
decay_card_D0_Kenu = <<~DECAYCARD
  Decay psi(3770)
  1.0 D0 anti-D0 PHSP;
  Enddecay

  Decay D0
  1.0 K- e+ nu_e PHOTOS PHSP;
  Enddecay

  Decay anti-D0
  1.0 K+ pi- PHSP;
  Enddecay

  End
DECAYCARD

decay_card_D0_Kmunu = <<~DECAYCARD
  Decay psi(3770)
  1.0 D0 anti-D0 PHSP;
  Enddecay

  Decay D0
  1.0 K- mu+ nu_mu PHOTOS PHSP;
  Enddecay

  Decay anti-D0
  1.0 K+ pi- PHSP;
  Enddecay

  End
DECAYCARD

decay_card_Dp_K0enu = <<~DECAYCARD
  Decay psi(3770)
  1.0 D+ D- PHSP;
  Enddecay

  Decay D+
  1.0 anti-K0 e+ nu_e PHOTOS PHSP;
  Enddecay

  Decay anti-K0
  1.0 K_S0 PHSP;
  Enddecay

  Decay K_S0
  1.0 pi+ pi- PHSP;
  Enddecay

  Decay D-
  1.0 K+ pi- pi- PHSP;
  Enddecay

  End
DECAYCARD

decay_card_Dp_K0munu = <<~DECAYCARD
  Decay psi(3770)
  1.0 D+ D- PHSP;
  Enddecay

  Decay D+
  1.0 anti-K0 mu+ nu_mu PHOTOS PHSP;
  Enddecay

  Decay anti-K0
  1.0 K_S0 PHSP;
  Enddecay

  Decay K_S0
  1.0 pi+ pi- PHSP;
  Enddecay

  Decay D-
  1.0 K+ pi- pi- PHSP;
  Enddecay

  End
DECAYCARD

# --------------------------------------------------------------
# Exclusive MC : 2 million events for each signal mode
# --------------------------------------------------------------
exMC_D0_Kenu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psipp_D0_Kenu"          # D0 -> K- e+ nu_e
  config.related_dataset = psipp_data
  config.events          = 2_000_000
  config.decay_card      = decay_card_D0_Kenu
  config.cross_section   = :default
end

exMC_D0_Kmunu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psipp_D0_Kmunu"         # D0 -> K- mu+ nu_mu
  config.related_dataset = psipp_data
  config.events          = 2_000_000
  config.decay_card      = decay_card_D0_Kmunu
  config.cross_section   = :default
end

exMC_Dp_K0enu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psipp_Dp_K0enu"         # D+ -> Kbar0 e+ nu_e, Kbar0 -> K_S0 -> pi+pi-
  config.related_dataset = psipp_data
  config.events          = 2_000_000
  config.decay_card      = decay_card_Dp_K0enu
  config.cross_section   = :default
end

exMC_Dp_K0munu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psipp_Dp_K0munu"        # D+ -> Kbar0 mu+ nu_mu, Kbar0 -> K_S0 -> pi+pi-
  config.related_dataset = psipp_data
  config.events          = 2_000_000
  config.decay_card      = decay_card_Dp_K0munu
  config.cross_section   = :default
end

# ==============================================================
# Double-tag selection : the opposite D is fully reconstructed in
# six hadronic tag modes; the other D is the semileptonic signal.
# (one TagAnalysis per independent signal final state)
# ==============================================================

# --------------------------------------------------------------
# (1) Signal D0 -> K- e+ nu_e ; tag Dbar0
# --------------------------------------------------------------
alg_D0_Kenu = TagAnalysis.new("D0KenuTag")
alg_D0_Kenu.set_header(["D0KenuTagAlg/D0KenuTag.h"])
           .set_constant({"ECMS" => [:double, 3.773]})   # sqrt(s) = 3.773 GeV
           .with_decay_card(decay_card_D0_Kenu)

alg_D0_Kenu.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toK3Pi, :D0toKsPiPi, :D0toKPiPi0Pi0, :D0toK3PiPi0   # six hadronic Dbar0 tag modes
  t.charm -1                                                     # tag the Dbar0 (opposite to the D0 signal)
  t.window :deltaE, abs: 0.05                                    # representative |DeltaE| < 50 MeV (per-mode windows in ROOT)
  t.window :mBC, min: 1.859, max: 1.873                          # tag M_BC window for Dbar0
end

alg_D0_Kenu.signal_side do |s|
  s.charged(km: 1, ep: 1)          # exactly one K- and one e+, no extra charged tracks
  s.require_charge 0               # net charge of the D0 signal side
  s.missing :nu_e                  # one missing neutrino
  s.min_photon_energy 0.025        # showers above 25 MeV
end

alg_D0_Kenu.fit do |f|
  f.constrain_four_momentum        # 4C kinematic fit
  f.chi2_cut 200                   # chi2 < 200
end

alg_D0_Kenu
  .note(:tag_deltae_windows, "per-tag-mode DELTA E windows of roughly +-25-60 MeV are optimised per tag mode and applied in ROOT; the DSL tag-side window is per-side, so a single representative |DeltaE| < 50 MeV window is declared here")
  .note(:signal_lepton_pid, "electron identification uses L_e > 0.001 and L_e > 0.8*(L_e+L_pi+L_K); the generated signal-side lepton PID uses fixed v1 thresholds, so these custom criteria are applied in the ROOT analysis")
  .note(:extra_photon_energy, "total energy of leftover (unused) photons required below 0.25 GeV; applied as an event-level cut in ROOT since no signal-side declaration exists")

alg_D0_Kenu.apply
alg_D0_Kenu.execute_on([psipp_data, psipp_incMC, exMC_D0_Kenu])

# --------------------------------------------------------------
# (2) Signal D0 -> K- mu+ nu_mu ; tag Dbar0
# --------------------------------------------------------------
alg_D0_Kmunu = TagAnalysis.new("D0KmunuTag")
alg_D0_Kmunu.set_header(["D0KmunuTagAlg/D0KmunuTag.h"])
            .set_constant({"ECMS" => [:double, 3.773]})
            .with_decay_card(decay_card_D0_Kmunu)

alg_D0_Kmunu.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toK3Pi, :D0toKsPiPi, :D0toKPiPi0Pi0, :D0toK3PiPi0
  t.charm -1
  t.window :deltaE, abs: 0.05
  t.window :mBC, min: 1.859, max: 1.873
end

alg_D0_Kmunu.signal_side do |s|
  s.charged(km: 1, mup: 1)         # exactly one K- and one mu+, no extra charged tracks
  s.require_charge 0
  s.missing :nu_mu                 # one missing neutrino
  s.min_photon_energy 0.025        # showers above 25 MeV
end

alg_D0_Kmunu.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_D0_Kmunu
  .note(:tag_deltae_windows, "per-tag-mode DELTA E windows of roughly +-25-60 MeV applied in ROOT; a single representative |DeltaE| < 50 MeV tag-side window is declared here")
  .note(:signal_lepton_pid, "muon identification requires L_mu > L_e, L_mu > 0.001 and 0.1 < E_EMC < 0.3 GeV; the generated signal-side lepton PID uses fixed v1 thresholds, so these custom criteria are applied in the ROOT analysis")
  .note(:extra_photon_energy, "total energy of leftover (unused) photons required below 0.25 GeV; applied as an event-level cut in ROOT")

alg_D0_Kmunu.apply
alg_D0_Kmunu.execute_on([psipp_data, psipp_incMC, exMC_D0_Kmunu])

# --------------------------------------------------------------
# (3) Signal D+ -> Kbar0 e+ nu_e (Kbar0 -> K_S0 -> pi+pi-) ; tag D-
# --------------------------------------------------------------
alg_Dp_K0enu = TagAnalysis.new("DpK0enuTag")
alg_Dp_K0enu.set_header(["DpK0enuTagAlg/DpK0enuTag.h"])
            .set_constant({"ECMS" => [:double, 3.773]})
            .with_decay_card(decay_card_Dp_K0enu)

alg_Dp_K0enu.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi   # six hadronic D- tag modes
  t.charm -1                                                     # tag the D- (opposite to the D+ signal)
  t.window :deltaE, abs: 0.05                                    # representative |DeltaE| < 50 MeV
  t.window :mBC, min: 1.863, max: 1.877                          # tag M_BC window for D-
end

alg_Dp_K0enu.signal_side do |s|
  s.charged(pip: 1, pim: 1, ep: 1)   # pi+ pi- from K_S0 plus e+, no extra charged tracks
  s.require_charge 1                 # net charge of the D+ signal side
  s.missing :nu_e                    # one missing neutrino
  s.min_photon_energy 0.025          # showers above 25 MeV
end

alg_Dp_K0enu.fit do |f|
  f.constrain_four_momentum                                                       # 4C kinematic fit
  f.invariant_mass_of(:pip, :pim).constrain_to_nominal_mass_of(:K_S0)             # K_S0 mass constraint on the pi+pi- pair
  f.chi2_cut 200                                                                  # chi2 < 200
end

alg_Dp_K0enu
  .note(:tag_deltae_windows, "per-tag-mode DELTA E windows of roughly +-25-60 MeV applied in ROOT; a single representative |DeltaE| < 50 MeV tag-side window is declared here")
  .note(:signal_lepton_pid, "electron identification uses L_e > 0.001 and L_e > 0.8*(L_e+L_pi+L_K); the generated signal-side lepton PID uses fixed v1 thresholds, so these custom criteria are applied in the ROOT analysis")
  .note(:ks0_reconstruction, "Kbar0 -> K_S0 -> pi+pi- rebuilt from two opposite-charge tracks without PID, requiring |cos(theta)| < 0.93, |Vz| < 20 cm, secondary-vertex chi2 < 100, decay length > 2 sigma and M(pi+pi-) in (0.487, 0.511) GeV; only the K_S0 nominal-mass constraint is expressible in the tag fit")
  .note(:extra_photon_energy, "total energy of leftover (unused) photons required below 0.25 GeV; applied as an event-level cut in ROOT")

alg_Dp_K0enu.apply
alg_Dp_K0enu.execute_on([psipp_data, psipp_incMC, exMC_Dp_K0enu])

# --------------------------------------------------------------
# (4) Signal D+ -> Kbar0 mu+ nu_mu (Kbar0 -> K_S0 -> pi+pi-) ; tag D-
# --------------------------------------------------------------
alg_Dp_K0munu = TagAnalysis.new("DpK0munuTag")
alg_Dp_K0munu.set_header(["DpK0munuTagAlg/DpK0munuTag.h"])
             .set_constant({"ECMS" => [:double, 3.773]})
             .with_decay_card(decay_card_Dp_K0munu)

alg_Dp_K0munu.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
  t.window :deltaE, abs: 0.05
  t.window :mBC, min: 1.863, max: 1.877
end

alg_Dp_K0munu.signal_side do |s|
  s.charged(pip: 1, pim: 1, mup: 1)  # pi+ pi- from K_S0 plus mu+, no extra charged tracks
  s.require_charge 1                 # net charge of the D+ signal side
  s.missing :nu_mu                   # one missing neutrino
  s.min_photon_energy 0.025          # showers above 25 MeV
end

alg_Dp_K0munu.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:pip, :pim).constrain_to_nominal_mass_of(:K_S0)   # K_S0 mass constraint on the pi+pi- pair
  f.chi2_cut 200
end

alg_Dp_K0munu
  .note(:tag_deltae_windows, "per-tag-mode DELTA E windows of roughly +-25-60 MeV applied in ROOT; a single representative |DeltaE| < 50 MeV tag-side window is declared here")
  .note(:signal_lepton_pid, "muon identification requires L_mu > L_e, L_mu > 0.001 and 0.1 < E_EMC < 0.3 GeV; the generated signal-side lepton PID uses fixed v1 thresholds, so these custom criteria are applied in the ROOT analysis")
  .note(:ks0_reconstruction, "Kbar0 -> K_S0 -> pi+pi- rebuilt from two opposite-charge tracks without PID, requiring |cos(theta)| < 0.93, |Vz| < 20 cm, secondary-vertex chi2 < 100, decay length > 2 sigma and M(pi+pi-) in (0.487, 0.511) GeV; only the K_S0 nominal-mass constraint is expressible in the tag fit")
  .note(:extra_photon_energy, "total energy of leftover (unused) photons required below 0.25 GeV; applied as an event-level cut in ROOT")

alg_Dp_K0munu.apply
alg_Dp_K0munu.execute_on([psipp_data, psipp_incMC, exMC_Dp_K0munu])