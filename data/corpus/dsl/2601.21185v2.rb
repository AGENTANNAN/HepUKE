# BESIII analysis: D^{0(+)} -> Kbar l+ nu_l (l = e, mu) at sqrt(s)=3.773 GeV.
# 20.3 fb^-1 psi(3770) data, double-tag (DT) technique.
# 4 signal channels:
#   D0  -> K-       e+  nu_e
#   D0  -> K-       mu+ nu_mu
#   D+  -> Kbar0    e+  nu_e   (Kbar0 -> K_S0 -> pi+ pi-)
#   D+  -> Kbar0    mu+ nu_mu

### Dataset ###
psi3770_data  = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

### Decay cards ###
# D0 -> K- e+ nu_e
decay_card_D0_Kenu = <<~DECAYCARD
  Decay D*0
  1.0000  D0                                  PHSP;
  Enddecay

  Decay D0
  1.0000  K-  e+  nu_e                        ISGW2;
  Enddecay

  End
DECAYCARD

# D0 -> K- mu+ nu_mu
decay_card_D0_Kmunu = <<~DECAYCARD
  Decay D*0
  1.0000  D0                                  PHSP;
  Enddecay

  Decay D0
  1.0000  K-  mu+  nu_mu                      ISGW2;
  Enddecay

  End
DECAYCARD

# D+ -> Kbar0 e+ nu_e
decay_card_Dp_K0enu = <<~DECAYCARD
  Decay D+
  1.0000  anti-K0  e+  nu_e                   ISGW2;
  Enddecay

  Decay anti-K0
  1.0000  K_S0                                PHSP;
  Enddecay

  Decay K_S0
  1.0000  pi+  pi-                            PHSP;
  Enddecay

  End
DECAYCARD

# D+ -> Kbar0 mu+ nu_mu
decay_card_Dp_K0munu = <<~DECAYCARD
  Decay D+
  1.0000  anti-K0  mu+  nu_mu                 ISGW2;
  Enddecay

  Decay anti-K0
  1.0000  K_S0                                PHSP;
  Enddecay

  Decay K_S0
  1.0000  pi+  pi-                            PHSP;
  Enddecay

  End
DECAYCARD

### Exclusive MC ###
[ [decay_card_D0_Kenu,   "D0_Kmenu"],
  [decay_card_D0_Kmunu,  "D0_Kmmunu"],
  [decay_card_Dp_K0enu,  "Dp_K0enu"],
  [decay_card_Dp_K0munu, "Dp_K0munu"] ].each do |card, name|
  mc = DatasetManager.create_exclusive_mc do |config|
    config.sample_name    = name
    config.related_dataset = psi3770_data
    config.events         = 1_000_000
    config.decay_card     = card
    config.cross_section  = :default
  end
  mc.save_to_config(format: :yaml, file_path: "exMC_#{name}")
end

# Reload the exMC handles (each algorithm attaches its own).
exMC_D0_Kenu = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "D0_Kmenu"; c.related_dataset = psi3770_data
  c.events = 1_000_000; c.decay_card = decay_card_D0_Kenu; c.cross_section = :default
end
exMC_D0_Kmunu = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "D0_Kmmunu"; c.related_dataset = psi3770_data
  c.events = 1_000_000; c.decay_card = decay_card_D0_Kmunu; c.cross_section = :default
end
exMC_Dp_K0enu = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "Dp_K0enu"; c.related_dataset = psi3770_data
  c.events = 1_000_000; c.decay_card = decay_card_Dp_K0enu; c.cross_section = :default
end
exMC_Dp_K0munu = DatasetManager.create_exclusive_mc do |c|
  c.sample_name = "Dp_K0munu"; c.related_dataset = psi3770_data
  c.events = 1_000_000; c.decay_card = decay_card_Dp_K0munu; c.cross_section = :default
end

########################################################################
# Algorithm I : D0 -> K- e+ nu_e (Dbar0 hadronic single tag)
########################################################################
alg_D0Kenu = TagAnalysis.new("D0toKmEpNue")
alg_D0Kenu.set_header(["D0toKmEpNueAlg/D0toKmEpNue.h"])
          .set_constant({ "ECMS" => [:double, 3.773] })

alg_D0Kenu.tag_side(:D0) do |t|
  t.modes :D0toKPi,
          :D0toKPiPi0,
          :D0toKPiPiPi,
          :D0toKsPiPi,
          :D0toKPiPi0Pi0,
          :D0toKPiPiPiPi0
  t.charm -1                                # Dbar0 tag
end

alg_D0Kenu.signal_side do |s|
  s.charged(km: 1, ep: 1)                   # K- e+ on the signal side
  s.require_charge 0                         # -1 + 1 = 0
  s.missing :nu_e                            # massless neutrino
end

alg_D0Kenu.fit do |f|
  f.constrain_four_momentum                  # tag + K + e + nu = ecms_lab
  f.chi2_cut 200
  f.store_fitted_momenta
end

alg_D0Kenu
  .note(:tag_selection,
        "Tag D-bar0 candidates selected via mBC in [1.859, 1.873] GeV/c^2 and " \
        "mode-dependent DeltaE windows (~+/-3.5 sigma). If multiple " \
        "combinations exist, the one with smallest |DeltaE| is chosen. Windows " \
        "and yields taken from the companion paper (Ref. [60]).")
  .note(:Umiss_signal_extraction,
        "Signal yield extracted from binned max-likelihood fit to " \
        "Umiss = Emiss - |pmiss|; the D-meson momentum is constrained using " \
        "the tag direction: p_D = -phat_Dbar * sqrt(Ebeam^2 - m_Dbar^2).")
  .with_decay_card(decay_card_D0_Kenu)
  .apply

alg_D0Kenu.execute_on([psi3770_data, psi3770_incMC, exMC_D0_Kenu])

########################################################################
# Algorithm II : D0 -> K- mu+ nu_mu
########################################################################
alg_D0Kmunu = TagAnalysis.new("D0toKmMupNumu")
alg_D0Kmunu.set_header(["D0toKmMupNumuAlg/D0toKmMupNumu.h"])
           .set_constant({ "ECMS" => [:double, 3.773] })

alg_D0Kmunu.tag_side(:D0) do |t|
  t.modes :D0toKPi,
          :D0toKPiPi0,
          :D0toKPiPiPi,
          :D0toKsPiPi,
          :D0toKPiPi0Pi0,
          :D0toKPiPiPiPi0
  t.charm -1
end

alg_D0Kmunu.signal_side do |s|
  s.charged(km: 1, mup: 1)                  # K- mu+ on the signal side
  s.require_charge 0
  s.missing :nu_mu                           # massless neutrino
end

alg_D0Kmunu.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
  f.store_fitted_momenta
end

alg_D0Kmunu
  .note(:tag_selection,
        "Tag D-bar0 selection as in Algorithm I (mBC [1.859,1.873] GeV/c^2, " \
        "mode-dependent DeltaE windows).")
  .note(:mkmuoncut,
        "Optimised background suppression: M(K- mu+) < 1.56 GeV/c^2 applied " \
        "in the ROOT stage as a cut on the stored invariant mass.")
  .note(:Umiss_signal_extraction,
        "Signal yield extracted from Umiss fit; D-meson momentum constrained " \
        "using the tag direction.")
  .with_decay_card(decay_card_D0_Kmunu)
  .apply

alg_D0Kmunu.execute_on([psi3770_data, psi3770_incMC, exMC_D0_Kmunu])

########################################################################
# Algorithm III : D+ -> Kbar0 e+ nu_e (Kbar0 -> K_S0 -> pi+ pi-)
########################################################################
alg_DpK0enu = TagAnalysis.new("DptoK0EpNue")
alg_DpK0enu.set_header(["DptoK0EpNueAlg/DptoK0EpNue.h"])
           .set_constant({ "ECMS" => [:double, 3.773] })

alg_DpK0enu.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi,
          :DptoKsPi,
          :DptoKPiPiPi0,
          :DptoKsPiPi0,
          :DptoKsPiPiPi,
          :DptoKKPi
  t.charm -1                                # D- hadronic tag
end

alg_DpK0enu.signal_side do |s|
  s.charged(pip: 1, pim: 1, ep: 1)          # pi+ pi- (from K_S0) + e+
  s.require_charge 1                         # 0 + 1 = +1 (D+)
  s.missing :nu_e
end

alg_DpK0enu.fit do |f|
  f.constrain_four_momentum
  # K_S0 mass constraint on the (pi+, pi-) system
  f.invariant_mass_of(:pip, :pim).constrain_to_nominal_mass_of(:K_S0)
  f.chi2_cut 200
  f.store_fitted_momenta
end

alg_DpK0enu
  .note(:tag_selection,
        "Tag D- selection: mBC in [1.863, 1.877] GeV/c^2, mode-dependent " \
        "DeltaE windows (~+/-3.5 sigma). If multiple combinations exist, the " \
        "one with smallest |DeltaE| is chosen.")
  .note(:Ks_reconstruction,
        "Kbar0 reconstructed via K_S0 -> pi+ pi- with standard BESIII K_S0 " \
        "vertex-fit selection; the signal pi+ pi- are the K_S0 daughters.")
  .note(:Umiss_signal_extraction,
        "Signal yield extracted from Umiss fit; D-meson momentum constrained " \
        "using the tag direction.")
  .with_decay_card(decay_card_Dp_K0enu)
  .apply

alg_DpK0enu.execute_on([psi3770_data, psi3770_incMC, exMC_Dp_K0enu])

########################################################################
# Algorithm IV : D+ -> Kbar0 mu+ nu_mu
########################################################################
alg_DpK0munu = TagAnalysis.new("DptoK0MupNumu")
alg_DpK0munu.set_header(["DptoK0MupNumuAlg/DptoK0MupNumu.h"])
            .set_constant({ "ECMS" => [:double, 3.773] })

alg_DpK0munu.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi,
          :DptoKsPi,
          :DptoKPiPiPi0,
          :DptoKsPiPi0,
          :DptoKsPiPiPi,
          :DptoKKPi
  t.charm -1
end

alg_DpK0munu.signal_side do |s|
  s.charged(pip: 1, pim: 1, mup: 1)
  s.require_charge 1
  s.missing :nu_mu
end

alg_DpK0munu.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:pip, :pim).constrain_to_nominal_mass_of(:K_S0)
  f.chi2_cut 200
  f.store_fitted_momenta
end

alg_DpK0munu
  .note(:tag_selection,
        "Tag D- selection as in Algorithm III (mBC [1.863,1.877], mode-" \
        "dependent DeltaE windows).")
  .note(:Ks_reconstruction,
        "Kbar0 -> K_S0 -> pi+ pi- standard vertex-fit selection.")
  .note(:mk0muoncut,
        "Optimised background suppression: M(K_S0 mu+) < 1.59 GeV/c^2 applied " \
        "in the ROOT stage.")
  .note(:Umiss_signal_extraction,
        "Signal yield extracted from Umiss fit; D-meson momentum constrained " \
        "using the tag direction.")
  .with_decay_card(decay_card_Dp_K0munu)
  .apply

alg_DpK0munu.execute_on([psi3770_data, psi3770_incMC, exMC_Dp_K0munu])
