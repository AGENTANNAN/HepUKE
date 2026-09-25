### Dataset preparation ###
data_3773 = DatasetManager.real_data.find("712_3773")        # psi(3770) real data at 3.773 GeV
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")    # corresponding inclusive MC

# Decay card for the D0 -> pi- pi0 e+ nu_e signal (hadronic tag side = anti-D0)
decay_card_d0 = <<~DECAYCARD
    Decay psi(3770)
    1.000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 pi- pi0 e+ nu_e PHSP;
    Enddecay

    Decay anti-D0
    0.200 K+ pi- PHSP;
    0.200 K+ pi- pi0 PHSP;
    0.200 K+ pi- pi0 pi0 PHSP;
    0.200 K+ pi- pi- pi+ PHSP;
    0.200 K+ pi- pi- pi+ pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for the D+ -> pi- pi+ e+ nu_e signal (hadronic tag side = D-)
decay_card_dplus = <<~DECAYCARD
    Decay psi(3770)
    1.000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.000 pi- pi+ e+ nu_e PHSP;
    Enddecay

    Decay D-
    0.1667 K+ pi- pi- PHSP;
    0.1667 K+ pi- pi- pi0 PHSP;
    0.1666 K_S0 pi- PHSP;
    0.1667 K_S0 pi- pi0 PHSP;
    0.1667 K_S0 pi- pi+ pi- PHSP;
    0.1666 K+ K- pi- PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 200k exclusive-MC events for each signal mode
exMC_d0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_D0_pipi0enu"
  config.related_dataset = data_3773
  config.events          = 200_000
  config.decay_card      = decay_card_d0
  config.cross_section   = :default
end

exMC_dplus = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_Dplus_pipipenu"
  config.related_dataset = data_3773
  config.events          = 200_000
  config.decay_card      = decay_card_dplus
  config.cross_section   = :default
end

### Event selection (BOSS) — tag-based (semileptonic ST + missing nu) ###
# --- Mode I: D0 -> pi- pi0 e+ nu_e, tagged by the opposite-charm anti-D0 ---
alg_d0 = TagAnalysis.new("D0SemilepTag")
alg_d0.set_header(["D0SemilepTagAlg/D0SemilepTag.h"])
      .set_constant({"ECMS" => [:double, 3.773]})
      .with_decay_card(decay_card_d0)

# Tag side: hadronic anti-D0 (charm -1)
alg_d0.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPi0Pi0, :D0toKPiPiPi, :D0toKPiPiPiPi0
  t.charm -1
end

# Signal side: 2 photons (pi0), pi-, e+, net charge 0, missing nu_e
alg_d0.signal_side do |s|
  s.photons 2                 # pi0 -> gamma gamma
  s.charged(pim: 1, ep: 1)
  s.require_charge 0
  s.missing :nu_e             # massless semileptonic neutrino
end

# 4C fit; the two signal photons are constrained to the nominal pi0 mass
alg_d0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

# BOSS-side procedures with no DSL counterpart
alg_d0.note(:extra_photon_veto, "veto events with an extra photon above E_gamma_max 0.25 GeV " \
                                 "(photons not used by the signal pi0); BOSS-level, not expressible in the tag DSL")
     .note(:bremsstrahlung_recovery, "Bremsstrahlung / FSR photon recovery applied to the signal e+ " \
                                      "candidate (photons emitted close to the e+ track are added back to " \
                                      "its four-momentum before the kinematic fit)")
     .note(:umiss_window, "|Umiss| < 0.06 GeV window (Umiss = Emiss - |Pmiss|, from the fitted missing " \
                           "nu_e) defines the PWA sample; m_Umiss is stored by the tag fit and windowed in ROOT")

alg_d0.apply
alg_d0.execute_on([data_3773, incMC_3773, exMC_d0])

# --- Mode II: D+ -> pi- pi+ e+ nu_e, tagged by D- ---
alg_dplus = TagAnalysis.new("DplusSemilepTag")
alg_dplus.set_header(["DplusSemilepTagAlg/DplusSemilepTag.h"])
         .set_constant({"ECMS" => [:double, 3.773]})
         .with_decay_card(decay_card_dplus)

# Tag side: hadronic D- (charm -1)
alg_dplus.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
end

# Signal side: pi-, pi+, e+, net charge 1, missing nu_e
alg_dplus.signal_side do |s|
  s.charged(pim: 1, pip: 1, ep: 1)
  s.require_charge 1
  s.missing :nu_e
end

# 4C fit
alg_dplus.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

# BOSS-side procedures with no DSL counterpart
alg_dplus.note(:extra_photon_veto, "veto events with an extra photon above E_gamma_max 0.25 GeV " \
                                    "(photons not used by the tag); BOSS-level, not expressible in the tag DSL")
         .note(:bremsstrahlung_recovery, "Bremsstrahlung / FSR photon recovery applied to the signal e+ " \
                                          "candidate before the kinematic fit")
         .note(:background_veto, "K_S0 veto for the D+ mode: events with m(pi+ pi-) within 70 MeV/c^2 of the " \
                                  "K_S0 nominal mass (0.4976 GeV/c^2) are vetoed")
         .note(:umiss_window, "|Umiss| < 0.06 GeV window (Umiss = Emiss - |Pmiss|, from the fitted missing " \
                               "nu_e) defines the PWA sample; m_Umiss is stored by the tag fit and windowed in ROOT")

alg_dplus.apply
alg_dplus.execute_on([data_3773, incMC_3773, exMC_dplus])