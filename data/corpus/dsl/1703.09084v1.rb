# BESIII: D+ -> Kbar0 e+ nu_e and D+ -> pi0 e+ nu_e
# Data: 2.93 fb^-1 at sqrt(s)=3.773 GeV (psi(3770)); D-tag technique.

### Dataset description ###
psi3770_data   = DatasetManager.real_data.find("712_3773")
psi3770_incMC  = DatasetManager.inclusive_mc.find("712_3773")

# Decay card for signal D+ -> Kbar0 e+ nu_e (K0bar reconstructed via K0S -> pi+ pi-)
decay_card_K0enu = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D+ D-                              PHSP;
  Enddecay

  Decay D+
  1.0000 anti-K0 e+ nu_e                    PHOTOS ISGW2;
  Enddecay

  Decay anti-K0
  1.0000 K_S0                               PHSP;
  Enddecay

  Decay K_S0
  1.0000 pi+ pi-                            PHSP;
  Enddecay

  Decay D-
  1.0000 K+ pi- pi-                         PHSP;
  Enddecay

  End
DECAYCARD

# Decay card for signal D+ -> pi0 e+ nu_e
decay_card_pi0enu = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D+ D-                              PHSP;
  Enddecay

  Decay D+
  1.0000 pi0 e+ nu_e                        PHOTOS ISGW2;
  Enddecay

  Decay pi0
  1.0000 gamma gamma                        PHSP;
  Enddecay

  Decay D-
  1.0000 K+ pi- pi-                         PHSP;
  Enddecay

  End
DECAYCARD

exMC_K0enu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "Dp_to_K0bar_ep_nue"
  config.related_dataset = psi3770_data
  config.events         = 500000
  config.decay_card     = decay_card_K0enu
  config.cross_section  = :default
end

exMC_pi0enu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "Dp_to_pi0_ep_nue"
  config.related_dataset = psi3770_data
  config.events         = 500000
  config.decay_card     = decay_card_pi0enu
  config.cross_section  = :default
end

### Event selection (BOSS) — D-tag with recoil semileptonic reconstruction ###

# ---------------- Algorithm 1: D+ -> Kbar0 e+ nu_e ----------------
alg_K0enu = TagAnalysis.new("DpToK0enu")
alg_K0enu.set_header(["DpToK0enuAlg/DpToK0enu.h"])
         .set_constant({ "ECMS" => [:double, 3.773] })
         .with_decay_card(decay_card_K0enu)

alg_K0enu.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKsK, :DptoKKPi,
          :DptoKPiPiPi0, :DptoPiPiPi, :DptoKsPiPi0,
          :DptoKPiPiPiPi, :DptoKsPiPiPi
  t.charm -1
end

alg_K0enu.signal_side do |s|
  # Signal side: K0S -> pi+ pi-, plus e+; charge conservation on signal side = +1
  s.charged(pip: 1, pim: 1, ep: 1)
  s.require_charge 1
  s.missing :nu_e            # massless neutrino (semileptonic)
  s.min_photon_angle 10.0
end

alg_K0enu.fit do |f|
  f.constrain_four_momentum
  # K0S mass constraint on the two signal-side pions
  f.invariant_mass_of(:pip, :pim).constrain_to_nominal_mass_of(:K_S0)
  f.chi2_cut 200
end

alg_K0enu.note(:fsr_recovery,
               "FSR recovery: 4-momenta of photons within 5 degrees of the positron " \
               "direction are added to the positron 4-momentum before kinematic fit.")
        .note(:e_gamma_max_veto,
               "Additional background suppression: reject events with any unused " \
               "photon energy E_gamma,max > 300 MeV (extra-shower veto).")

alg_K0enu.apply
alg_K0enu.execute_on([psi3770_data, psi3770_incMC, exMC_K0enu])

# ---------------- Algorithm 2: D+ -> pi0 e+ nu_e ----------------
alg_pi0enu = TagAnalysis.new("DpToPi0enu")
alg_pi0enu.set_header(["DpToPi0enuAlg/DpToPi0enu.h"])
          .set_constant({ "ECMS" => [:double, 3.773] })
          .with_decay_card(decay_card_pi0enu)

alg_pi0enu.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKsK, :DptoKKPi,
          :DptoKPiPiPi0, :DptoPiPiPi, :DptoKsPiPi0,
          :DptoKPiPiPiPi, :DptoKsPiPiPi
  t.charm -1
end

alg_pi0enu.signal_side do |s|
  # Signal side: e+ track + at least 2 photons for pi0 -> gamma gamma
  s.charged(ep: 1)
  s.require_charge 1
  s.photons 2
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
  s.missing :nu_e
end

alg_pi0enu.fit do |f|
  f.constrain_four_momentum
  # 1-C mass constraint of gamma gamma to pi0 nominal mass
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

alg_pi0enu.note(:fsr_recovery,
                "FSR recovery: 4-momenta of photons within 5 degrees of the positron " \
                "direction are added to the positron 4-momentum before kinematic fit.")
          .note(:e_gamma_max_veto,
                "Additional background suppression: reject events with any unused " \
                "photon energy E_gamma,max > 300 MeV.")
          .note(:pi0_photon_endcap_energy,
                "Endcap photon energy threshold 50 MeV; barrel 25 MeV; shower time " \
                "within 700 ns of event start; M(gamma gamma) in (0.110, 0.150) GeV/c^2 " \
                "before 1-C fit; best pi0 by minimum 1-C chi^2.")

alg_pi0enu.apply
alg_pi0enu.execute_on([psi3770_data, psi3770_incMC, exMC_pi0enu])
