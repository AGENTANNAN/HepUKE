### Dataset description ###
# ψ(3770) at 3.773 GeV — real data and matching inclusive MC
data_3773 = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# Signal-side decay cards (the tag D0 -> K+pi-, K+pi-pi0, K+pi-pi+pi- modes are
# taken from the pre-stored DTag candidates; only the signal side is generated here).

# Sub-chain 1: D0 -> K_L0 phi, phi -> K+ K-
decay_card_KLS_phi = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 K_L0 phi PHSP;
    Enddecay

    Decay phi
    1.0000 K+ K- VSS;
    Enddecay

    End
DECAYCARD

# Sub-chain 2: D0 -> K_L0 eta, eta -> gamma gamma
decay_card_KLS_eta_gg = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 K_L0 eta PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Sub-chain 3: D0 -> K_L0 eta, eta -> pi+ pi- pi0
decay_card_KLS_eta_3pi = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 K_L0 eta PHSP;
    Enddecay

    Decay eta
    1.0000 pi+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Sub-chain 4: D0 -> K_L0 omega, omega -> pi+ pi- pi0
decay_card_KLS_omega = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 K_L0 omega PHSP;
    Enddecay

    Decay omega
    1.0000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Sub-chain 5: D0 -> K_L0 eta', eta' -> pi+ pi- eta, eta -> gamma gamma
decay_card_KLS_etap_eta = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 K_L0 eta' PHSP;
    Enddecay

    Decay eta'
    1.0000 pi+ pi- eta PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Sub-chain 6: D0 -> K_L0 eta', eta' -> gamma rho0, rho0 -> pi+ pi-
decay_card_KLS_etap_rho = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 K_L0 eta' PHSP;
    Enddecay

    Decay eta'
    1.0000 gamma rho0 PHSP;
    Enddecay

    Decay rho0
    1.0000 pi+ pi- VSS;
    Enddecay

    End
DECAYCARD

# Dedicated signal-side exclusive MC for each sub-decay chain
exMC_KLS_phi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0toKLSphi"
  config.related_dataset = data_3773
  config.events          = 100000
  config.decay_card      = decay_card_KLS_phi
  config.cross_section   = :default
end

exMC_KLS_eta_gg = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0toKLSeta_gg"
  config.related_dataset = data_3773
  config.events          = 100000
  config.decay_card      = decay_card_KLS_eta_gg
  config.cross_section   = :default
end

exMC_KLS_eta_3pi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0toKLSeta_3pi"
  config.related_dataset = data_3773
  config.events          = 100000
  config.decay_card      = decay_card_KLS_eta_3pi
  config.cross_section   = :default
end

exMC_KLS_omega = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0toKLSomega"
  config.related_dataset = data_3773
  config.events          = 100000
  config.decay_card      = decay_card_KLS_omega
  config.cross_section   = :default
end

exMC_KLS_etap_eta = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0toKLSetap_eta"
  config.related_dataset = data_3773
  config.events          = 100000
  config.decay_card      = decay_card_KLS_etap_eta
  config.cross_section   = :default
end

exMC_KLS_etap_rho = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0toKLSetap_rho"
  config.related_dataset = data_3773
  config.events          = 100000
  config.decay_card      = decay_card_KLS_etap_rho
  config.cross_section   = :default
end

### Event selection (BOSS) — double-tag analysis ###

# ---- Channel 1: D0 -> K_L0 phi (phi -> K+ K-) ----
alg_phi = TagAnalysis.new("D0ToKLSphi")
alg_phi.set_header(["D0ToKLSphiAlg/D0ToKLSphi.h"])
       .set_constant({"ECMS" => [:double, 3.773]})
       .with_decay_card(decay_card_KLS_phi)
# Tag side: D0 -> K+pi-, K+pi-pi0, K+pi-pi+pi-; both charm assignments, ranked by invariant mass
alg_phi.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.rank_by :inv
end
# Signal side: K+ K- with K_L0 missing
alg_phi.signal_side do |s|
  s.charged(kp: 1, km: 1)
  s.missing :K_L0
end
# 4C fit, no resonance constraint
alg_phi.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end
alg_phi.apply
alg_phi.execute_on([data_3773, incMC_3773, exMC_KLS_phi])

# ---- Channel 2: D0 -> K_L0 eta (eta -> gamma gamma) ----
alg_eta_gg = TagAnalysis.new("D0ToKLSeta2Gam")
alg_eta_gg.set_header(["D0ToKLSeta2GamAlg/D0ToKLSeta2Gam.h"])
          .set_constant({"ECMS" => [:double, 3.773]})
          .with_decay_card(decay_card_KLS_eta_gg)
alg_eta_gg.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.rank_by :inv
end
alg_eta_gg.signal_side do |s|
  s.photons 2
  s.min_photon_angle 10.0
  s.missing :K_L0
end
alg_eta_gg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 200
end
alg_eta_gg.apply
alg_eta_gg.execute_on([data_3773, incMC_3773, exMC_KLS_eta_gg])

# ---- Channel 3: D0 -> K_L0 eta (eta -> pi+ pi- pi0) ----
alg_eta_3pi = TagAnalysis.new("D0ToKLSeta3Pi")
alg_eta_3pi.set_header(["D0ToKLSeta3PiAlg/D0ToKLSeta3Pi.h"])
           .set_constant({"ECMS" => [:double, 3.773]})
           .with_decay_card(decay_card_KLS_eta_3pi)
alg_eta_3pi.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.rank_by :inv
end
alg_eta_3pi.signal_side do |s|
  s.charged(pip: 1, pim: 1)
  s.photons 2
  s.min_photon_angle 10.0
  s.missing :K_L0
end
alg_eta_3pi.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end
alg_eta_3pi.apply
alg_eta_3pi.execute_on([data_3773, incMC_3773, exMC_KLS_eta_3pi])

# ---- Channel 4: D0 -> K_L0 omega (omega -> pi+ pi- pi0) ----
alg_omega = TagAnalysis.new("D0ToKLSomega")
alg_omega.set_header(["D0ToKLSomegaAlg/D0ToKLSomega.h"])
         .set_constant({"ECMS" => [:double, 3.773]})
         .with_decay_card(decay_card_KLS_omega)
alg_omega.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.rank_by :inv
end
alg_omega.signal_side do |s|
  s.charged(pip: 1, pim: 1)
  s.photons 2
  s.min_photon_angle 10.0
  s.missing :K_L0
end
alg_omega.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end
alg_omega.apply
alg_omega.execute_on([data_3773, incMC_3773, exMC_KLS_omega])

# ---- Channel 5: D0 -> K_L0 eta' (eta' -> pi+ pi- eta, eta -> gamma gamma) ----
alg_etap_eta = TagAnalysis.new("D0ToKLSetapToPiPiEta")
alg_etap_eta.set_header(["D0ToKLSetapToPiPiEtaAlg/D0ToKLSetapToPiPiEta.h"])
            .set_constant({"ECMS" => [:double, 3.773]})
            .with_decay_card(decay_card_KLS_etap_eta)
alg_etap_eta.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.rank_by :inv
end
alg_etap_eta.signal_side do |s|
  s.charged(pip: 1, pim: 1)
  s.photons 2
  s.min_photon_angle 10.0
  s.missing :K_L0
end
alg_etap_eta.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 200
end
alg_etap_eta.apply
alg_etap_eta.execute_on([data_3773, incMC_3773, exMC_KLS_etap_eta])

# ---- Channel 6: D0 -> K_L0 eta' (eta' -> gamma rho0, rho0 -> pi+ pi-) ----
alg_etap_rho = TagAnalysis.new("D0ToKLSetapToGammaRho")
alg_etap_rho.set_header(["D0ToKLSetapToGammaRhoAlg/D0ToKLSetapToGammaRho.h"])
            .set_constant({"ECMS" => [:double, 3.773]})
            .with_decay_card(decay_card_KLS_etap_rho)
alg_etap_rho.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.rank_by :inv
end
alg_etap_rho.signal_side do |s|
  s.charged(pip: 1, pim: 1)
  s.photons 1
  s.min_photon_angle 10.0
  s.missing :K_L0
end
alg_etap_rho.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end
alg_etap_rho.apply
alg_etap_rho.execute_on([data_3773, incMC_3773, exMC_KLS_etap_rho])