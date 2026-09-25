# Paper: 2202.04232v2 — Ds+ → π+π0η′ amplitude analysis & BF measurement (TagAnalysis DT)
# Data: BOSS 703 at 4.178–4.226 GeV (6 energy points, 6.32 fb−1)
# TagAnalysis surface — no Selection object, apply takes no argument
# Store-not-cut: tag mBC/ΔE stored unconditionally

### Dataset preparation — 6 energy points ###
data_4180 = DatasetManager.real_data.find("703_4180")
data_4190 = DatasetManager.real_data.find("703_4190")
data_4200 = DatasetManager.real_data.find("703_4200")
data_4210 = DatasetManager.real_data.find("703_4210")
data_4220 = DatasetManager.real_data.find("703_4220")
data_4230 = DatasetManager.real_data.find("703_4230")

data_points = [data_4180, data_4190, data_4200, data_4210, data_4220, data_4230]

incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")
incMC_4190 = DatasetManager.inclusive_mc.find("703_4190")
incMC_4200 = DatasetManager.inclusive_mc.find("703_4200")
incMC_4210 = DatasetManager.inclusive_mc.find("703_4210")
incMC_4220 = DatasetManager.inclusive_mc.find("703_4220")
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")

incMC_points = [incMC_4180, incMC_4190, incMC_4200, incMC_4210, incMC_4220, incMC_4230]

### Decay card for signal MC — e+e- → Ds*± Ds∓ → γ Ds+ Ds- ###
decay_card_for_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 D_s*- D_s+ PHSP;
    Enddecay

    Decay D_s*-
    1.0000 gamma D_s- PHSP;
    Enddecay

    Decay D_s+
    1.0000 pi+ pi0 eta' PHSP;
    Enddecay

    Decay D_s-
    1.0000 K+ K- pi- PHSP;
    Enddecay

    Decay eta'
    1.0000 pi+ pi- eta PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_DsToPiPi0EtaP"
  config.related_dataset = data_4180
  config.events = 500000
  config.decay_card = decay_card_for_signal
  config.cross_section = :default
end

### TagAnalysis ###
alg = TagAnalysis.new("DsDTagPiPi0EtaP")

alg.set_header(["DsDTagPiPi0EtaPAlg/DsDTagPiPi0EtaP.h"])
   .set_constant({ "ECMS" => [:double, 4.226] })
   .with_decay_card(decay_card_for_signal)

# Tag side 1: Ds- reconstructed via 12 hadronic tag modes
alg.tag_side(:Ds) do |t|
  t.modes :DstoKsK, :DstoKKPi, :DstoKsKPi0,
          :DstoKsKKPiPi, :DstoKsKPiPiPi0,
          :Dstopi3pi0, :Dstoppi0eta,
          :DstoKKpipi, :DstoKsKpipi0,
          :Dstopi3pi, :Dstopi3pi0eta, :Dstopieta
  t.charm -1
end

# Signal side: Ds+ → π+ π0 η′, η′ → π+ π− η, η → γγ
# Total signal daughters: 2 π+, 1 π−, ≥4 photons (2 for π0 → γγ, 2 for η → γγ)
alg.signal_side do |s|
  s.photons 4
  s.charged(pip: 2, pim: 1)
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

# 9C kinematic fit: 4C + π0 mass + η mass + η′ mass + tag Ds mass + Ds* mass
# Note: η′ mass constraint over (pip, pim, gamma, gamma) and Ds* mass constraint
# are not fully expressible in the current TagAnalysis DSL v1.
# The DSL provides: 4C + π0 mass (γγ) + η mass (γγ) + tag Ds mass.
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:Ds)
  f.chi2_cut 200
end

alg.apply
root_files = alg.execute_on(data_points + incMC_points + [exMC_signal])

# Note: Full 9C fit in the paper constrains π0, η, η′, tag Ds−, and Ds*± masses
# in addition to 4-momentum conservation. The η′ and Ds* constraints are not
# directly expressible in TagAnalysis DSL v1 (derived participants don't
# include intermediate particles like :pi0 or :eta as named symbols,
# and the Ds* soft photon is bundled in the tag reconstruction).
# Post-fit analysis: tag mBC/ΔE cuts, signal Ds+ mass window [1.92, 2.00] GeV/c2,
# fake-η veto, χ2_9C ranking — all applied at ROOT analysis stage.