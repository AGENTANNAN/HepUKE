# Paper 2404.09219v1: Amplitude analysis of D0(+) -> pi+ pi-(0) eta at psi(3770)
# Uses DT tag method at 3.773 GeV, 7.9 fb^-1

# --- Datasets ---
psi3770_data = DatasetManager.real_data.find("712_3773")
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")

# --- Decay card for D0 -> pi+ pi- eta (signal) ---
decay_card_d0 = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.000 pi+ pi- eta PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# --- Decay card for D+ -> pi+ pi0 eta (signal) ---
decay_card_dp = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.000 pi+ pi0 eta PHSP;
    Enddecay

    Decay D-
    1.000 K+ pi- pi- PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# --- Exclusive MC ---
exMC_d0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_d0_pimeta"
  config.related_dataset = psi3770_data
  config.events = 200000
  config.decay_card = decay_card_d0
  config.cross_section = :default
end

exMC_dp = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_dp_pipi0eta"
  config.related_dataset = psi3770_data
  config.events = 200000
  config.decay_card = decay_card_dp
  config.cross_section = :default
end

# ==========================================================================
# Algorithm 1: D0 -> pi+ pi- eta
# Tag side: anti-D0 modes; signal side: pi+ pi- + gamma gamma (eta)
# ==========================================================================
alg_d0 = TagAnalysis.new("D0toPiPiEtaDTag")
alg_d0.set_header(["D0toPiPiEtaDTagAlg/D0toPiPiEtaDTag.h"])
       .set_constant({ "ECMS" => [:double, 3.773] })
       .note(:ks_veto, "events with |M(pi+pi-)-m(K_S0)| < 0.03 GeV/c^2 vetoed to suppress D0->K_S0 eta background; veto applied in ROOT analysis")
       .note(:bdtd_eta_cut, "BDTG classifier applied on eta candidate (M(gammagamma), chi2(eta), photon helicity angle) for background suppression; BDTG cut retains 83% signal and rejects 78% background for D0 channel")

alg_d0.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPi0Pi0, :D0toKPiPiPi
end

alg_d0.signal_side do |s|
  s.charged(pip: 1, pim: 1)
  s.photons 2
  s.min_photon_angle 10.0
end

alg_d0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 200
end

alg_d0.with_decay_card(decay_card_d0).apply
alg_d0.execute_on([psi3770_data, psi3770_incMC, exMC_d0])

# ==========================================================================
# Algorithm 2: D+ -> pi+ pi0 eta
# Tag side: D- modes; signal side: pi+ + gamma gamma (pi0) + gamma gamma (eta)
# ==========================================================================
alg_dp = TagAnalysis.new("DptoPiPi0EtaDTag")
alg_dp.set_header(["DptoPiPi0EtaDTagAlg/DptoPiPi0EtaDTag.h"])
       .set_constant({ "ECMS" => [:double, 3.773] })
       .note(:bdtd_eta_cut, "BDTG classifier applied on eta candidate (M(gammagamma), chi2(eta), photon helicity angle) for background suppression; BDTG cut retains 77% signal and rejects 84% background for D+ channel")

alg_dp.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi, :DptoKsPi0, :DptoKsPiPiPi, :DptoKKPi
end

alg_dp.signal_side do |s|
  s.charged(pip: 1)
  s.photons 4..48
  s.min_photon_angle 10.0
end

alg_dp.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 200
end

alg_dp.with_decay_card(decay_card_dp).apply
alg_dp.execute_on([psi3770_data, psi3770_incMC, exMC_dp])