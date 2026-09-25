require 'dsl4bes'

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# ============================================================
# Channel 1: J/psi -> gamma eta' -> pi+ pi- eta -> gamma mu+ mu-
# 5C kinematic fit: 4-momentum conservation + eta mass constraint
# The full 6C fit (adding eta' mass constraint) is applied in ROOT
# ============================================================

decay_card_mumu = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta_prime HELAMP;
  Enddecay
  Decay eta_prime
  1.0000 pi+ pi- eta D_DALITZ;
  Enddecay
  Decay eta
  1.0000 gamma mu+ mu- TFF;
  Enddecay
  End
DECAYCARD

exMC_mumu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_eta_tff_mumu"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_mumu
  config.cross_section   = :default
end

algo_mumu = Algorithm.new("EtaTFFMuMu")
algo_mumu
  .set_header(["EvtRecEvent/EvtRecTrack.h", "EventModel/Event.h", "EtaTFFMuMu/EtaTFFMuMu.h"])
  .set_constant(ECMS: 3.097)

sel_mumu = Selection.new
  .select_track do
    cos_theta   0.93
    Vz          10.0
    Vr          1.0
    nChrp       "==2"
    nChrn       "==2"
  end
  .select_photon do
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=2"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :pion, against: [:muon, :electron]
    identify :muon, against: [:pion, :electron]
    npip "==1"; npim "==1"
    nmup "==1"; nmum "==1"
  end
  .kinematic_fit([:gamma, :gamma, :pip, :pim, :mup, :mum]) do
    invariant_mass_of(:gamma, :mup, :mum).constrain_to_nominal_mass_of(:eta)
    constrain_four_momentum
    chi2_cut 40
    nominal
  end

algo_mumu
  .note(:full_6c_kinematic_fit, "Full 6C kinematic fit adds eta' mass constraint M(pip, pim, gamma_eta, mu+, mu-) = m(eta') (1C) on top of the 5C (4C + eta mass). The eta' constraint and final chi2 optimization are applied in the ROOT analysis stage.")
  .with_decay_card(decay_card_mumu)
  .apply(sel_mumu)

algo_mumu.execute_on([jpsi_data, jpsi_incMC, exMC_mumu])

# ============================================================
# Channel 2: J/psi -> gamma eta' -> pi+ pi- eta -> gamma e+ e-
# 5C kinematic fit: 4-momentum conservation + eta mass constraint
# The full 6C fit (adding eta' mass constraint) is applied in ROOT
# Photon conversion veto (PCF) applied in ROOT analysis
# ============================================================

decay_card_ee = <<~DECAYCARD
  Decay J/psi
  1.0000 gamma eta_prime HELAMP;
  Enddecay
  Decay eta_prime
  1.0000 pi+ pi- eta D_DALITZ;
  Enddecay
  Decay eta
  1.0000 gamma e+ e- TFF;
  Enddecay
  End
DECAYCARD

exMC_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_eta_tff_ee"
  config.related_dataset = jpsi_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_ee
  config.cross_section   = :default
end

algo_ee = Algorithm.new("EtaTFFEE")
algo_ee
  .set_header(["EvtRecEvent/EvtRecTrack.h", "EventModel/Event.h", "EtaTFFEE/EtaTFFEE.h"])
  .set_constant(ECMS: 3.097)

sel_ee = Selection.new
  .select_track do
    cos_theta   0.93
    Vz          10.0
    Vr          1.0
    nChrp       "==2"
    nChrn       "==2"
  end
  .select_photon do
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=2"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :pion, against: [:electron]
    identify :electron, against: [:pion]
    npip "==1"; npim "==1"
    nep "==1"; nem "==1"
  end
  .kinematic_fit([:gamma, :gamma, :pip, :pim, :ep, :em]) do
    invariant_mass_of(:gamma, :ep, :em).constrain_to_nominal_mass_of(:eta)
    constrain_four_momentum
    chi2_cut 200
    nominal
  end

algo_ee
  .note(:full_6c_kinematic_fit, "Full 6C kinematic fit adds eta' mass constraint M(pip, pim, gamma_eta, e+, e-) = m(eta') (1C) on top of the 5C (4C + eta mass). The eta' constraint and final chi2 optimization are applied in the ROOT analysis stage.")
  .note(:photon_conversion_veto, "Photon conversion veto using PCF: events with 2.0 cm < R_xy < 8.0 cm and cos(theta_eg) > 0.5 are rejected to suppress gamma -> e+e- backgrounds. Applied in ROOT analysis.")
  .with_decay_card(decay_card_ee)
  .apply(sel_ee)

algo_ee.execute_on([jpsi_data, jpsi_incMC, exMC_ee])