# Paper: 1811.00752v3 — Observation of D_s+ → p n_bar
# BESIII, √s = 4.178 GeV, 3.19 fb⁻¹
# Tag-based analysis: D_s ST + missing (anti-neutron) + isolated photon from D_s*+

### Dataset preparation ###
data_4180 = DatasetManager.real_data.find("703_4180")
incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")

decay_card = <<~DECAYCARD
  Decay psi(4260)
  1.000 D_s*+ D_s- PHSP;
  Enddecay

  Decay D_s*+
  1.000 D_s+ gamma VSP_PWAVE;
  Enddecay

  Decay D_s+
  1.000 p+ anti-n- PHSP;
  Enddecay

  End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Ds_p_nbar"
  config.related_dataset = data_4180
  config.events = 4_000_000
  config.decay_card = decay_card
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg = TagAnalysis.new("DspNbar")
alg.set_header(["DspNbarAlg/DspNbar.h"])
  .set_constant({ "ECMS" => [:double, 4.178] })
  .with_decay_card(decay_card)
  .note(:background_veto, "kinematic fit with constraints on D_s-, D_s+, D_s* masses and initial four-momentum; anti-neutron treated as missing particle; fit chi2 < 200 selects best photon hypothesis")

# 11 ST tag modes for D_s- reconstruction
alg.tag_side(:Ds) do |t|
  t.modes :DstoKsK, :DstoKKPi, :DstoKsKPi0, :DstoKKPiPi0,
          :DstoKsKminusPiPi, :DstoPiPiPi, :DstoPiEta,
          :DstoPiPi0Eta, :DstoPiEPPiPiEta, :DstoPiEPRhoGam,
          :DstoKPiPi
end

alg.signal_side do |s|
  s.photons 1
  s.min_photon_angle 10.0
  s.charged(prp: 1)
  s.require_charge 1
  s.missing :n_bar
end

alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg.apply
alg.execute_on([data_4180, incMC_4180, exMC_signal])