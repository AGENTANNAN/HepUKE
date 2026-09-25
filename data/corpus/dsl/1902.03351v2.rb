# 1902.03351v2: Search for Ds+ -> gamma e+ nu_e at 4.178 GeV
# Modified double-tag technique; ST on Ds- then signal Ds+ -> gamma e+ nu_e

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

data_4178 = DatasetManager.real_data.find("703_4180")
incMC_4178 = DatasetManager.inclusive_mc.find("703_4180")

decay_card = <<~DECAYCARD
  Decay psi(4260)
  1.000 D_s+ D_s*- PHSP;
  Enddecay

  Decay D_s*-
  1.000 D_s- gamma PHSP;
  Enddecay

  Decay D_s+
  1.000 gamma e+ nu_e PHSP;
  Enddecay

  Decay D_s-
  1.000 K+ K- pi- PHSP;
  Enddecay

  End
DECAYCARD

exMC = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_DstoGammaENu"
  config.related_dataset = data_4178
  config.events = 100000
  config.decay_card = decay_card
  config.cross_section = :default
end

alg = TagAnalysis.new("DstoGammaENu")
alg.set_header(["DstoGammaENuAlg/DstoGammaENu.h"])
    .set_constant({"ECMS" => [:double, 4.178]})
    .with_decay_card(decay_card)

# Single tag on Ds- with all hadronic modes
alg.tag_side(:Ds) do |t|
  t.modes :DstoKKPi, :DstoKKPiPi0, :DstoKsK,
          :DstoPiEta, :DstoPiPiPi, :DstoKPiPi,
          :DstoKsKminusPiPi, :DstoKsKplusPiPi,
          :DstoKsKsPi, :DstoKsKPi0, :DstoKsPi,
          :DstoPiPiPiEta, :DstoPiPi0Eta
  t.charm -1
end

# Signal side: radiative photon + positron + missing neutrino
alg.signal_side do |s|
  s.photons 1
  s.charged(ep: 1)
  s.require_charge 1
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
  s.missing :nu_e
end

alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

# Some tag modes have non-trivial intermediate states not directly mappable.
# The DTagAlg reconstructs eta/eta'/rho/pion intermediate states from their
# dominant decay channels internally.
alg.note(:tag_modes_eta_submodes,
  "Ds tag modes with eta sub-decays (gamma gamma, pi+pi-pi0) and eta' sub-decays " \
  "(eta_gammagamma pi+pi-, gamma rho0) reconstructed internally by DTagTool; " \
  "exact sub-mode used follows the DTagAlg default reconstruction")
alg.note(:Ds_star_soft_gamma,
  "soft gamma/pi0 from Ds*+- -> Ds+- gamma/pi0 handled by DTagAlg tag-cleaning; " \
  "four signal hypotheses (Ds+ from Ds*+ or direct) resolved in ROOT-level kinfit")
alg.note(:blind_analysis,
  "analysis performed as blind analysis; signal region unblinded after " \
  "selection optimization on inclusive MC")

alg.apply
alg.execute_on([data_4178, incMC_4178, exMC])