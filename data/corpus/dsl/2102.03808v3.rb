DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

data_4178 = DatasetManager.real_data.find("703_4180")
data_4189 = DatasetManager.real_data.find("703_4190")
data_4199 = DatasetManager.real_data.find("703_4200")
data_4209 = DatasetManager.real_data.find("703_4210")
data_4219 = DatasetManager.real_data.find("703_4220")
data_4226 = DatasetManager.real_data.find("703_4230")

incMC_4178 = DatasetManager.inclusive_mc.find("703_4180")
incMC_4189 = DatasetManager.inclusive_mc.find("703_4190")
incMC_4199 = DatasetManager.inclusive_mc.find("703_4200")
incMC_4209 = DatasetManager.inclusive_mc.find("703_4210")
incMC_4219 = DatasetManager.inclusive_mc.find("703_4220")
incMC_4226 = DatasetManager.inclusive_mc.find("703_4230")

datasets     = [data_4178, data_4189, data_4199, data_4209, data_4219, data_4226]
incMC_datasets = [incMC_4178, incMC_4189, incMC_4199, incMC_4209, incMC_4219, incMC_4226]

decay_card = <<~DECAYCARD
  Decay psi(4260)
  1.000 D_s*+ D_s- PHSP;
  Enddecay

  Decay D_s*+
  1.000 gamma D_s+ PHSP;
  Enddecay

  Decay D_s+
  1.000 K_S0 K- pi+ pi+ PHSP;
  Enddecay

  Decay D_s-
  1.000 K+ K- pi- PHSP;
  Enddecay

  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_Ds_KsKpipi"
  config.related_dataset = data_4178
  config.events          = 1_000_000
  config.decay_card      = decay_card
  config.cross_section   = :default
end

# TagAnalysis for Ds double-tag: Ds- ST + Ds+ → K_S0 K- π+ π+
alg = TagAnalysis.new("DsDTagKsKpipi")
alg.set_header(["DsDTagKsKpipiAlg/DsDTagKsKpipi.h"])
   .note(:tag_mode_unavailable, "Ds- → π- η' (DstoPiEtap) mode not available in authorized tag-mode list;
     only 8 of 9 ST modes used")
   .note(:tag_mode_unavailable, "Ds- → K_S0 K- π+ π- mode not available; closest is :DstoKsKminusPiPi")
   .note(:background_veto, "Background from D0 → K- π+ π+ π- vs D0bar → K_S0 K+ K- (K_S0 π+ π-) vetoed
     via D0/D0bar invariant-mass checks")
   .note(:helix_correction, "Helix parameter correction applied before kinematic fit")

alg.tag_side(:Ds_minus) {
  modes :DstoKsK, :DstoKKPi, :DstoKsKminusPiPi, :DstoKPiPi, :DstoPiEta,
        :DstoKsKPi0, :DstoKKPiPiPi, :DstoPiPiPi
  charm -1
}

alg.signal_side {
  charged(km: 1, pip: 2)
  require_charge 1
  min_photon_energy 0.025
}

alg.fit {
  constrain_four_momentum
  chi2_cut 200
}

alg.with_decay_card(decay_card).apply
alg.execute_on(datasets + incMC_datasets + [exMC_signal])