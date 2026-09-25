# Search for D_s+ → p pbar e+ nu_e
# arXiv:1910.13097v1 — tag-based D_s ST analysis with semileptonic signal
# 3.19 fb⁻¹ at √s = 4.178 GeV

### Dataset description ###
data_4180 = DatasetManager.real_data.find("703_4180")
incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")

### Decay card for signal process ###
decay_card = <<~DECAYCARD
  Decay psi(4260)
  1.0000 D_s+ D_s*- PHSP;
  Enddecay

  Decay D_s*-
  1.0000 gamma D_s- PHSP;
  Enddecay

  Decay D_s-
  1.0000 K+ K- pi- PHSP;
  Enddecay

  Decay D_s+
  1.0000 p+ anti-p- e+ nu_e PHSP;
  Enddecay

  End
DECAYCARD

### Exclusive MC for signal ###
exMC = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_Ds_ppbar_enu"
  config.related_dataset = data_4180
  config.events = 100_000
  config.decay_card = decay_card
  config.cross_section = :default
end

### Tag-based analysis ###
alg = TagAnalysis.new("DsTagPPbarEnu")
alg.set_header(["DsTagPPbarEnuAlg/DsTagPPbarEnu.h"])
   .set_constant({ "ECMS" => [:double, 4.178] })

# Tag side: D_s- reconstructed via three ST modes
alg.tag_side(:Ds) do |t|
  t.modes :DstoKKPi, :DstoKsK, :DstoKsKminusPiPi
end

# Signal side: D_s+ → p pbar e+ nu_e
alg.signal_side do |s|
  s.charged(prp: 1, prm: 1, ep: 1)
  s.require_charge 1       # +1 -1 +1 = +1 (D_s+ charge)
  s.missing :nu_e           # massless neutrino (semileptonic)
end

# 4C kinematic fit: tag + signal + nu_e = ecms_lab
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg.note(:electron_reconstruction, "electron momentum typically very low; in ~95% of events the electron is not reconstructed. The paper classifies events as signal based on p pbar pair presence even without electron track. Electron momentum required < 0.09 GeV/c when a third track is found.")
   .note(:mm2_requirement, "missing mass squared M M^2 > 0 GeV^2/c^4 required to reduce continuum q qbar background")
   .note(:tag_pion_momentum, "pion momenta from D_s- decay required > 0.1 GeV/c to suppress D* → Dπ background")
   .note(:tag_recoil_mass, "recoil mass against D_s- candidate required in [2.06, 2.18] GeV/c²; best candidate chosen by closest recoil mass to D_s*+ nominal mass")
   .note(:tag_Ks_mass, "K_S0 invariant mass required in (0.487, 0.511) GeV/c²; K_S0 decay length > 2× uncertainty")
   .note(:unused_track_limit, "fewer than four unused charged tracks required after tag reconstruction")

alg.with_decay_card(decay_card).apply
alg.execute_on([data_4180, incMC_4180, exMC])