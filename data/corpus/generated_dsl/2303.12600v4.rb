### Dataset preparation ###
# Real data and inclusive MC at the 8 c.m. energy points spanning 4.128-4.226 GeV
data_4130 = DatasetManager.real_data.find("705_4130")        # 4.128 GeV
incMC_4130 = DatasetManager.inclusive_mc.find("705_4130")
data_4160 = DatasetManager.real_data.find("705_4160")        # 4.157 GeV
incMC_4160 = DatasetManager.inclusive_mc.find("705_4160")
data_4180 = DatasetManager.real_data.find("703_4180")        # 4.178 GeV
incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")
data_4190 = DatasetManager.real_data.find("703_4190")        # 4.189 GeV
incMC_4190 = DatasetManager.inclusive_mc.find("703_4190")
data_4200 = DatasetManager.real_data.find("703_4200")        # 4.199 GeV
incMC_4200 = DatasetManager.inclusive_mc.find("703_4200")
data_4210 = DatasetManager.real_data.find("703_4210")        # 4.209 GeV
incMC_4210 = DatasetManager.inclusive_mc.find("703_4210")
data_4220 = DatasetManager.real_data.find("703_4220")        # 4.219 GeV
incMC_4220 = DatasetManager.inclusive_mc.find("703_4220")
data_4230 = DatasetManager.real_data.find("703_4230")        # 4.226 GeV
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")

# Decay card for the signal process (EvtGen format, EvtGen particle names).
# e+e- -> D_s*+ D_s- (no explicit top-meson given -> BESIII KKMC convention psi(4260)),
# D_s*+ -> gamma D_s+, D_s+ -> tau+ nu_tau, tau+ -> pi+ anti-nu_tau, D_s- -> K+ K- pi-
decay_card_signal = <<~DECAYCARD
  Decay psi(4260)
  1.000 D_s*+ D_s-
  Enddecay

  Decay D_s*+
  1.000 gamma D_s+
  Enddecay

  Decay D_s+
  1.000 tau+ anti-nu_tau
  Enddecay

  Decay tau+
  1.000 pi+ anti-nu_tau
  Enddecay

  Decay D_s-
  1.000 K+ K- pi-
  Enddecay

  End
DECAYCARD

# 500k-event exclusive signal MC for D_s*+ D_s- production
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_dsstar_ds_taunu"
  config.related_dataset = data_4180        # associated real dataset near the D_s* threshold
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (tag-based, BOSS) ###
alg = TagAnalysis.new("DsTagTauNu")
alg.set_header(["DsTagTauNuAlg/DsTagTauNu.h"])
   .set_constant({"ECMS" => [:double, 4.178]})

# Tag side: the D_s- reconstructed in the 9 implemented hadronic tag modes, charm = -1
alg.tag_side(:Ds) do |t|
  t.modes :DstoKKPi,               # D_s- -> K+ K- pi-
          :DstoKKPiPi0,            # D_s- -> K+ K- pi- pi0
          :DstoPiPiPi,             # D_s- -> pi+ pi- pi-
          :DstoKsK,                # D_s- -> K_S0 K-
          :DstoKsKPiPi,            # D_s- -> K_S0 K+ pi- pi-
          :DstoPiEta,              # D_s- -> pi- eta
          :DstoPiPi0Eta,           # D_s- -> pi- pi0 eta
          :DstoPiEtaPrimeToPiPiEta,# D_s- -> pi- eta', eta' -> pi+ pi- eta
          :DstoPiEtaPrimeToGammaRho# D_s- -> pi- eta', eta' -> gamma rho0
  t.charm -1
end

# Signal side (what the tag did not use): the D_s*+ -> gamma D_s+ transition photon,
# exactly one pi+ of positive charge (opposite the charm -1 tag), and a missing
# massless neutrino system (nu_tau + anti-nu_tau from D_s+ -> tau+ nu_tau, tau+ -> pi+ anti-nu_tau).
# `charged(pip: 1)` also enforces "no extra charged tracks" and applies pion PID.
alg.signal_side do |s|
  s.photons 1
  s.charged(pip: 1)
  s.require_charge(1)
  s.missing :nu
end

# 4C kinematic fit over the derived participants (tag + signal photon + pi+ + missing nu)
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

# BOSS-side procedures that the DSL cannot express
alg.note(:transition_photon_selection, "the D_s*+ -> gamma D_s+ transition photon is retained with energy 0.114-0.149 GeV in the D_s* rest frame; among candidate combinations the one minimizing |DeltaE| (direct/indirect tag) and giving the D_s* mass closest to nominal is kept, and finally the candidate whose recoil mass against the tagged D_s- is closest to the nominal D_s* mass is selected")
   .note(:background_veto, "further double-tag selection applied on the signal side: no extra charged track and no extra pi0; E/p < 0.9 for the signal pion; maximum extra photon energy < 0.3 GeV; |cos(theta_miss)| < 0.9; M_miss^2 within [-0.2, 0.6] GeV^2/c^4")

alg.with_decay_card(decay_card_signal).apply
alg.execute_on([data_4130, incMC_4130,
                data_4160, incMC_4160,
                data_4180, incMC_4180,
                data_4190, incMC_4190,
                data_4200, incMC_4200,
                data_4210, incMC_4210,
                data_4220, incMC_4220,
                data_4230, incMC_4230,
                exMC_signal])