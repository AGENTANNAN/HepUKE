### Dataset description ###
# Three data-set groups / six energy points: 4.178, 4.189-4.219, 4.225-4.230 GeV
data_4178 = DatasetManager.real_data.find("703_4180")   # 4.178 GeV
data_4190 = DatasetManager.real_data.find("703_4190")   # 4.1888 GeV
data_4200 = DatasetManager.real_data.find("703_4200")   # 4.1989 GeV
data_4210 = DatasetManager.real_data.find("703_4210")   # 4.2092 GeV
data_4220 = DatasetManager.real_data.find("703_4220")   # 4.2187 GeV
data_4230 = DatasetManager.real_data.find("703_4230")   # 4.2263 GeV
data_points = [data_4178, data_4190, data_4200, data_4210, data_4220, data_4230]

incMC_4178 = DatasetManager.inclusive_mc.find("703_4180")
incMC_4190 = DatasetManager.inclusive_mc.find("703_4190")
incMC_4200 = DatasetManager.inclusive_mc.find("703_4200")
incMC_4210 = DatasetManager.inclusive_mc.find("703_4210")
incMC_4220 = DatasetManager.inclusive_mc.find("703_4220")
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")
incMC_points = [incMC_4178, incMC_4190, incMC_4200, incMC_4210, incMC_4220, incMC_4230]

# Decay card (EvtGen syntax) for e+e- -> D_s+ D_s*-
#   signal : D_s+ -> phi e+ nu_e  (phi -> K+ K-)
#   tag    : D_s- -> K+ K- pi-
decay_card_signal = <<~DECAYCARD
  Decay psi(4260)
  1.0000 D_s+ D_s*-  PHSP;
  Enddecay

  Decay D_s*-
  1.0000 D_s- gamma  PHSP;
  Enddecay

  Decay D_s+
  1.0000 phi e+ nu_e  PHSP;
  Enddecay

  Decay phi
  1.0000 K+ K-  VSS;
  Enddecay

  Decay D_s-
  1.0000 K+ K- pi-  PHSP;
  Enddecay

  End
DECAYCARD

# 1M-event exclusive signal MC, one sample per energy point
exMCs_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_Ds_Xenu_phi"
  config.events        = 1_000_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Tag-based event selection (BOSS) ###
alg_name = "DsTagSemiElec"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({ "ECMS" => [:double, 4.178] })

# Tag side: single tag D_s- -> K+ K- pi- (charm -1)
alg.tag_side(:Ds) do |t|
  t.modes :DstoKKPi
  t.charm -1
end

# Signal side of the single tag: D_s+ -> X e+ nu_e
alg.signal_side do |s|
  s.photons 0                       # require zero photons on the signal side
  s.charged(ep: 1, at_least: true)  # exactly one positron (X hadrons allowed as extra charged)
  s.missing :nu_e                   # missing neutrino (massless form)
end

# 4C kinematic fit against the measured CMS four-momentum
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

# The recoil-mass windows, the unbinned ST-yield fit, the 18-bin positron momentum
# binning with 3x3 PID unfolding, RS-WS subtraction, tau nu subtraction, the
# extrapolation with six exclusive semileptonic modes and the BF extraction
# (N_DT / (n_ST * b_tag)) are all post-fit ROOT-level procedures.
alg.note(:positron_momentum_cut,
         "signal-side positron required to have p > 200 MeV/c; the signal_side charged " \
         "multiset carries no momentum requirement, so this cut is applied as an event " \
         "selection at BOSS level")

alg.apply
alg.execute_on(data_points + incMC_points + exMCs_signal)