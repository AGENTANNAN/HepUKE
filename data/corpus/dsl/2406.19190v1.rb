# BESIII Analysis: Ds+ -> K0S e+ nu_e branching fraction
# Paper: 2406.19190v1 — Double-tag method at sqrt(s)=4.128-4.226 GeV
# TagAnalysis: ST Ds- tagged via 14 hadronic modes, signal side: K0S e+ + missing nu_e

### Dataset preparation ###
# Energy points 4.128 – 4.226 GeV
data_4130 = DatasetManager.real_data.find("705_4130")     # 4.128 GeV
data_4160 = DatasetManager.real_data.find("705_4160")     # 4.157 GeV
data_4180 = DatasetManager.real_data.find("703_4180")     # 4.178 GeV
data_4190 = DatasetManager.real_data.find("703_4190")     # 4.189 GeV
data_4200 = DatasetManager.real_data.find("703_4200")     # 4.199 GeV
data_4210 = DatasetManager.real_data.find("703_4210")     # 4.209 GeV
data_4220 = DatasetManager.real_data.find("703_4220")     # 4.219 GeV
data_4230 = DatasetManager.real_data.find("703_4230")     # 4.226 GeV

scan_points = [
  data_4130, data_4160, data_4180, data_4190,
  data_4200, data_4210, data_4220, data_4230
]

# Decay card for Ds+ -> K0S e+ nu_e (signal process)
# KKMC + psi(4260) top mother; Ds*+ Ds- production at scan energies
decay_card_DsKeNu = <<~DECAYCARD
  Decay psi(4260)
  1.0 D_s*- D_s+ PHSP;
  Enddecay

  Decay D_s*-
  1.0 gamma D_s- VSP_PWAVE;
  Enddecay

  Decay D_s+
  1.0 K_S0 e+ nu_e ISGW2;
  Enddecay

  Decay K_S0
  1.0 pi+ pi- PHSP;
  Enddecay

  Decay D_s-
  1.0 K+ K- pi- PHSP;
  Enddecay

  End
DECAYCARD

# Exclusive MC for signal
exMC_DsKeNu = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_DsKeNu"
  config.events        = 500_000
  config.decay_card    = decay_card_DsKeNu
  config.cross_section = :default
end

### TagAnalysis — double-tag Ds -> K0S e+ nu_e ###

alg = TagAnalysis.new("DsTagK0eNu")
alg.set_header(["DsTagK0eNuAlg/DsTagK0eNu.h"])
   .set_constant({ "ECMS" => [:double, 4.178] })   # central value; per-run energy from MeasuredEcmsSvc

# Tag side: ST Ds- via 14 hadronic modes
# Ds modes from the DTagAlg channel-name table
alg.tag_side(:Ds) do |t|
  t.modes :DstoKKPi, :DstoKKPiPi0, :DstoPiPiPi,
          :DstoKsK, :DstoKsKPi0, :DstoKPiPi,
          :DstoKsKsPi, :DstoKsKPiPi, :DstoKsKPiPiPi,
          :DstoEtaGamGamPi, :DstoEtapPiPiEtaRho,
          :DstoEtapEtaGamGamPiPiPi, :DstoEtapGamRhoPi,
          :DstoEtaGamGamRho
  t.charm -1               # tag Ds- (anti-Ds+)
end

# Signal side: K0S -> pip pim + e+ + missing nu_e
# Tracks unused by the tag; lepton key :ep drives electron PID
alg.signal_side do |s|
  s.photons 0                # no extra photons needed (transition gamma handled by DTagAlg)
  s.charged(pip: 1, pim: 1, ep: 1)
  s.require_charge 1         # pip(+) + pim(-) + ep(+) = +1 (Ds+ charge)
  s.missing :nu_e            # massless neutrino; triggers q2 auto-storage
end

# Kinematic fit: 4C fit with optional K0S mass constraint
# The transition gamma + tag Ds- enter the fit automatically
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:pip, :pim).constrain_to_nominal_mass_of(:K_S0)   # K0S mass constraint
  f.chi2_cut 200
end

alg.with_decay_card(decay_card_DsKeNu)

# Notes for inexpressible BOSS-side procedures
alg.note(:transition_gamma_kmfit,
  "Two kinematic fit hypotheses tested for transition gamma: " \
  "(1) gamma + ST Ds- -> Ds*-, (2) gamma + semileptonic -> Ds*+. " \
  "Hypothesis with smallest chi2 (and <100) kept; handled in BOSS " \
  "DTagAlg via the Ds* reconstruction step")
  .note(:k0s_mass_window,
  "K0S candidates: |M(pi+pi-) - m_K0S| < 12 MeV/c2; decay length > 2x vertex resolution; " \
  "handled by DTagAlg V0 reconstruction on signal side")
  .note(:k0s_to_pipi_br,
  "K0 -> K0S branching fraction 50% applied; K0S -> pi+pi- PDG BF used for absolute BF calculation")
  .note(:positron_pid,
  "Positron PID: CLe' > 0.001 and CLe'/(CLe'+CLpi'+CLK') > 0.8, combining dE/dx+TOF+EMC; " \
  "signal-side SimplePIDSvc thresholds are fixed DTagAlg defaults")
  .note(:extra_photon_veto,
  "Largest unused shower energy < 0.2 GeV to reject hadronic Ds decays with fake photons; " \
  "applied in BOSS after DTagAlg signal-side reconstruction")
  .note(:k0e_invmass_veto,
  "M(K0 e+) < 1.78 GeV/c2 veto to suppress Ds+ -> K0 K+ background; ROOT-level cut")
  .note(:transition_pi_veto,
  "Pion momentum > 0.1 GeV/c for pions not from K0S/eta/eta' to suppress D* transition pions; " \
  "applied in BOSS-level selection")

alg.apply
alg.execute_on(scan_points + [exMC_DsKeNu].flatten)