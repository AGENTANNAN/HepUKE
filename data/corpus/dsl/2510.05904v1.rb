# BESIII TagAnalysis: First measurement of Ds+ -> K0 mu+ nu_mu
# Paper: arXiv:2510.05904v1
# Dataset: e+e- at 4.128-4.226 GeV, 7.33 fb^-1
# Tag pattern: ST+missing -- Ds- tagged via 14 hadronic modes,
#   signal side: Ds+ -> K0(->K_S0->pi+pi-) mu+ nu_mu
#
# Per Rule T3: loose chi2_cut 200 in BOSS; the paper's chi2 < 40 is ROOT-level.

dm = DatasetManager.new
dm.real_data.find("703_4180")
dm.inclusive_mc.find("703_4180")

tag = TagAnalysis.new("Ds_K0MuNu")
tag.set_header("First measurement of Ds+ -> K0 mu+ nu_mu",
  journal: "arXiv:2510.05904v1",
  dataset: "e+e- at 4.128-4.226 GeV, 7.33 fb^-1 (703_4180 primary)")
tag.set_constant(:ECMS, 4178.0)

# Ds- tag side: 14 hadronic modes (from Fig 1 of the paper,
# mapped to authoritative DTagAlg channel names)
tag.tag_side(:Ds) do
  modes :DstoKKPi, :DstoKKPiPi0, :DstoKsK, :DstoKsKPi0,
        :DstoKsPi, :DstoKsPiPi0, :DstoKPiPi, :DstoKPiPiPi0,
        :DstoPiPiPi, :DstoPiPiPiPi0, :DstoPiPiPiPiPi,
        :DstoPiEta, :DstoPiPi0, :DstoPiEtaPiPiPi0
end

# Signal side: K0 -> K_S0 -> pi+ pi-, plus mu+ and missing neutrino
tag.signal_side do
  charged pip: 1, pim: 1, mup: 1
  missing :nu_mu
end

tag.fit do
  constrain_four_momentum :Ds
  chi2_cut 200
end

# Muon PID, background rejection, and other cuts applied in ROOT
tag.note(:muon_pid, "Muon: L_mu > 0.001, L_mu > L_e, L_mu > L_K; EMC deposited: 0.1 < E < 0.3 GeV")
tag.note(:K0mu_mass_cut, "M(K0 mu+) < 1.70 GeV/c^2 rejects Ds+ -> K0 pi+ background")
tag.note(:extra_photon_veto, "Max unused photon energy < 0.15 GeV suppresses pi0 backgrounds")
tag.note(:DstarDs_kinematic_fit, "D_s*+- D_s-+ pair constrained; chi2 < 40 per paper, applied in ROOT")
tag.note(:MM2_fit, "Signal DT yield from unbinned ML fit to MM^2 distribution")
tag.note(:multi_energy, "Data at 4.128-4.226 GeV spanning multiple energy points; 703_4180 used as primary")
tag.note(:tag_mode_count, "14 Ds hadronic tag modes from paper Fig 1; check exact mode set against paper")

tag.with_decay_card(%{
  Decay D_s+
  1.0 K0 mu+ nu_mu PHSP;
  CDecay K0
  1.0 K_S0 PHSP;
  EndCDecay
  CDecay K_S0
  1.0 pi+ pi- PHSP;
  EndCDecay
  EndDecay
})

tag.execute_on(["703_4180"])