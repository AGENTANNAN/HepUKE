### Dataset description ###
# sqrt(s) = 4.128-4.226 GeV: seven energy points, 7.33 fb^-1 in total
ds_data_points = [
  DatasetManager.real_data.find("705_4130"),   # 4.128 GeV
  DatasetManager.real_data.find("705_4160"),   # 4.158 GeV
  DatasetManager.real_data.find("703_4180"),   # 4.178 GeV
  DatasetManager.real_data.find("703_4190"),   # 4.189 GeV
  DatasetManager.real_data.find("703_4200"),   # 4.199 GeV
  DatasetManager.real_data.find("703_4220"),   # 4.219 GeV
  DatasetManager.real_data.find("703_4230")    # 4.226 GeV
]

ds_incMC_points = [
  DatasetManager.inclusive_mc.find("705_4130"),
  DatasetManager.inclusive_mc.find("705_4160"),
  DatasetManager.inclusive_mc.find("703_4180"),
  DatasetManager.inclusive_mc.find("703_4190"),
  DatasetManager.inclusive_mc.find("703_4200"),
  DatasetManager.inclusive_mc.find("703_4220"),
  DatasetManager.inclusive_mc.find("703_4230")
]

# Decay card for the signal process e+e- -> Ds*+ Ds-, Ds+ -> K+K-pi+, Ds- -> K+K-pi-
# (EvtGen format, EvtGen particle names; the Ds* -> Ds gamma transition is not reconstructed)
decay_card_signal = <<~DECAYCARD
  Decay psi(4260)
  1.0000 Ds*+ Ds- PHSP;
  Enddecay

  Decay Ds*+
  1.0000 Ds+ gamma VSS;
  Enddecay

  Decay Ds+
  1.0000 K+ K- pi+ PHSP;
  Enddecay

  Decay Ds-
  1.0000 K+ K- pi- PHSP;
  Enddecay

  End
DECAYCARD

# 200k-event exclusive signal MC, one sample per energy point of the scan
exMC_signal = DatasetManager.create_exclusive_mc_for(ds_data_points) do |config|
  config.sample_name   = "exmc_dsstar_ds_kkpi"
  config.events        = 200000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) — tag-based double-tag analysis ###
alg_name = "DsTagKKPi"
ds_tag_alg = TagAnalysis.new(alg_name)
ds_tag_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
          .set_constant({"ECMS" => [:double, 4.178]})  # representative point of the 4.128-4.226 GeV scan

# BOSS-side procedures with no dedicated DSL construct
ds_tag_alg
  .note(:ds_star_transition, "the Ds*+ -> Ds+ gamma / pi0 transition is not reconstructed; its photon / pi0 remains unreconstructed in the event")
  .note(:best_candidate_selection, "double-tag best candidate chosen by minimum |M_rec - m(Ds*+)| from the Ds recoil mass")
  .note(:mrec_window, "Ds recoil mass window [2.10, 2.13] GeV/c^2 applied to the pi+pi+pi-eta, pi+pi0eta'(gamma rho) and K_S pi+pi0 tag modes")
  .note(:photon_selection, "good photons: E > 25 MeV (barrel) / 50 MeV (endcap), angle to nearest charged track > 10 deg, 0 < t < 700 ns")
  .note(:ks_selection, "K_S0 -> pi+pi-: |M(pi+pi-) - m(K_S)| < 12 MeV/c^2, secondary vertex chi2 < 100, decay length > 2 sigma")
  .note(:pi0_eta_1c_fit, "pi0/eta from gamma gamma use a 1C mass-constrained fit (chi2 < 30) with masses in [115,150] and [490,580] MeV/c^2 and at least one barrel photon")
  .note(:resonance_windows, "eta->pi+pi-pi0 [530,560], eta'->pi+pi-eta [943,973], eta'->gamma rho [946,970], rho0->pi+pi- [570,970] MeV/c^2")
  .note(:pip_pim_veto, "K+pi+pi- mode vetoes pi+pi- pairs with invariant mass in [487,511] MeV/c^2")
  .note(:dt_yield_extraction, "double-tag yields counted in |m_bar - m(Ds)| < 15 MeV/c^2 and |Delta_m| < 30 MeV/c^2, with sideband subtraction over 80 < |Delta_m| < 140 MeV/c^2")

# Tag BOTH charm sides with the hadronic tag modes -> double tag
ds_tag_alg.tag_side(:Ds) do |t|
  t.modes :DstoKKPi, :DstoKsK, :DstoKsKPi0, :DstoKsKsPi, :DstoKKPiPi0, :DstoKsKPiPi
end

ds_tag_alg.tag_side(:Ds) do |t|
  t.modes :DstoKKPi, :DstoKsK, :DstoKsKPi0, :DstoKsKsPi, :DstoKKPiPi0, :DstoKsKPiPi
end

# Signal side (everything the tag did not use) requires zero good photons
ds_tag_alg.signal_side do |s|
  s.photons 0
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025   # EMC barrel floor; endcap floor 0.050 handled downstream
end

# Four-momentum-constrained kinematic fit
ds_tag_alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

ds_tag_alg.with_decay_card(decay_card_signal).apply

# Execute on the scan data, inclusive MC and the signal MC samples
root_files = ds_tag_alg.execute_on(ds_data_points + ds_incMC_points + exMC_signal)