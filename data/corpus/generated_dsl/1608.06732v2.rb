# BOSS-part DSL for the tagged analyses Ds+ -> mu+ nu_mu and Ds+ -> tau+ nu_tau (tau+ -> pi+ nu_tau)
# at sqrt(s) = 4.009 GeV.  Both signal modes share the same hadronic single-tag Ds- reconstruction.

### Dataset description ###
ds_data  = DatasetManager.real_data.find("703_4009")     # 4.009 GeV real data (482 pb^-1)
ds_incMC = DatasetManager.inclusive_mc.find("703_4009")  # matching inclusive MC sample

# ---------------------------------------------------------------------------
# Signal decay cards (EvtGen syntax; top mother psi(4260), tag side left hadronic)
# ---------------------------------------------------------------------------
decay_card_mu = <<~DECAYCARD
    Decay psi(4260)
    1.000 D_s+ D_s- PHSP;
    Enddecay

    Decay D_s+
    1.000 mu+ nu_mu PHSP;
    Enddecay

    Decay D_s-
    1.000 K+ K- pi- PHSP;
    Enddecay

    End
DECAYCARD

decay_card_tau = <<~DECAYCARD
    Decay psi(4260)
    1.000 D_s+ D_s- PHSP;
    Enddecay

    Decay D_s+
    1.000 tau+ nu_tau PHSP;
    Enddecay

    Decay tau+
    1.000 pi+ nu_tau PHSP;
    Enddecay

    Decay D_s-
    1.000 K+ K- pi- PHSP;
    Enddecay

    End
DECAYCARD

# 500k-event exclusive MC for each of the two signal modes
exMC_mu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "DsToMuNu_exclusive_mc"
  config.related_dataset = ds_data
  config.events          = 500_000
  config.decay_card      = decay_card_mu
  config.cross_section   = :default
end

exMC_tau = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "DsToTauNu_exclusive_mc"
  config.related_dataset = ds_data
  config.events          = 500_000
  config.decay_card      = decay_card_tau
  config.cross_section   = :default
end

### Tag-based event selection ###
# Nine hadronic single-tag modes reconstructed on the recoiling Ds- side.
tag_modes = [
  :DstoKsK,                   # Ds- -> K_S0 K-
  :DstoKKPi,                  # Ds- -> K+ K- pi-
  :DstoKKPiPi0,               # Ds- -> K+ K- pi- pi0
  :DstoKsKPiPi,               # Ds- -> K_S0 K+ pi- pi-
  :DstoPiPiPi,                # Ds- -> pi+ pi- pi-
  :DstoEtaPi,                 # Ds- -> eta pi-, eta -> gamma gamma
  :DstoPi0EtaPi,              # Ds- -> pi0 eta pi-
  :DstoEtaPrimePiToPiPiEta,   # Ds- -> eta' pi-, eta' -> pi+ pi- eta
  :DstoEtaPrimePiToPiPiGamma  # Ds- -> eta' pi-, eta' -> pi+ pi- gamma
]

# ============================ Muonic mode: Ds+ -> mu+ nu_mu ============================
alg_mu = TagAnalysis.new("DsToMuNu")
alg_mu.set_header(["DsToMuNuAlg/DsToMuNu.h"])
      .set_constant({"ECMS" => [:double, 4.009]})
      .with_decay_card(decay_card_mu)

# Single tag: hadronic Ds- modes (one tag_side call => ST)
alg_mu.tag_side(:Ds) do |t|
  t.modes(*tag_modes)
  t.charm  -1                              # tag the Ds- (recoiling) side
  t.window :mBC, min: 1.962, max: 1.982    # tag Ds- mass window
end

# Signal side: one charged track opposite the tag + the undetected neutrino, no extra shower
alg_mu.signal_side do |s|
  s.charged(mup: 1)          # single signal-side track (mu+ mass hypothesis)
  s.require_charge(1)        # charge opposite the Ds- tag
  s.photons 0..0             # no unassociated neutral shower on the signal side
  s.min_photon_angle 10.0
  s.min_photon_energy 0.300  # 300 MeV shower-veto threshold
  s.missing :nu_mu           # undetected (massless) neutrino
end

# Four-momentum conservation only (no mass-constrained fit because of the missing neutrino)
alg_mu.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_mu
  .note(:tag_track_selection, "tag-side charged tracks: |cos(theta)|<0.93, |Vz|<10 cm, Vr<1 cm (applied inside DTagAlg)")
  .note(:tag_pid, "tag-side PID: pion if CL_pi>0 and CL_pi>CL_K; kaon if CL_K>0 and CL_K>CL_pi (DTagAlg internal)")
  .note(:tag_photon_selection, "photons: E>25 MeV in the barrel (|cos(theta)|<0.80) or E>50 MeV in the endcap (0.86<|cos(theta)|<0.92), TDC in [0,700] ns, angle >10 deg to any charged track")
  .note(:tag_mass_windows, "pi0: 0.115<M(gg)<0.150; eta: 0.510<M(gg)<0.570; eta'(pi+pi-eta): 0.943<M(pi+pi-eta)<0.973; eta'(gamma rho0): 0.932<M(g rho0)<0.980 with 0.570<M(pi+pi-)<0.970; K_S0: pi+pi- pairs with |Vz|<20 cm, 0.487<M(pi+pi-)<0.511 and secondary-vertex decay length >= 2x its resolution; pi0/eta mass-constrained (1C)")
  .note(:signal_side_track, "single signal-side track treated as mu+ WITHOUT PID - only the muon mass hypothesis is assigned for the missing-mass calculation")
  .note(:missing_mass, "MM2 = (E_beam - E_mu)^2/c^4 - (-p_Ds- - p_mu)^2/c^2 stored as a variable; signal region -0.15 < MM2 < 0.20 (GeV/c^2)^2 applied downstream in ROOT")

alg_mu.apply
alg_mu.execute_on([ds_data, ds_incMC, exMC_mu])

# ====================== Tauonic mode: Ds+ -> tau+ nu_tau, tau+ -> pi+ nu_tau ======================
alg_tau = TagAnalysis.new("DsToTauNu")
alg_tau.set_header(["DsToTauNuAlg/DsToTauNu.h"])
       .set_constant({"ECMS" => [:double, 4.009]})
       .with_decay_card(decay_card_tau)

# Same nine hadronic single-tag modes
alg_tau.tag_side(:Ds) do |t|
  t.modes(*tag_modes)
  t.charm  -1
  t.window :mBC, min: 1.962, max: 1.982
end

# Signal side: one charged track opposite the tag (pi+ mass hypothesis) + undetected neutrinos
alg_tau.signal_side do |s|
  s.charged(pip: 1)          # single signal-side track (pi+ mass hypothesis)
  s.require_charge(1)        # charge opposite the Ds- tag
  s.photons 0..0             # no unassociated neutral shower on the signal side
  s.min_photon_angle 10.0
  s.min_photon_energy 0.300  # 300 MeV shower-veto threshold
  s.missing :nu_tau          # undetected (massless) neutrinos from Ds+ and tau+ decays, summed
end

# Four-momentum conservation only (no mass-constrained fit because of the missing neutrinos)
alg_tau.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_tau
  .note(:tag_track_selection, "tag-side charged tracks: |cos(theta)|<0.93, |Vz|<10 cm, Vr<1 cm (applied inside DTagAlg)")
  .note(:tag_pid, "tag-side PID: pion if CL_pi>0 and CL_pi>CL_K; kaon if CL_K>0 and CL_K>CL_pi (DTagAlg internal)")
  .note(:tag_photon_selection, "photons: E>25 MeV in the barrel (|cos(theta)|<0.80) or E>50 MeV in the endcap (0.86<|cos(theta)|<0.92), TDC in [0,700] ns, angle >10 deg to any charged track")
  .note(:tag_mass_windows, "pi0: 0.115<M(gg)<0.150; eta: 0.510<M(gg)<0.570; eta'(pi+pi-eta): 0.943<M(pi+pi-eta)<0.973; eta'(gamma rho0): 0.932<M(g rho0)<0.980 with 0.570<M(pi+pi-)<0.970; K_S0: pi+pi- pairs with |Vz|<20 cm, 0.487<M(pi+pi-)<0.511 and secondary-vertex decay length >= 2x its resolution; pi0/eta mass-constrained (1C)")
  .note(:signal_side_track, "single signal-side track treated as pi+ WITHOUT PID - only the pion mass hypothesis is assigned for the missing-mass calculation")
  .note(:missing_mass, "MM2 = (E_beam - E_pi)^2/c^4 - (-p_Ds- - p_pi)^2/c^2 stored as a variable; signal region -0.15 < MM2 < 0.20 (GeV/c^2)^2 applied downstream in ROOT")
  .note(:efficiency_curve, "detection efficiency referenced from the description: 91.4 +/- 0.5 % (mu nu) and 41.0 +/- 0.3 % (tau nu, tau -> pi nu)")

alg_tau.apply
alg_tau.execute_on([ds_data, ds_incMC, exMC_tau])