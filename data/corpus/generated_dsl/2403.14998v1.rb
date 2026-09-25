# =====================================================================
# e+e- -> Ds+ Ds- Born cross-section scan (3.94 - 4.95 GeV)
# Single-tag method: tag Ds- -> K+ K- pi-, infer the recoil Ds+ via
# its recoil mass.  Tag-based layer (TagAnalysis), not Algorithm+Selection.
# =====================================================================

### Dataset preparation ###
# Full real-data scan: 138 energy points inside [3.94, 4.95] GeV (~22.9 fb^-1)
scan_data  = DatasetManager.find_by_energy(3940, 4950)
# Matching inclusive MC at the same scan points
scan_incMC = DatasetManager.inclusive_mc.where(cms_energy: { value: 3940..4950 })

# Signal decay card (EvtGen): psi(4260) -> Ds+ Ds-, Ds -> K+ K- pi
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000  Ds+  Ds-        PHSP;
    Enddecay

    Decay Ds+
    1.0000  K+   K-   pi+   PHSP;
    Enddecay

    Decay Ds-
    1.0000  K+   K-   pi-   PHSP;
    Enddecay

    End
DECAYCARD

# Peaking-background decay card: psi(4260) -> Ds+ Ds*-, Ds*- -> Ds- gamma
decay_card_bkg = <<~DECAYCARD
    Decay psi(4260)
    1.0000  Ds+  Ds*-       PHSP;
    Enddecay

    Decay Ds*-
    1.0000  Ds-  gamma      PHSP;
    Enddecay

    Decay Ds-
    1.0000  K+   K-   pi-   PHSP;
    Enddecay

    Decay Ds+
    1.0000  K+   K-   pi+   PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event signal MC at every scan point (flat Born cross section as first input)
exMC_signal = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name   = "DsDs_signal_mc"   # auto-suffixed per scan point
  config.events        = 100_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

# Peaking-background (Ds+ Ds*-) MC, needed for the Ecms > 4.6 GeV points
exMC_bkg = DatasetManager.create_exclusive_mc_for(DatasetManager.find_by_energy(4600, 4950)) do |config|
  config.sample_name   = "DsDsstar_bkg_mc"
  config.events        = 100_000
  config.decay_card    = decay_card_bkg
  config.cross_section = :default
end

### Event selection (BOSS) — single tag ###
alg = TagAnalysis.new("DsST")
alg.set_header(["DsSTAlg/DsST.h"])
   .set_constant({ "ECMS" => [:double, 4.260] })   # representative; per-run energy comes from MeasuredEcmsSvc
   .with_decay_card(decay_card_signal)
   # inexpressible / not representable in the store-not-cut tag layer:
   .note(:recoil_mass_window, "recoil-mass window 1.945 < RM(Ds-) < 1.990 GeV/c^2 applied on the " \
        "mass-constrained combination M_recoil(K+K-pi-) + M(K+K-pi-) - m(Ds-); missing mass is stored " \
        "unconditionally at BOSS level and the window is imposed in ROOT")
   .note(:resonance_windows, "tag sub-mode resonance windows inside DTagAlg: 1.005 < M(K+K-) < 1.035 GeV/c^2 " \
        "for phi, and 0.832 < M(K+pi-) < 0.928 GeV/c^2 with helicity angle |cos theta(K+)| < 0.52 for K*(892)0")
   .note(:background_veto, "peaking background e+e- -> Ds+ Ds*- subtracted for Ecms > 4.6 GeV, estimated from " \
        "exclusive MC normalised to luminosity and cross sections")

# Tag side: Ds- reconstructed from the pre-stored DTag candidates (single tag),
# through the two intermediate modes phi pi- and K*0 K-
alg.tag_side(:Ds) do |t|
  t.modes :DstoPhiPi, :DstoKstarK   # Ds- -> phi pi- (phi -> K+K-), Ds- -> K*0 K- (K*0 -> K+ pi-)
  t.charm -1
end

# Signal side: the recoiling Ds+ is not reconstructed -> inferred from the recoil (missing) mass
alg.signal_side do |s|
  s.missing :Ds
end

# Kinematic fit: 4-momentum conservation against the measured CMS 4-vector, chi2 < 200
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg.apply    # tag spec: takes no Selection argument
alg.execute_on(scan_data + scan_incMC + exMC_signal + exMC_bkg)