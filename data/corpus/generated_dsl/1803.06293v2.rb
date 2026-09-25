# ============================================================================
# Dataset preparation — psi(3770), sqrt(s) = 3.773 GeV
# ============================================================================
data_3773  = DatasetManager.real_data.find("712_3773")        # 3.773 GeV real data
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")     # corresponding inclusive MC

# Decay card: e+e- -> D0 anti-D0
decay_card_D0D0bar = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D0 anti-D0 PHSP;
  Enddecay

  Decay D0
  1.0000 K- pi+ PHSP;
  Enddecay

  Decay anti-D0
  1.0000 K+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# Decay card: e+e- -> D+ D-  (pi0 and K_S0 decay rules used by the tag modes)
decay_card_DpDm = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D+ D- PHSP;
  Enddecay

  Decay D+
  1.0000 K- pi+ pi+ PHSP;
  Enddecay

  Decay D-
  1.0000 K+ pi- pi- PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  Decay K_S0
  1.0000 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# 100k-event exclusive MC samples for the two signal processes
exMC_D0D0bar = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0D0bar"
  config.related_dataset = data_3773
  config.events          = 100_000
  config.decay_card      = decay_card_D0D0bar
  config.cross_section   = :default
end

exMC_DpDm = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_DpDm"
  config.related_dataset = data_3773
  config.events          = 100_000
  config.decay_card      = decay_card_DpDm
  config.cross_section   = :default
end

# ============================================================================
# e+e- -> D0 anti-D0 : double tag (one D0 with charm +1, the other with charm -1)
# ============================================================================
alg_D0D0bar = TagAnalysis.new("D0D0barXS")
alg_D0D0bar.set_header(["D0D0barXSAlg/D0D0barXS.h"])
           .set_constant({"ECMS" => [:double, 3.773]})
           .set_alias({"std::vector<double>" => "Vdouble"})
           .with_decay_card(decay_card_D0D0bar)
           # tag-side charged tracks: quality cuts (|cos(theta)|<0.93, Vr<10 mm,
           # Vz<100 mm) are applied inside the DTagAlg pre-selection, not by this spec
           .note(:tag_track_selection, "tag-side charged tracks come from the pre-stored "
                 "DTagAlg collection; |cos(theta)|<0.93, Vr<10 mm, Vz<100 mm applied inside "
                 "DTagAlg/DTagTool")
           # tag-side PID: combined TOF + dE/dx, highest-probability hypothesis;
           # K_S0 daughters exempt from IP and PID requirements
           .note(:pid_correction_method, "tag-side PID combines TOF and dE/dx information and "
                 "takes the highest-probability particle hypothesis; K_S0 daughters are exempt "
                 "from the impact-parameter and PID requirements — all inside DTagAlg")
           # cosmic / QED vetoes for the D0 -> K pi tag mode
           .note(:background_veto, "the D0 -> K pi tag mode is subject to cosmic-ray and QED "
                 "vetoes built from TOF, EMC and muon-counter information")
           # photon selection used by the tag-side pi0 reconstruction
           .note(:tag_photon_selection, "photons must have E>25 MeV in the barrel "
                 "(|cos(theta)|<0.8) or E>50 MeV in the endcap (0.84<|cos(theta)|<0.92), "
                 "EMC TDC in 0-700 ns, and no associated charged track; applied inside DTagAlg")
           # pi0 -> gamma gamma tag-side reconstruction
           .note(:tag_pi0_reconstruction, "pi0 -> gamma gamma candidates require 115-150 MeV "
                 "invariant mass, at least one barrel photon, and are refit with a pi0 mass "
                 "constraint (1-C Kalman fit) inside DTagAlg")
           # K_S0 -> pi+ pi- tag-side reconstruction
           .note(:tag_ks_reconstruction, "K_S0 -> pi+ pi- candidates require a secondary-vertex "
                 "fit chi2 < 100 and 487-511 MeV invariant mass inside DTagAlg")
           # mode-dependent tag windows are stored, not cut, at BOSS level
           .note(:tag_deltaE_mbc_windows, "tag DeltaE windows are mode-dependent and asymmetric "
                 "(±3 sigma about the mean, extended to -4 sigma on the low side for pi0 modes); "
                 "D0 MBC signal region 1.858-1.874 GeV/c^2; mBC and DeltaE of the tag are stored "
                 "unconditionally and windowed in the ROOT analysis")
           # double-tag charge / charm correlation and track sharing
           .note(:tag_charge_correlation, "double-tag events require opposite net charge, "
                 "opposite charm and no tracks shared between the two tags "
                 "(enforced by DTagTool::findDTag)")

# tag side 1 : D0 with charm -1
alg_D0D0bar.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
end

# tag side 2 : D0 with charm +1  (two tag sides => double tag)
alg_D0D0bar.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm 1
end

# signal side: whatever the tag did not use — no photons
alg_D0D0bar.signal_side do |s|
  s.photons 0
end

# kinematic fit: 4-momentum conservation + both tag masses constrained to m(D0)
alg_D0D0bar.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D0)
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

alg_D0D0bar.apply
alg_D0D0bar.execute_on([data_3773, incMC_3773, exMC_D0D0bar])

# ============================================================================
# e+e- -> D+ D- : double tag (one D+ with charm +1, the other with charm -1)
# ============================================================================
alg_DpDm = TagAnalysis.new("DpDmXS")
alg_DpDm.set_header(["DpDmXSAlg/DpDmXS.h"])
        .set_constant({"ECMS" => [:double, 3.773]})
        .set_alias({"std::vector<double>" => "Vdouble"})
        .with_decay_card(decay_card_DpDm)
        .note(:tag_track_selection, "tag-side charged tracks come from the pre-stored DTagAlg "
              "collection; |cos(theta)|<0.93, Vr<10 mm, Vz<100 mm applied inside DTagAlg/DTagTool")
        .note(:pid_correction_method, "tag-side PID combines TOF and dE/dx information and takes "
              "the highest-probability particle hypothesis; K_S0 daughters are exempt from the "
              "impact-parameter and PID requirements — all inside DTagAlg")
        .note(:background_veto, "the D0 -> K pi tag mode is subject to cosmic-ray and QED vetoes "
              "built from TOF, EMC and muon-counter information")
        .note(:tag_photon_selection, "photons must have E>25 MeV in the barrel (|cos(theta)|<0.8) "
              "or E>50 MeV in the endcap (0.84<|cos(theta)|<0.92), EMC TDC in 0-700 ns, and no "
              "associated charged track; applied inside DTagAlg")
        .note(:tag_pi0_reconstruction, "pi0 -> gamma gamma candidates require 115-150 MeV "
              "invariant mass, at least one barrel photon, and are refit with a pi0 mass "
              "constraint (1-C Kalman fit) inside DTagAlg")
        .note(:tag_ks_reconstruction, "K_S0 -> pi+ pi- candidates require a secondary-vertex fit "
              "chi2 < 100 and 487-511 MeV invariant mass inside DTagAlg")
        .note(:tag_deltaE_mbc_windows, "tag DeltaE windows are mode-dependent and asymmetric "
              "(±3 sigma about the mean, extended to -4 sigma on the low side for pi0 modes); "
              "D+ MBC signal region 1.8628-1.8788 GeV/c^2; mBC and DeltaE of the tag are stored "
              "unconditionally and windowed in the ROOT analysis")
        .note(:tag_charge_correlation, "double-tag events require opposite net charge, opposite "
              "charm and no tracks shared between the two tags (enforced by DTagTool::findDTag)")

# tag side 1 : D+ with charm -1
alg_DpDm.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1
end

# tag side 2 : D+ with charm +1  (two tag sides => double tag)
alg_DpDm.tag_side(:Dplus) do |t|
  t.modes :DptoKPiPi, :DptoKPiPiPi0, :DptoKsPi, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm 1
end

# signal side: whatever the tag did not use — no photons
alg_DpDm.signal_side do |s|
  s.photons 0
end

# kinematic fit: 4-momentum conservation + both tag masses constrained to m(D+)
alg_DpDm.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:Dplus)
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:Dplus)
  f.chi2_cut 200
end

alg_DpDm.apply
alg_DpDm.execute_on([data_3773, incMC_3773, exMC_DpDm])