# ============================================================
# Paper: Amplitude Analysis and Branching Fraction Measurement of
#        Ds+ → K+ K- pi+ at BESIII
# arXiv: 2011.08041v4
# ============================================================
# Tag-based (Double Tag) analysis: e+e- → Ds* Ds → γ Ds(tag) Ds(signal)
# Tag side: Ds- through 8 hadronic tag modes
# Signal side: Ds+ → K+ K- pi+  (also a tag mode: :DstoKKPi)
# Extra: γ from Ds* → Ds γ
# ============================================================

###
### Dataset: √s = 4.178 GeV, BOSS 703, 3.19 fb-1
###
data_4178 = DatasetManager.real_data.find("703_4180")
incMC_4178 = DatasetManager.inclusive_mc.find("703_4180")

# Decay card for the overall process
decay_card = <<~DECAYCARD
    Decay psi(4260)
    1.000 D_s+ D_s*-  PHSP;
    Enddecay

    Decay D_s*-
    1.000 gamma D_s-  VSP_PWAVE;
    Enddecay

    Decay D_s+
    1.000 K+ K- pi+  PHSP;
    Enddecay

    End
DECAYCARD

exMC_dt = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_DsDT_KKPi"
  config.related_dataset = data_4178
  config.events          = 500_000
  config.decay_card      = decay_card
  config.cross_section   = :default
end

###
### Double-Tag TagAnalysis:
###   Tag side 1: Ds- through 8 hadronic tag modes (charm -1)
###   Tag side 2: Ds+ through K+K-pi+ signal mode (charm +1)
###   Signal side: γ from Ds* → Ds γ
###
alg_dt = TagAnalysis.new("DsDTKKPi")
alg_dt.set_header(["DsDTKKPiAlg/DsDTKKPi.h"])
alg_dt.set_constant({"ECMS" => [:double, 4.178]})

# Tag side 1: Ds- (tag side) through 8 hadronic modes
alg_dt.tag_side(:Ds) do |t|
  t.modes :DstoKKPi, :DstoKsK, :DstoKsKplusPiPi, :DstoKPiPi,
          :DstoKsKminusPiPi, :DstoPiPiPi, :DstoPiPiPiEta, :DstoKKPiPi0
  t.charm -1
end

# Tag side 2: Ds+ (signal side) through K+K-pi+ mode
alg_dt.tag_side(:Ds) do |t|
  t.modes :DstoKKPi
  t.charm 1
  t.rank_by :inv
end

# Signal side: the γ from Ds* → Ds γ
alg_dt.signal_side do |s|
  s.photons 1
  s.min_photon_angle 10.0
end

# 6C kinematic fit: 4C + both Ds mass constraints
alg_dt.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:Ds)
  f.invariant_mass_of(:tag2).constrain_to_nominal_mass_of(:Ds)
  f.chi2_cut 200
end

alg_dt.note(:tag_mass_windows, "Tag Ds mass windows per mode (see Table IV): [1.940-1.996] GeV/c2 depending on mode; signal Ds mass window: [1.950, 1.986] GeV/c2")
alg_dt.note(:recoil_mass_window, "Recoil mass M_rec in [2.051, 2.180] GeV/c2 to retain both direct Ds production and Ds*→Dsγ")
alg_dt.note(:amplitude_analysis, "7C kinematic fit (5C + Ds mass constraints on both sides) performed in ROOT-level amplitude fit; BOSS level provides DT sample via TagAnalysis")
alg_dt.note(:mipwa_method, "MIPWA analysis path uses separate ST reconstruction with 1C kinematic fit and BDTG classifier; not expressible via TagAnalysis — requires separate Algorithm with partial reconstruction")
alg_dt.note(:background_veto, "Soft pi± and pi0 with momentum < 0.1 GeV/c vetoed to suppress D* decay products")
alg_dt.note(:photon_conversion_veto, "Ks mass window 0.487-0.511 GeV/c2; secondary vertex L/sigma_L > 2; pi0 mass window 0.115-0.150 GeV/c2 with chi2<30; eta mass window 0.490-0.580 GeV/c2 with chi2<30")

alg_dt.with_decay_card(decay_card)
        .apply

alg_dt.execute_on([data_4178, incMC_4178, exMC_dt])