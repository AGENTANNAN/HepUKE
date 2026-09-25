# BESIII first amplitude analysis and BF measurement of Ds+ -> pi+ pi0 pi0 eta.
# Double-tag method with Ds- tag; 7.33 fb^-1 at Ds threshold (sqrt(s) ~ 4.178 GeV).

# ==== Dataset (Ds sample near 4.178 GeV) ====
ds_data  = DatasetManager.real_data.find("706_4178")
ds_incMC = DatasetManager.inclusive_mc.find("706_4178")

# ---------------- Decay card for signal Ds+ -> pi+ pi0 pi0 eta ----------------
dc_signal = <<~DC
  Decay D_s*+
  1.0000 gamma D_s+ VSP_PWAVE;
  Enddecay

  Decay D_s+
  1.0000 pi+ pi0 pi0 eta PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DC

exmc_signal = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "Dsp_pipi0pi0eta"
  c.related_dataset = ds_data
  c.events          = 500000
  c.decay_card      = dc_signal
  c.cross_section   = :default
end

# ============ Tag analysis: Ds- tag, signal Ds+ -> pi+ pi0 pi0 eta ============
DS_TAG_MODES = [
  :DstoKsK,
  :DstoKKPi,
  :DstoKsKPi0,
  :DstoKKPiPi0,
  :DstoPiPiPi,
  :DstoPiEta,
  :DstoPiPi0Eta,
  :DstoPiEtaPrime,
  :DstoKPiPi
]

alg = TagAnalysis.new("DsTagPipPi0Pi0Eta")
alg.set_header(["DsTagPipPi0Pi0EtaAlg/DsTagPipPi0Pi0Eta.h"])
   .set_constant({ "ECMS" => [:double, 4.178] })
   .with_decay_card(dc_signal)
   .note(:transition_photon,
         "Ds*+ -> gamma Ds+ transition photon selected on the signal side; " \
         "photon energy in lab frame required < 0.18 GeV based on MC")
   .note(:mrec_dsstar,
         "recoiling mass against Ds*+ (signal Ds+ + transition photon) in [1.93,1.98] GeV/c2, " \
         "consistent with Ds- mass")
   .note(:fake_pi0_veto,
         "veto photon combinations (from pi0/eta/Ds*+ transition) with invariant mass in " \
         "[0.10,0.14] GeV/c2 that could form a fake pi0")
   .note(:multi_cand_selection,
         "for multiple candidates keep the one with minimum chi2 of the 9C kinematic fit " \
         "(4-momentum conservation + pi0/eta/Ds- mass constraints + Ds*+ mass constraint on " \
         "either Ds+gamma or Ds-gamma combination)")

alg.tag_side(:Ds) do |t|
  t.modes(*DS_TAG_MODES)
  t.charm(-1)
end

# Signal side: pi+, two pi0 (each -> gamma gamma), one eta (-> gamma gamma), plus transition gamma
alg.signal_side do |s|
  s.photons 7          # 2*pi0 (4 gammas) + eta (2 gammas) + 1 transition gamma
  s.charged(pip: 1)
  s.require_charge 1
  s.min_photon_angle 10.0
end

# 9C-like kinematic fit: 4-momentum conservation plus mass constraints on pi0, eta and Ds tag.
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:D_s_minus)
  f.chi2_cut 200
end

alg.apply
alg.execute_on([ds_data, ds_incMC, exmc_signal])
