# =============================================================================
# e+e- -> D*0 D*- pi+  —  partial reconstruction with two tagging methods
#   Tag 1 (D0-tag) : reconstruct D*0 -> D0 pi0 (+ bachelor pi+), leave D*- missing
#   Tag 2 (D--tag) : reconstruct D*- -> D- pi0 (+ bachelor pi+), leave D*0 missing
# =============================================================================

### ---------------------- Dataset preparation ---------------------- ###
# XYZ (and scan) data points between 4.189 and 4.951 GeV
data_names = %w[
  703_4190 703_4200 703_4210 703_4220 703_4230 703_4237 703_4245 703_4246
  703_4260 703_4270 703_4280 703_4310 703_4360 703_4390 703_4420 703_4470
  703_4530 703_4575 703_4600
  703_4190scan 703_4210scan 703_4220scan 703_4230scan
  705_4290 705_4315 705_4340 705_4380 705_4400 705_4440
  706_4610 706_4620 706_4640 706_4660 706_4680 706_4700
  707_4740 707_4750 707_4780 707_4840 707_4914 707_4946
]
real_data = data_names.map { |n| DatasetManager.real_data.find(n) }

# Matching inclusive MC samples at the same energy points
inc_mc_names = %w[
  703_4190 703_4200 703_4210 703_4220 703_4230 703_4237 703_4246 703_4260
  703_4270 703_4280 703_4360 703_4420 703_4600
  703_4190scan 703_4210scan 703_4220scan 703_4230scan
  706_4610 706_4620 706_4640 706_4660 706_4680 706_4700
  707_4740 707_4750 707_4780 707_4840 707_4914 707_4946
]
inc_mc = inc_mc_names.map { |n| DatasetManager.inclusive_mc.find(n) }

# Decay card for the signal process e+e- -> D*0 D*- pi+ (EvtGen syntax)
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.000 D*0 D*- pi+ PHSP;
    Enddecay
    Decay D*0
    1.000 D0 pi0 PHSP;
    Enddecay
    Decay D*-
    1.000 D- pi0 PHSP;
    Enddecay
    Decay D0
    1.000 K- pi+ PHSP;
    Enddecay
    Decay D-
    1.000 K+ pi- pi- PHSP;
    Enddecay
    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

# 200k-event exclusive signal MC per energy point (energy scan -> create_exclusive_mc_for)
exMC_signals = DatasetManager.create_exclusive_mc_for(real_data) do |config|
  config.sample_name   = "exmc_Dstar0Dstarmpi"
  config.events        = 200_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### ---------------------- Event selection (BOSS) ---------------------- ###
# ---------------- Tag 1 : D0-tag (D*0 reconstructed, D*- left missing) ----------------
alg_d0 = Algorithm.new("D0TagPartialRec")
alg_d0.set_header(["D0TagPartialRecAlg/D0TagPartialRec.h"])
      .set_constant({"ECMS" => [:double, 4.26]})          # representative CMS energy of the scan
      .set_alias({"std::vector<double>" => "Vdouble"})

sel_d0 = Selection.new
  .select_track {                                          # charged track quality / multiplicity
      cos_theta 0.93
      Vz        10.0
      Vr        1.0
      nChrp     ">=2"                                      # at least two positive tracks
      nChrn     ">=2"                                      # at least two negative tracks
      nTot      ">=5"                                      # at least five tracks in total
  }
  .select_photon {                                         # good photons
      tdc_emc_start     0
      tdc_emc_end       14
      angle_to_track    10.0
      energyThreshold_b 0.025                              # > 25 MeV in the barrel
      energyThreshold_e 0.050                              # > 50 MeV in the endcap
      nGam              ">=2"                              # at least two photons
  }
  .pid(method: :probability) {                             # probability PID, K/pi separation
      prob_cut 0.001
      identify :kaon, against: [:pion, :proton]            # K+ and K-
      identify :pion, against: [:kaon, :proton]            # pi+ and pi-
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {                # 1C Kalman fit -> pi0
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 25
      npi0     ">=1"
  }
  .partial_miss([2]) {                                     # miss recID 2 = D*- ; reconstruct D*0 + pi+
      best_combination_by_mass :D0,     1.861              # D0 -> K-pi+ window 1.835-1.887 GeV
      best_combination_by_mass :Dstar0, 2.0065             # M(D0 pi0) window 2.004-2.009 GeV
      require_recoil_mass 1.80, 2.20                       # missing-D* recoil mass window
  }

alg_d0
  .note(:background_veto, "P*(pi0) in the D*0 rest frame required outside 0.025-0.050 GeV (veto against mis-tagged D*0).")
  .note(:supplementary_cut, "M(pi+ D0) > 2.02 GeV required to reject combinatorial background.")
  .note(:tag_arbitration, "When both the D0-tag and the D--tag selection survive, the D0-tag event is retained.")
  .note(:mass_constrained_fit, "3C mass-constrained fit over pi0, D and D* nominal masses (chi2 < 50) performed on the reconstructed tag; the partial-reconstruction step already keeps the mass-closest combination.")
  .note(:additional_decay_modes, "D0 -> K-pi+pi0 (window 1.827-1.882 GeV) and D0 -> K-pi+pi+pi- (window 1.855-1.874 GeV) reconstructed in dedicated Algorithms with their own decay cards.")
  .with_decay_card(decay_card_signal)
  .apply(sel_d0)

# ---------------- Tag 2 : D--tag (D*- reconstructed, D*0 left missing) ----------------
alg_dm = Algorithm.new("DmTagPartialRec")
alg_dm.set_header(["DmTagPartialRecAlg/DmTagPartialRec.h"])
      .set_constant({"ECMS" => [:double, 4.26]})
      .set_alias({"std::vector<double>" => "Vdouble"})

sel_dm = Selection.new
  .select_track {
      cos_theta 0.93
      Vz        10.0
      Vr        1.0
      nChrp     ">=2"
      nChrn     ">=2"
      nTot      ">=5"
  }
  .select_photon {
      tdc_emc_start     0
      tdc_emc_end       14
      angle_to_track    10.0
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      nGam              ">=2"
  }
  .pid(method: :probability) {
      prob_cut 0.001
      identify :kaon, against: [:pion, :proton]
      identify :pion, against: [:kaon, :proton]
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 25
      npi0     ">=1"
  }
  .partial_miss([1]) {                                     # miss recID 1 = D*0 ; reconstruct D*- + pi+
      best_combination_by_mass :Dm,     1.8695             # D- -> K+pi-pi- window 1.856-1.883 GeV
      best_combination_by_mass :Dstarm, 2.0105             # M(D- pi0) window 2.008-2.013 GeV
      require_recoil_mass 1.80, 2.20
  }

alg_dm
  .note(:background_veto, "P*(pi0) in the D*- rest frame required outside 0.030-0.055 GeV (veto against mis-tagged D*-).")
  .note(:supplementary_cut, "M(pi+ D0) > 2.02 GeV required to reject combinatorial background.")
  .note(:mass_constrained_fit, "3C mass-constrained fit over pi0, D and D* nominal masses (chi2 < 50) performed on the reconstructed tag.")
  .with_decay_card(decay_card_signal)
  .apply(sel_dm)

# ---------------- Execute on real data + inclusive MC + signal MC ----------------
root_files_d0 = alg_d0.execute_on(real_data + inc_mc + exMC_signals)
root_files_dm = alg_dm.execute_on(real_data + inc_mc + exMC_signals)