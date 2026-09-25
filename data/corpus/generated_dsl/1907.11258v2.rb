# ============================================================
# ψ(3770) at 3.773 GeV : D -> φ + pseudoscalar (φ -> K+K-)
#   Mode I   : D+ -> φ π+                 (tag D- -> anti-K0 π-)
#   Mode II  : D+ -> φ K+                 (tag D- -> anti-K0 π-)
#   Mode III : D0 -> φ π0 (π0 -> γγ)      (tag anti-D0 -> K+ π-)
#   Mode IV  : D0 -> φ η  (η  -> γγ)      (tag anti-D0 -> K+ π-)
# ============================================================

### Dataset description ###
data_3773  = DatasetManager.real_data.find("712_3773")      # ψ(3770) real data at 3.773 GeV
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")   # corresponding inclusive MC

# ---- Decay cards for the four signal processes (EvtGen syntax) ----
# Mode I : ψ(3770) -> D+ D-, D+ -> φ π+, tag D- -> anti-K0 π-
decay_card_I = <<~DECAYCARD
  Decay psi(3770)
  1.0 D+ D- PHSP;
  Enddecay

  Decay D+
  1.0 phi pi+ PHSP;
  Enddecay

  Decay D-
  1.0 anti-K0 pi- PHSP;
  Enddecay

  Decay phi
  1.0 K+ K- VSS;
  Enddecay

  Decay anti-K0
  1.0 K_S0 PHSP;
  Enddecay

  Decay K_S0
  1.0 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# Mode II : ψ(3770) -> D+ D-, D+ -> φ K+, tag D- -> anti-K0 π-
decay_card_II = <<~DECAYCARD
  Decay psi(3770)
  1.0 D+ D- PHSP;
  Enddecay

  Decay D+
  1.0 phi K+ PHSP;
  Enddecay

  Decay D-
  1.0 anti-K0 pi- PHSP;
  Enddecay

  Decay phi
  1.0 K+ K- VSS;
  Enddecay

  Decay anti-K0
  1.0 K_S0 PHSP;
  Enddecay

  Decay K_S0
  1.0 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# Mode III : ψ(3770) -> D0 anti-D0, D0 -> φ π0, tag anti-D0 -> K+ π-
decay_card_III = <<~DECAYCARD
  Decay psi(3770)
  1.0 D0 anti-D0 PHSP;
  Enddecay

  Decay D0
  1.0 phi pi0 PHSP;
  Enddecay

  Decay anti-D0
  1.0 K+ pi- PHSP;
  Enddecay

  Decay phi
  1.0 K+ K- VSS;
  Enddecay

  Decay pi0
  1.0 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Mode IV : ψ(3770) -> D0 anti-D0, D0 -> φ η, tag anti-D0 -> K+ π-
decay_card_IV = <<~DECAYCARD
  Decay psi(3770)
  1.0 D0 anti-D0 PHSP;
  Enddecay

  Decay D0
  1.0 phi eta PHSP;
  Enddecay

  Decay anti-D0
  1.0 K+ pi- PHSP;
  Enddecay

  Decay phi
  1.0 K+ K- VSS;
  Enddecay

  Decay eta
  1.0 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# ---- Exclusive MC : 100k events for each of the four modes ----
exMC_I = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_DpToPhiPi"
  config.related_dataset = data_3773
  config.events          = 100000
  config.decay_card      = decay_card_I
  config.cross_section   = :default
end

exMC_II = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_DpToPhiK"
  config.related_dataset = data_3773
  config.events          = 100000
  config.decay_card      = decay_card_II
  config.cross_section   = :default
end

exMC_III = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_D0ToPhiPi0"
  config.related_dataset = data_3773
  config.events          = 100000
  config.decay_card      = decay_card_III
  config.cross_section   = :default
end

exMC_IV = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_D0ToPhiEta"
  config.related_dataset = data_3773
  config.events          = 100000
  config.decay_card      = decay_card_IV
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# Four independent decay modes -> four Algorithm objects (Rule T1).

alg_I   = Algorithm.new("DpToPhiPi")
alg_II  = Algorithm.new("DpToPhiK")
alg_III = Algorithm.new("D0ToPhiPi0")
alg_IV  = Algorithm.new("D0ToPhiEta")

alg_I.set_header(["DpToPhiPiAlg/DpToPhiPi.h"])
     .set_constant({"ECMS" => [:double, 3.773]})
     .set_alias({"std::vector<double>" => "Vdouble"})

alg_II.set_header(["DpToPhiKAlg/DpToPhiK.h"])
      .set_constant({"ECMS" => [:double, 3.773]})
      .set_alias({"std::vector<double>" => "Vdouble"})

alg_III.set_header(["D0ToPhiPi0Alg/D0ToPhiPi0.h"])
       .set_constant({"ECMS" => [:double, 3.773]})
       .set_alias({"std::vector<double>" => "Vdouble"})

alg_IV.set_header(["D0ToPhiEtaAlg/D0ToPhiEta.h"])
      .set_constant({"ECMS" => [:double, 3.773]})
      .set_alias({"std::vector<double>" => "Vdouble"})

# ---- Common charged-track selection ----
# Charged-D+ modes (I, II) : >=2 positive, >=1 negative, net charge +1
charged_common = Selection.new
    .select_track {
      cos_theta 0.93      # |cosθ| < 0.93
      Vz        100.0     # |Vz| < 100 cm
      Vr        10.0      # Vr < 10 cm
      nChrp     ">=2"     # at least two positive tracks
      nChrn     ">=1"     # at least one negative track
      nNet      "==1"     # net charge +1
    }

# Neutral-D0 modes (III, IV) : >=1 positive, >=1 negative, net charge 0
# plus photon selection (25 MeV barrel / 50 MeV endcap, TDC 0-14, >=2 photons,
# isolated by more than 10° from any charged track)
neutral_common = Selection.new
    .select_track {
      cos_theta 0.93
      Vz        100.0
      Vr        10.0
      nChrp     ">=1"     # at least one positive track
      nChrn     ">=1"     # at least one negative track
      nNet      "==0"     # net charge 0
    }
    .select_photon {
      tdc_emc_start     0
      tdc_emc_end       14        # EMC TDC window (0-700 ns)
      energyThreshold_b 0.025     # 25 MeV in the barrel
      energyThreshold_e 0.050     # 50 MeV in the endcap
      angle_to_track    10.0      # isolated by > 10° from any charged track
      nGam              ">=2"     # at least two photons
    }

# ---- Mode I : D+ -> φ π+ (K+K-π+) ----
sel_I = charged_common.dup
    .pid(method: :probability) {
      prob_cut 0.001                                  # probability > 0.001
      identify :kaon, against: [:pion, :proton]       # K+ and K-
      identify :pion, against: [:kaon, :proton]       # π+
      nkp  "==1"
      nkm  "==1"
      npip "==1"
    }
    .kinematic_fit([:kp, :km, :pip]) {
      nominal                  # nominal fit: corrected 4-momenta are saved
      constrain_four_momentum  # 4C energy-momentum constraint
      chi2_cut 200             # loose χ² < 200 (tight cut applied in ROOT)
    }

# ---- Mode II : D+ -> φ K+ (K+K+K-) ----
sel_II = charged_common.dup
    .pid(method: :probability) {
      prob_cut 0.001
      identify :kaon, against: [:pion, :proton]       # K-vs-π/proton separation
      nkp "==2"
      nkm "==1"
    }
    .kinematic_fit([:kp, :kp, :km]) {
      nominal
      constrain_four_momentum
      chi2_cut 200
    }

# ---- Mode III : D0 -> φ π0 (K+K-γγ, π0 reconstructed by Kalman fit) ----
sel_III = neutral_common.dup
    .pid(method: :probability) {
      prob_cut 0.001
      identify :kaon, against: [:pion, :proton]
      nkp "==1"
      nkm "==1"
    }
    .kalman_kinematic_fit([:gamma, :gamma]) {
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # γγ mass -> M(π0)
      chi2_cut 25
      npi0 ">=1"
    }
    .kinematic_fit([:kp, :km, :pi0]) {
      nominal
      constrain_four_momentum
      chi2_cut 200
    }

# ---- Mode IV : D0 -> φ η (K+K-γγ, η reconstructed by Kalman fit) ----
sel_IV = neutral_common.dup
    .pid(method: :probability) {
      prob_cut 0.001
      identify :kaon, against: [:pion, :proton]
      nkp "==1"
      nkm "==1"
    }
    .kalman_kinematic_fit([:gamma, :gamma]) {
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)  # γγ mass -> M(η)
      chi2_cut 25
      neta ">=1"
    }
    .kinematic_fit([:kp, :km, :eta]) {
      nominal
      constrain_four_momentum
      chi2_cut 200
    }

# ---- Generate the complete algorithm for each mode and run it ----
alg_I.note(:tag_side_reconstruction, "events are required to be tagged by the opposite-side D- -> anti-K0 pi-; the tag-side reconstruction is not part of this signal-side selection")
     .with_decay_card(decay_card_I).apply(sel_I)
root_files_I = alg_I.execute_on([data_3773, incMC_3773, exMC_I])

alg_II.note(:tag_side_reconstruction, "events are required to be tagged by the opposite-side D- -> anti-K0 pi-; the tag-side reconstruction is not part of this signal-side selection")
      .with_decay_card(decay_card_II).apply(sel_II)
root_files_II = alg_II.execute_on([data_3773, incMC_3773, exMC_II])

alg_III.note(:tag_side_reconstruction, "events are required to be tagged by the opposite-side anti-D0 -> K+ pi-; the tag-side reconstruction is not part of this signal-side selection")
       .with_decay_card(decay_card_III).apply(sel_III)
root_files_III = alg_III.execute_on([data_3773, incMC_3773, exMC_III])

alg_IV.note(:tag_side_reconstruction, "events are required to be tagged by the opposite-side anti-D0 -> K+ pi-; the tag-side reconstruction is not part of this signal-side selection")
      .with_decay_card(decay_card_IV).apply(sel_IV)
root_files_IV = alg_IV.execute_on([data_3773, incMC_3773, exMC_IV])