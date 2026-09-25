# ============================================================================
# BESIII  ψ(3770) @ √s = 3.773 GeV  (712_3773 real data + inclusive MC)
# D+ → D0 e+νe   —   double-tag (semileptonic-tag) analysis, BOSS part.
#
#   Tag side    : D− reconstructed from pre-stored tag candidates (DTagTool)
#                 in the six standard hadronic modes.
#   Signal side : the D0 recoiling against the tag; the soft positron
#                 (p < 5 MeV/c) is not reconstructed and is merged with the
#                 neutrino into ONE massless missing four-vector.
# ============================================================================

### Datasets ###
data_3773  = DatasetManager.real_data.find("712_3773")     # ψ(3770) real data
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")  # ψ(3770) inclusive MC

### Decay cards (EvtGen format) - one per signal D0 mode ###
# Signal mode I :  D0 → K−π+
decay_card_D0_KPi = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D-    PHSP;
    Enddecay

    Decay D+
    1.0000 D0 e+ nu_e    PHSP;
    Enddecay

    Decay D0
    1.0000 K- pi+    PHSP;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi-    PHSP;
    Enddecay

    End
DECAYCARD

# Signal mode II :  D0 → K−π+π+π−
decay_card_D0_K3Pi = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D-    PHSP;
    Enddecay

    Decay D+
    1.0000 D0 e+ nu_e    PHSP;
    Enddecay

    Decay D0
    1.0000 K- pi+ pi+ pi-    PHSP;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi-    PHSP;
    Enddecay

    End
DECAYCARD

# Signal mode III :  D0 → K−π+π0
decay_card_D0_KPiPi0 = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D-    PHSP;
    Enddecay

    Decay D+
    1.0000 D0 e+ nu_e    PHSP;
    Enddecay

    Decay D0
    1.0000 K- pi+ pi0    PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma    PHSP;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi-    PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive MC : 200k events for each D0 signal mode ###
exMC_D0_KPi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0toKPi_eNu"
  config.related_dataset = data_3773
  config.events          = 200000
  config.decay_card      = decay_card_D0_KPi
  config.cross_section   = :default
end

exMC_D0_K3Pi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0toKPiPiPi_eNu"
  config.related_dataset = data_3773
  config.events          = 200000
  config.decay_card      = decay_card_D0_K3Pi
  config.cross_section   = :default
end

exMC_D0_KPiPi0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0toKPiPi0_eNu"
  config.related_dataset = data_3773
  config.events          = 200000
  config.decay_card      = decay_card_D0_KPiPi0
  config.cross_section   = :default
end

# The six standard D− tag modes.
# Channel names are given in the D+ convention; `charm -1` pins the D− tag.
tag_modes_Dminus = [
  :DptoKPiPi,      # D− → K+ π− π−
  :DptoKPiPiPi0,   # D− → K+ π− π− π0
  :DptoKsPi,       # D− → K_S0 π−
  :DptoKsPiPi0,    # D− → K_S0 π− π0
  :DptoKsPiPiPi,   # D− → K_S0 π+ π− π−
  :DptoKKPi,       # D− → K+ K− π−
]

# ----------------------------------------------------------------------------
# Signal mode I :  D+ → D0(→ K−π+) e+νe
# ----------------------------------------------------------------------------
alg_KPi = TagAnalysis.new("DplusToD0eNu_KPi")
alg_KPi.set_header(["DplusToD0eNu_KPiAlg/DplusToD0eNu_KPi.h"])
       .set_constant({ "ECMS" => [:double, 3.773] })
       .with_decay_card(decay_card_D0_KPi)
       .note(:deltaE_window, "The D− tag ΔE window is mode-dependent (~3σ) and " \
                            "DTagTool keeps, per tag mode, the candidate with the " \
                            "smallest |ΔE|. Only a single representative symmetric " \
                            "window can be declared at BOSS level; the final " \
                            "per-mode ΔE windows are applied in the downstream ROOT stage.")

alg_KPi.tag_side(:Dm) do |t|          # one tag side → single D− tag per event
  t.modes(*tag_modes_Dminus)          # the six standard D− tag modes
  t.charm -1                          # pin the tagged D to be D−
  t.window :mBC,    min: 1.83         # MBC > 1.83 GeV/c²
  t.window :deltaE, abs: 0.03         # representative ~3σ |ΔE| window
end

alg_KPi.signal_side do |s|
  s.charged(km: 1, pip: 1)            # D0 → K− π+
  s.require_charge 0                  # signal-side net charge zero
  s.min_photon_angle 10.0             # 10° minimum photon angle
  s.missing :nu_e                     # soft e+ merged with ν → one massless missing 4-vector
end

alg_KPi.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:km, :pip).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

alg_KPi.apply
root_files_KPi = alg_KPi.execute_on([data_3773, incMC_3773, exMC_D0_KPi])

# ----------------------------------------------------------------------------
# Signal mode II :  D+ → D0(→ K−π+π+π−) e+νe
# ----------------------------------------------------------------------------
alg_K3Pi = TagAnalysis.new("DplusToD0eNu_K3Pi")
alg_K3Pi.set_header(["DplusToD0eNu_K3PiAlg/DplusToD0eNu_K3Pi.h"])
        .set_constant({ "ECMS" => [:double, 3.773] })
        .with_decay_card(decay_card_D0_K3Pi)
        .note(:deltaE_window, "The D− tag ΔE window is mode-dependent (~3σ) and " \
                             "DTagTool keeps, per tag mode, the candidate with the " \
                             "smallest |ΔE|. Only a single representative symmetric " \
                             "window can be declared at BOSS level; the final " \
                             "per-mode ΔE windows are applied in the downstream ROOT stage.")

alg_K3Pi.tag_side(:Dm) do |t|
  t.modes(*tag_modes_Dminus)
  t.charm -1
  t.window :mBC,    min: 1.83
  t.window :deltaE, abs: 0.03
end

alg_K3Pi.signal_side do |s|
  s.charged(km: 1, pip: 2, pim: 1)    # D0 → K− π+ π+ π−
  s.require_charge 0
  s.min_photon_angle 10.0
  s.missing :nu_e
end

alg_K3Pi.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:km, :pip, :pip, :pim).constrain_to_nominal_mass_of(:D0)
  f.chi2_cut 200
end

alg_K3Pi.apply
root_files_K3Pi = alg_K3Pi.execute_on([data_3773, incMC_3773, exMC_D0_K3Pi])

# ----------------------------------------------------------------------------
# Signal mode III :  D+ → D0(→ K−π+π0) e+νe
# ----------------------------------------------------------------------------
alg_KPiPi0 = TagAnalysis.new("DplusToD0eNu_KPiPi0")
alg_KPiPi0.set_header(["DplusToD0eNu_KPiPi0Alg/DplusToD0eNu_KPiPi0.h"])
          .set_constant({ "ECMS" => [:double, 3.773] })
          .with_decay_card(decay_card_D0_KPiPi0)
          .note(:deltaE_window, "The D− tag ΔE window is mode-dependent (~3σ) and " \
                                "DTagTool keeps, per tag mode, the candidate with the " \
                                "smallest |ΔE|. Only a single representative symmetric " \
                                "window can be declared at BOSS level; the final " \
                                "per-mode ΔE windows are applied in the downstream ROOT stage.")
          .note(:pi0_1c_fit, "π0 candidates (π0 → γγ) must first pass a standalone 1C " \
                             "kinematic fit constraining m(γγ) to the nominal π0 mass " \
                             "with χ² < 20. This separate mass-constrained reconstruction " \
                             "is not expressible inside the TagAnalysis fit block and is " \
                             "performed by the BOSS-side π0 routine before the tag fit.")
          .note(:endcap_photon_veto, "Events are rejected when both photons of the π0 " \
                                     "candidate are reconstructed in the EMC endcap.")

alg_KPiPi0.tag_side(:Dm) do |t|
  t.modes(*tag_modes_Dminus)
  t.charm -1
  t.window :mBC,    min: 1.83
  t.window :deltaE, abs: 0.03
end

alg_KPiPi0.signal_side do |s|
  s.photons 2..2                      # exactly two photons (π0 → γγ)
  s.min_photon_energy 0.025           # shower energy > 25 MeV
  s.min_photon_angle 10.0             # 10° minimum photon angle
  s.charged(km: 1, pip: 1)            # D0 → K− π+
  s.require_charge 0
  s.missing :nu_e
end

alg_KPiPi0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)   # π0 mass
  f.invariant_mass_of(:gamma, :gamma).between(0.110, 0.155)               # diphoton mass window
  f.invariant_mass_of(:km, :pip, :gamma, :gamma).constrain_to_nominal_mass_of(:D0)  # D0 mass
  f.chi2_cut 200
end

alg_KPiPi0.apply
root_files_KPiPi0 = alg_KPiPi0.execute_on([data_3773, incMC_3773, exMC_D0_KPiPi0])