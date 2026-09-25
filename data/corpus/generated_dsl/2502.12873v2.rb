# =====================================================================
# ψ(3770) → D0 D̄0 strong-phase difference measurement between
# D0 → K+K−π+π− and D̄0 → K+K−π+π− (quantum-correlated DD̄, C = −1)
# Double-tag technique: tag side from DTagAlg candidates,
# signal side = fully reconstructed D0 → K+K−π+π−.
# =====================================================================

### Dataset preparation ###
dd_data  = DatasetManager.real_data.find("712_3773")      # ψ(3770) real data @ 3.773 GeV
dd_incMC = DatasetManager.inclusive_mc.find("712_3773")   # corresponding inclusive MC

# Signal decay card: ψ(3770) → D0 D̄0, D0 → K+K−π+π−, D̄0 → K+π− (phase space)
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.000  D0  anti-D0   PHSP;
    Enddecay

    Decay D0
    1.000  K+  K-  pi+  pi-   PHSP;
    Enddecay

    Decay anti-D0
    1.000  K+  pi-   PHSP;
    Enddecay

    End
DECAYCARD

# 2M-event exclusive signal MC
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_d0d0bar_kkpipi_kpi"
  config.related_dataset = dd_data
  config.events          = 2_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Tag analysis (BOSS) ###
alg_name = "D0StrongPhase"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.773]})

# ---- Tag side: the tagging D0 from DTagAlg pre-stored candidates ----
#      flavor / CP-even / CP-odd / mixed-CP / semileptonic tags,
#      including the partially reconstructed modes with undetected
#      K_L0 or ν_e (handled inside the DTag reconstruction itself).
alg.tag_side(:D0) do |t|
  t.modes :D0toKPi,          # K- pi+                          (flavor)
          :D0toKPiPi0,       # K- pi+ pi0                      (flavor)
          :D0toKPiPiPi,      # K- pi+ pi- pi+                  (flavor)
          :D0toKK,           # K+ K-                           (CP-even)
          :D0toPiPi,         # pi+ pi-                         (CP-even)
          :D0toKsPi0Pi0,     # K_S0 pi0 pi0                    (CP-even)
          :D0toPiPiPi0,      # pi+ pi- pi0                     (CP-even)
          :D0toKL0Pi0,       # K_L0 pi0                        (CP-even, partially reco.)
          :D0toKsPi0,        # K_S0 pi0                        (CP-odd)
          :D0toKsEta,        # K_S0 eta                        (CP-odd)
          :D0toKsEtaPrime,   # K_S0 eta' (eta'->pi+pi-eta)     (CP-odd)
          :D0toKsPiPiPi0,    # K_S0 pi+ pi- pi0                (CP-odd)
          :D0toKsPiPi,       # K_S0 pi+ pi-                    (mixed-CP)
          :D0toKL0PiPi,      # K_L0 pi+ pi-                    (mixed-CP, partially reco.)
          :D0toKeNu          # K- e+ nu_e                      (semileptonic, partially reco.)
end

# ---- Signal side: the other D0 → K+K−π+π− (fully reconstructed) ----
alg.signal_side do |s|
  s.charged(kp: 1, km: 1, pip: 1, pim: 1)   # one K+, one K-, one pi+, one pi-
  s.require_charge(0)                       # net charge zero
end

# ---- 4C kinematic fit (no missing particle), χ² < 200 ----
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

# BOSS-side procedures that have no dedicated DSL construct
alg.note(:background_veto,
         "any pi+pi- pair in the signal-side D0 -> K+K-pi+pi- with invariant mass in
          [0.477, 0.507] GeV/c^2 is vetoed to suppress K_S0 -> pi+pi- contamination;
          the window is +-3 sigma around the K_S0 peak")
      .note(:best_candidate_selection,
            "for the double-tag selection both D candidates must be reconstructed and the
             combination whose M_BC is closest to the D0 nominal mass is retained")

alg.with_decay_card(decay_card_signal).apply   # no Selection argument for a TagAnalysis
alg.execute_on([dd_data, dd_incMC, exMC_signal])