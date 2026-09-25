# =============================================================================
# Semileptonic Λc+ → p K− e+ νe  —  tag-based (ST Λc− tag + DT signal side)
# =============================================================================

### Dataset description ###
# 7 c.m. energy points spanning 4.600 – 4.699 GeV (total 4.5 fb−1)
data_points = [
  DatasetManager.real_data.find("703_4600"),   # 4599.53 MeV
  DatasetManager.real_data.find("706_4610"),   # 4611.86 MeV
  DatasetManager.real_data.find("706_4620"),   # 4628.00 MeV
  DatasetManager.real_data.find("706_4640"),   # 4640.91 MeV
  DatasetManager.real_data.find("706_4660"),   # 4661.24 MeV
  DatasetManager.real_data.find("706_4680"),   # 4681.92 MeV
  DatasetManager.real_data.find("706_4700")    # 4698.82 MeV
]
incMC_points = [
  DatasetManager.inclusive_mc.find("703_4600"),
  DatasetManager.inclusive_mc.find("706_4610"),
  DatasetManager.inclusive_mc.find("706_4620"),
  DatasetManager.inclusive_mc.find("706_4640"),
  DatasetManager.inclusive_mc.find("706_4660"),
  DatasetManager.inclusive_mc.find("706_4680"),
  DatasetManager.inclusive_mc.find("706_4700")
]

# Decay card for the signal process (EvtGen format; ψ(4260) as KKMC top mother)
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Lambda_c+ anti-Lambda_c-         PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000 p+ K- e+ nu_e                    PHSP;
    Enddecay

    End
DECAYCARD

# 500k signal MC events, the same Λc+ → pK−e+νe signal generated at each energy point
exMC_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_LcpKepnu"
  config.events        = 500_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Tag analysis (BOSS) ###
alg_name = "LcpKepnuTag"
tag_analysis = TagAnalysis.new(alg_name)
tag_analysis.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.600]})  # per-run measured beam energy is used by the fit at every scan point
            .note(:dt_signal_reconstruction,
                  "the 'double tag' of this analysis is expressed as a single hadronic " \
                  "Λc− tag plus a signal side of exactly three recoiling charged tracks " \
                  "(p K− e+) and a missing νe; DTagTool has no Λc double-tag support " \
                  "(findDTag rejects decayMode() >= 1000), so the DT information is " \
                  "obtained from the ST tag plus the recoiling signal tracks")
            .note(:mode_dependent_delta_e,
                  "the ST Λc− selection also requires a mode-dependent ΔE window; only " \
                  "one window per observable is expressible on the tag side, so the " \
                  "explicit M_BC window is applied in BOSS and the 14 per-mode ΔE " \
                  "windows are applied in the ROOT analysis")
            .note(:background_veto,
                  "π0 veto applied in the DT selection to suppress backgrounds with a " \
                  "π0 in the signal side (e.g. Λc+ → pK−π+ with a misidentified π+)")
            .note(:pid_correction_method,
                  "signal-side e+ PID requires CL_e > 0.001, CL_e > CL_π and CL_e > CL_K; " \
                  "the signal-side electron key uses the fixed PID recipe of the tag layer, " \
                  "the residual difference is treated as a systematic uncertainty")

# ---------------------------------------------------------------------------
# Tag side: single tag Λc− reconstructed in 14 hadronic tag modes
# ---------------------------------------------------------------------------
tag_analysis.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKPiP,          # Λc− → K+ π− p̄  (pK−π+)
          :LambdacPtoKPiPPi0,       # pK−π+π0
          :LambdacPtoKPiPPiPi,      # pK−π+π+π−
          :LambdacPtoKPiPPiPiPi0,   # pK−π+π+π−π0
          :LambdacPtoKsP,           # pK_S
          :LambdacPtoKsPPi0,        # pK_Sπ0
          :LambdacPtoKsPPiPi,       # pK_Sπ+π−
          :LambdacPtoLambdaPi,      # Λπ+
          :LambdacPtoLambdaPiPi0,   # Λπ+π0
          :LambdacPtoLambdaPiPiPi,  # Λπ+π+π−
          :LambdacPtoSigma0Pi,      # Σ0π+
          :LambdacPtoSigmaPPi0,     # Σ+π0
          :LambdacPtoSigmaPPiPi     # Σ+π+π−
  t.charm -1                        # the tagged side is Λc−
  # explicit ST M_BC window asked for in the description
  t.window :mBC, min: 2.28, max: 2.30
end

# ---------------------------------------------------------------------------
# Signal side: the tracks the tag did not use  →  p K− e+ with a missing νe
# ---------------------------------------------------------------------------
tag_analysis.signal_side do |s|
  s.charged(kp: 1, km: 1, ep: 1)    # exactly 3 charged tracks recoiling against the ST Λc−
  s.require_charge(1)               # Λc+ → p K− e+ νe
  s.missing :nu_e                   # massless missing neutrino (semileptonic tag mode)
end

# ---------------------------------------------------------------------------
# Kinematic fit: 4-momentum conservation with a missing neutrino
# (U_miss = E_miss − |p_miss|, m_Umiss, m_q2 are stored automatically)
# ---------------------------------------------------------------------------
tag_analysis.fit do |f|
  f.constrain_four_momentum                              # 4C fit with the missing νe
  f.invariant_mass_of(:kp, :km, :ep).between(0.0, 2.15)  # M(pK−e+) < 2.15 GeV
  f.chi2_cut 200
end

# Generate the algorithm for the signal decay card and run on all datasets
tag_analysis.with_decay_card(decay_card_signal).apply
root_files = tag_analysis.execute_on(data_points + incMC_points + exMC_signal)