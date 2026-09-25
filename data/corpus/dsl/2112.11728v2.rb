# 2112.11728v2: Measurement of R = σ(e+e-→hadrons) / σ(e+e-→μ+μ-)
# at 14 CMS energies from 2.2324 to 3.6710 GeV.
# Inclusive hadronic event selection — counting experiment.
# The BOSS DSL is designed for exclusive channel reconstruction; almost all
# steps of this analysis are captured via algorithm.note().
#
# Energy points: 2.2324, 2.4000, 2.8000, 3.0500, 3.0600, 3.0800, 3.4000,
#   3.5000, 3.5424, 3.5538, 3.5611, 3.6002, 3.6500, 3.6710 GeV.
# Represented below by a single example energy; the full analysis requires
# separate jobOptions per energy point.
# Dataset sample names follow BESIII convention: [BOSS]_[Energy_in_MeV].

# Placeholder decay card (inclusive hadronic events; used only to satisfy the
# Algorithm#apply requirement that a decay card is attached)
decay_card_placeholder = <<~DECAYCARD
    Decay vpho
    1.000 hadrons PHSP;
    Enddecay
    End
DECAYCARD

# ── Algorithm ────────────────────────────────────────────────────
alg = Algorithm.new("RValueHadronic")
alg.set_header(["RValueHadronicAlg/RValueHadronic.h"])
   .set_constant({ "ECMS" => [:double, 3.400] })
   .with_decay_card(decay_card_placeholder)

# Minimal event selection — only track selection expressed in DSL;
# the bulk of the selection logic lives in .note() entries below.
event_selection = Selection.new
event_selection.select_track {
  cos_theta 0.93
  Vz 5.0
  Vr 0.5
  nChrp ">=2"
  nChrn ">=2"
  nTot  ">=2"
}
.select_photon {
  tdc_emc_start 0
  tdc_emc_end 700
  energyThreshold_b 0.100
  energyThreshold_e 0.100
  nGam ">=1"
}
.pid(method: :probability) {
  prob_cut 0.001
  identify :pion, against: [:kaon, :proton]
}
# Minimal kinematic fit to satisfy the DSL requirement; the actual analysis
# does NOT use a kinematic fit — event counting is purely topology-based.
.kinematic_fit([:pip, :pim]) {
  nominal
  constrain_four_momentum
  chi2_cut 200
}

alg
  .note(:inclusive_analysis, "This is an inclusive hadronic counting experiment. No exclusive final state is reconstructed. The BOSS DSL does not natively model inclusive hadronic event selection, so almost all selection logic is captured via these notes.")
  .note(:multi_energy, "Analysis performed at 14 CMS energies: 2.2324, 2.4000, 2.8000, 3.0500, 3.0600, 3.0800, 3.4000, 3.5000, 3.5424, 3.5538, 3.5611, 3.6002, 3.6500, 3.6710 GeV. Separate jobOptions per energy point. Dataset sample names follow [BOSS]_[Energy_in_MeV] convention; look up exact names in BES3_dataset.md.")
  .note(:qed_veto, "Reject e+e-→e+e- and e+e-→γγ: require at least 2 EMC showers with |Δθ| = |θ1+θ2−180°| < 10° and 2nd-most-energetic shower E > 0.65·E_beam.")
  .note(:deuteron_removal, "Remove non-collision charged tracks (deuterons) via χ_p < 10 where χ_p = (dE/dx − dE/dx_p)/σ_p.")
  .note(:momentum_cuts, "Charged track p < 0.94·p_beam; remove if E/(pc) > 0.8 AND p > 0.65·p_beam; remove gamma-conversion pairs: oppositely charged, both E/(pc) > 0.8, M_inv < 0.1 GeV/c², opening angle < 15°.")
  .note(:two_prong_selection, "2-prong events: tracks must NOT be back-to-back (|Δθ| < 10° AND ||Δφ|−180°| < 15°); N_iso_2prg > 1 (isolated photon: E > 100 MeV, angle to charged track > 20°).")
  .note(:three_prong_selection, "3-prong events: two highest-p tracks NOT back-to-back; n_tracks with E/(pc) > 0.8 must be < 2; n_tracks with r_PID > 0.25 must be < 2.")
  .note(:four_plus_prong, "Events with > 3 prongs directly counted as hadronic without additional requirements.")
  .note(:beam_background, "Beam-associated background estimated via V_z_evt sideband (5,10) cm vs signal (0,5) cm, cross-checked with separated-beam data and double-Gaussian fit method.")
  .note(:r_value_formula, "R = (N_obs_had − N_bkg) / [σ⁰(μμ)·L_int·ε_trig·ε_had·(1+δ)]. L_int from large-angle Bhabha. ε_trig ≈ 100%. (1+δ) ISR factor via FD scheme with vacuum polarization from dispersion relations. ε_had from luarlw generator (tuned to data). Iterative procedure for self-consistent (1+δ) and σ⁰_had.")
  .note(:background_sources, "QED backgrounds (e+e-, γγ, μ+μ-, τ+τ-, two-photon processes) estimated via babayaga, kkmc, diag36, ekhara, galuga generators. Residual contamination subtracted from N_obs_had.")
  .note(:systematic_uncertainties, "Dominant systematic uncertainties from inclusive hadronic simulation model (luarlw vs hybrid conexc+phokhara+luarlw) and ISR correction scheme (FD vs structure function), each at the ~1-2% level. Total uncertainty < 3.0%.")
  .note(:no_kinematic_fit_actual, "The kinematic_fit block above is a DSL-requirement placeholder. The actual analysis performs no kinematic fit — event counting is purely based on track/shower topology criteria. The fit participants [:pip, :pim] are not physically meaningful for this analysis.")
  .apply(event_selection)

# Example execution for a single energy point; repeat for all 14 energies
# alg.execute_on([DatasetManager.real_data.find("703_3400"),
#                 DatasetManager.inclusive_mc.find("703_3400")])
# Note: exclusive MC not used in this analysis (inclusive hadronic generator instead).