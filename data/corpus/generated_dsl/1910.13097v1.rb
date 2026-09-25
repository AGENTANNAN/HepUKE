# ============================================================================
# Dataset preparation — √s = 4.178 GeV (3.19 fb⁻¹)
# ============================================================================
data_4180  = DatasetManager.real_data.find("703_4180")      # 3.19 fb⁻¹ real data at 4.178 GeV
incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")   # matching inclusive MC sample

# Decay card for the signal process e+e- → D_s+ D_s*- (EvtGen syntax)
#   tag side   : D_s*- → γ D_s-   (D_s- reconstructed in three hadronic modes)
#   signal side: D_s+ → p p̄ e+ ν_e
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 D_s+ D_s*- PHSP;
    Enddecay

    Decay D_s*-
    1.0000 gamma D_s- PHSP;
    Enddecay

    Decay D_s+
    1.0000 p+ anti-p- e+ anti-nu_e PHSP;
    Enddecay

    Decay D_s-
    0.3333 K+ K- pi-                  PHSP;
    0.3333 K_S0 K-                    PHSP;
    0.3334 K_S0 K- pi+ pi-            PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive MC for the signal mode
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4180_dsp_dsst_ppepnu"
  config.related_dataset = data_4180
  config.events          = 100000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

# ============================================================================
# Tag analysis (D_s tag): ST D_s*- → γ D_s-  vs  D_s+ → p p̄ e+ ν_e
# ============================================================================
alg_name = "DsTagPPbarEnu"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 4.178]})
   # --- Inexpressible BOSS-side tag-side procedures, preserved for systematics ---
   .note(:tag_pion_momentum, "tag-side pions (D_s- → K+K-π-, K_S K-, K_S K-π+π-) required to " \
                             "have momentum > 0.1 GeV/c; no DSL track-momentum primitive exists " \
                             "for tag-internal tracks")
   .note(:tag_ks_selection, "K_S candidates on the tag side required to have M(π+π-) in " \
                            "(0.487, 0.511) GeV/c² and a decay length greater than twice its " \
                            "uncertainty; these are secondary-vertex quantities of the tag, not " \
                            "covered by the tag-side mBC/ΔE window vocabulary")
   .note(:tag_recoil_mass_window, "D_s- recoil mass required in [2.06, 2.18] GeV/c²; the " \
                                  "observable is M(D_s- γ) built from the tag candidate and the " \
                                  "radiative photon (windowed in the ROOT analysis)")
   .note(:tag_best_candidate, "when several tag candidates survive, the one whose M(D_s- γ) is " \
                              "closest to the D_s*+ nominal mass (2.1122 GeV/c²) is kept; " \
                              "candidate ranking is performed in the ROOT analysis")
   .note(:unused_track_multiplicity, "fewer than four unused charged tracks are allowed opposite " \
                                     "the tag; this is an event-level tag-side track-count " \
                                     "requirement with no dedicated DSL primitive")
   .note(:electron_reconstruction_efficiency, "D_s+ → p p̄ e+ ν_e: the e+ is not reconstructed in " \
                                              "~95% of signal events, so events are classified as " \
                                              "signal by the p p̄ pair alone; the residual electron " \
                                              "efficiency is handled by re-weighting the signal MC")
   .note(:soft_third_track_momentum, "when a third signal-side track is found, it is required to " \
                                     "have momentum < 0.09 GeV/c; this soft-track cut is applied " \
                                     "in the ROOT analysis")
   .note(:missing_mass_squared, "missing-mass squared > 0 required to suppress continuum q q̄ " \
                                "background; the missing-mass observables are stored " \
                                "unconditionally (m_P4_miss_fit, m_Umiss, m_Umiss2) and the cut " \
                                "is applied in the ROOT analysis")
   .with_decay_card(decay_card_signal)

# Tag side — D_s*- → γ D_s-, tagged D_s- (charge pinned to the negative side)
alg.tag_side(:Ds) do |t|
  t.modes :DstoKKPi, :DstoKsK, :DstoKsKPiPi   # D_s- → K+K-π-, K_S K-, K_S K-π+π-
  t.charm(-1)
end

# Signal side — everything the tag did not use:
#   the radiative photon from D_s*- → γ D_s-, and D_s+ → p p̄ e+ ν_e
alg.signal_side do |s|
  s.photons 1                                   # γ from D_s*- → γ D_s-
  s.charged(prp: 1, prm: 1, ep: 1)              # one p, one p̄, one e+
  s.require_charge 1                            # net charge of the signal side = +1
  s.missing :nu_e                               # massless missing ν_e
end

# Kinematic fit — 4C: tag + signal + neutrino four-momenta constrained to the CMS energy
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

# No Selection object — the tag blocks are the selection; apply takes no argument
alg.apply

# Execute on real data, inclusive MC and the signal exclusive MC
root_files = alg.execute_on([data_4180, incMC_4180, exMC_signal])