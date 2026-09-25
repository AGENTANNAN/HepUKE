# ============================================================================
# ψ(3770) → D+D− ,  D+ → μ+ν_μ   (double-tag method)
#   tag side   : single-tag D− in eight hadronic modes
#   signal side: D+ → μ+ν_μ reconstructed as one μ+ plus a missing ν_μ
# ============================================================================

### Dataset description ###
psipp_data  = DatasetManager.real_data.find("712_3773")     # ψ(3770) real data (~20.3 fb⁻¹)
psipp_incMC = DatasetManager.inclusive_mc.find("712_3773")  # matching inclusive MC

# Decay card for the exclusive signal process (EvtGen format)
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- PHSP;
    Enddecay

    Decay D+
    1.0000 mu+ nu_mu PHSP;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi- PHSP;
    Enddecay

    End
DECAYCARD

# 1,000,000-event exclusive signal MC
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_DpToMuNu"
  config.related_dataset = psipp_data
  config.events          = 1_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Tag-based event selection (TagAnalysis) ###
alg = TagAnalysis.new("DpToMuNu")
alg.set_header(["DpToMuNuAlg/DpToMuNu.h"])
   .set_constant({"ECMS" => [:double, 3.773]})            # ψ(3770) c.m. energy
   .with_decay_card(decay_card_signal)

# ---- tag side: single-tag D− in eight hadronic modes ----
alg.tag_side(:Dm) do |t|
  t.modes :DptoKPiPi,     # D− → K+π−π−
          :DptoKsPi,      # D− → K_Sπ−
          :DptoKPiPiPi0,  # D− → K+π−π−π0
          :DptoKsPiPi0,   # D− → K_Sπ−π0
          :DptoKsPiPiPi,  # D− → K_Sπ+π−π−
          :DptoKKPi,      # D− → K+K−π−
          :DptoPiPiPi,    # D− → π+π−π−
          :DptoKPiPiPiPi  # D− → K+π−π+π−π−
  t.charm -1              # pin the tagged particle to D−
end

# ---- signal side: exactly one μ+ (Q = +1), no photons, missing ν_μ ----
alg.signal_side do |s|
  s.charged(mup: 1)     # exactly one μ+ ; default exact total → no extra charged track
  s.require_charge(1)   # total signal-side charge = +1
  s.photons 0           # no photon participates in the fit
  s.missing :nu_mu      # undetected (massless) neutrino
end

# ---- kinematic fit: tag + μ+ + missing ν constrained to the lab-frame 4-momentum ----
alg.fit do |f|
  f.constrain_four_momentum   # tag + μ+ + ν_μ = measured CMS four-momentum (ecms_lab)
  f.chi2_cut 200              # χ² < 200
end

# ---- BOSS-side procedures with no dedicated DSL method ----
alg.note(:muon_identification,
         "signal-side μ+ identified by MUC hit depth versus momentum and cos(theta), " \
         "rather than the default tag-layer ParticleID probability recipe")
   .note(:signal_muon_emc_energy,
         "EMC energy of the μ+ candidate required to lie within [0.00, 0.35] GeV")
   .note(:extra_shower_veto,
         "maximum energy of any shower not used in the reconstruction required to be below 0.3 GeV")

# ---- render and run ----
alg.apply                                     # takes no Selection argument
root_files = alg.execute_on([psipp_data, psipp_incMC, exMC_signal])