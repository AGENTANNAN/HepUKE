# ============================================================================
# Search for a massless dark photon γ' in Λc+ → p γ' (double-tag method)
#   ST side : ten hadronic Λc− tag modes (pKπ, Λπ, Λπππ, pKπππ families)
#   DT side : one tight proton recoiling against the ST Λc− + invisible γ'
# Tag-based analysis → TagAnalysis surface (no Selection object).
# ============================================================================

### Dataset preparation ###
# Seven c.m. energy points of the 4.5 fb^-1 data set (BOSS 703 @ 4.600 GeV, BOSS 706 for the rest)
data_4600 = DatasetManager.real_data.find("703_4600")   # 4.600 GeV
data_4612 = DatasetManager.real_data.find("706_4610")   # 4.612 GeV
data_4628 = DatasetManager.real_data.find("706_4620")   # 4.628 GeV
data_4641 = DatasetManager.real_data.find("706_4640")   # 4.641 GeV
data_4661 = DatasetManager.real_data.find("706_4660")   # 4.661 GeV
data_4682 = DatasetManager.real_data.find("706_4680")   # 4.682 GeV
data_4699 = DatasetManager.real_data.find("706_4700")   # 4.699 GeV
data_points = [data_4600, data_4612, data_4628, data_4641, data_4661, data_4682, data_4699]

# Corresponding inclusive MC samples at the same seven energy points
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")
incMC_4612 = DatasetManager.inclusive_mc.find("706_4610")
incMC_4628 = DatasetManager.inclusive_mc.find("706_4620")
incMC_4641 = DatasetManager.inclusive_mc.find("706_4640")
incMC_4661 = DatasetManager.inclusive_mc.find("706_4660")
incMC_4682 = DatasetManager.inclusive_mc.find("706_4680")
incMC_4699 = DatasetManager.inclusive_mc.find("706_4700")
incMC_points = [incMC_4600, incMC_4612, incMC_4628, incMC_4641, incMC_4661, incMC_4682, incMC_4699]

# Decay card for the signal process Λc+ → p γ' (EvtGen syntax; psi(4260) used as
# the KKMC top mother since the signal is produced directly in e+e- annihilation).
# The tag side Λc− decays are left to the generator default (their decay modes are
# fixed by the DTagAlg tag-mode selection at reconstruction time).
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Lambda_c+ anti-Lambda_c- PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000 p+ gamma_prime PHSP;
    Enddecay

    End
DECAYCARD

# 500k exclusive MC events for Λc+ → p γ' (one sample per energy point)
exMC_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_lambdac_to_p_gammaprime"
  config.events        = 500_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) — tag-based ###
alg_name = "LcPGammaPrimeDT"
lcp_alg = TagAnalysis.new(alg_name)
lcp_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
       .set_constant({ "ECMS" => [:double, 4.600] })  # beam energy is taken per run from the measured CMS
       .with_decay_card(decay_card_signal)

# ---------------------------------------------------------------------------
# ST tag side: Λc− reconstructed in the ten hadronic tag modes
# ---------------------------------------------------------------------------
lcp_alg.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKPiP,             # pKπ
          :LambdacPtoKPiPPi0,          # pKππ0 family
          :LambdacPtoKPiPPiPi,         # pKπππ family
          :LambdacPtoKPiPPiPiPi0,      # pKππππ0 family
          :LambdacPtoLambdaPi,         # Λπ
          :LambdacPtoLambdaPiPi0,      # Λππ0 family
          :LambdacPtoLambdaPiPiPi,     # Λπππ
          :LambdacPtoLambdaPiPiPiPi0,  # Λππππ0 family
          :LambdacPtoKPiPPi0Pi0,       # pKππ0π0 family
          :LambdacPtoLambdaPiPiPiPi    # Λππππ family
  t.charm(-1)                                  # the tagged side is the Λc− (c̄)
  t.window :mBC, min: 2.275, max: 2.310       # ST M_BC ∈ (2.275, 2.310) GeV/c²
end

# ---------------------------------------------------------------------------
# DT signal side: exactly one tight track (PID: proton) + invisible massless γ'
# ---------------------------------------------------------------------------
lcp_alg.signal_side do |s|
  s.charged prp: 1                      # exactly one signal-side charged track identified as a proton
  s.require_charge 1                    # proton charge (+1) recoiling against the ST Λc−
  s.missing :gamma_prime, mass: nil     # massless dark photon: p4_miss form with a per-combination residual seed
end

# ---------------------------------------------------------------------------
# 4C kinematic fit: tag + proton + missing γ' constrained to the measured CMS
# ---------------------------------------------------------------------------
lcp_alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

# ---------------------------------------------------------------------------
# BOSS-side procedures that cannot be expressed in formal DSL syntax
# ---------------------------------------------------------------------------
lcp_alg
  .note(:background_veto, "signal-side extra showers are vetoed by E_max < 0.3 GeV
    and E_sum < 0.5 GeV to suppress the π0 → γγ background; this veto also removes
    the peaking Λc+ → p K_L0 background, in which the K_L0 fakes the missing γ'
    through its shower activity")
  .note(:mode_dependent_delta_e, "the ST ΔE requirement is asymmetric and
    mode-dependent; ΔE is stored per tag mode by the tag scan and the
    mode-dependent window is applied in the ROOT analysis")
  .note(:loose_track_veto, "events containing any loose track in addition to the
    single signal-side tight track are rejected; the tight/loose track quality
    flagging is done by the DTagTool track selection")

# Render the tag spec (apply takes no Selection argument)
lcp_alg.apply

# Execute on real data, inclusive MC and the signal exclusive MC samples
root_files = lcp_alg.execute_on(data_points + incMC_points + exMC_signal)