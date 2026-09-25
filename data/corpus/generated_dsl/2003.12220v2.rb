### Dataset preparation ###
data_3773  = DatasetManager.real_data.find("712_3773")     # ψ(3770) real data (2.93 fb⁻¹)
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")  # corresponding inclusive MC

# Exclusive MC: ψ(3770) → D⁺D⁻ phase space (1M events)
decay_card_dd = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D+ D- PHSP;
    Enddecay
    End
DECAYCARD

exMC_DD = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_DpDm_phsp"
  config.related_dataset = data_3773
  config.events          = 1_000_000
  config.decay_card      = decay_card_dd
  config.cross_section   = :default
end

### Tag-based event selection (BOSS) ###
alg_name = "DpToEtaMuNu"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({ "ECMS" => [:double, 3.773] })   # √s = 3.773 GeV
   .with_decay_card(decay_card_dd)
   # Signal-muon PID: fixed v1 DTag lepton-PID recipe (CL_μ > 0.001, CL_μ > CL_e, CL_μ > CL_K)
   # plus a calorimeter requirement 0.105 < E_EMC(μ) < 0.275 GeV — no dedicated tag-DSL form.
   .note(:pid_correction_method,
         "signal μ⁺ identified with the fixed v1 DTag lepton-PID recipe (CL_mu > 0.001, " \
         "CL_mu > CL_e, CL_mu > CL_K); an additional muon EMC-energy requirement " \
         "0.105 < E_EMC(mu) < 0.275 GeV is applied on the signal muon")
   # Extra-shower background suppression on the signal side.
   .note(:background_veto,
         "extra-shower suppression on the signal side: maximum extra photon energy < 0.30 GeV " \
         "and zero extra π⁰ (no additional γγ pair consistent with the π⁰ mass)")
   # Mode-dependent ΔΕ windows: DSL emits only one window per observable, so the inclusive
   # union is declared here; per-mode tight windows are applied later in ROOT.
   .note(:tag_deltaE_mode_dependent_window,
         "mode-dependent tag ΔΕ windows (non-π⁰ modes: -0.025 < ΔΕ < 0.025 GeV; " \
         "π⁰ modes: -0.055 < ΔΕ < 0.045 GeV) are applied in ROOT; BOSS emits the inclusive " \
         "union window (-0.055, 0.045) GeV")

# Tag the D⁻ in six hadronic modes; best candidate per mode/charge (smallest |ΔΕ|) is
# picked internally by DTagAlg. Tag windows are opt-in (store-not-cut default).
alg.tag_side(:Dm) do |t|
  t.modes :DptoKPiPi, :DptoKsPi, :DptoKPiPiPi0, :DptoKsPiPi0, :DptoKsPiPiPi, :DptoKKPi
  t.charm -1                                  # pin the tagged side to D⁻
  t.window :deltaE, min: -0.055, max: 0.045   # inclusive union of the mode-dependent ΔΕ windows
  t.window :mBC,    min: 1.863,  max: 1.877   # M_BC ∈ (1.863, 1.877) GeV/c²
end

# Signal side: D⁺ → η μ⁺ ν_μ, η → γγ (tag's unused showers/tracks + missing ν_μ)
alg.signal_side do |s|
  s.photons 2       # two photons for η → γγ
  s.charged(mup: 1) # one μ⁺ (lepton key → SimplePID/ParticleID muon recipe, fixed v1 thresholds)
  s.missing :nu_mu  # massless missing ν_μ (semileptonic tag; auto-stores U_miss, q²)
end

# 4C kinematic fit with the two photons constrained to the nominal η mass;
# M(ημ⁺) < 1.74 GeV/c² veto to reject D⁺ → η π⁺ (pion–muon misidentification).
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.invariant_mass_of(:gamma, :gamma, :mup).between(0.0, 1.74)
  f.chi2_cut 200
end

alg.apply   # TagAnalysis#apply takes NO Selection argument
alg.execute_on([data_3773, incMC_3773, exMC_DD])