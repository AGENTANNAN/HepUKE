# Tag-based analysis at sqrt(s) = 3.773 GeV on the psi(3770):
#   tag (anti-D0) -> K+ pi-, K+ pi- pi0, K+ pi- pi+ pi-   (hadronic single tag)
#   signal D0 -> K_S0 pi- e+ nu_e, K_S0 -> pi+ pi-

### Dataset description ###
data_3773  = DatasetManager.real_data.find("712_3773")      # psi(3770) real data
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")   # corresponding inclusive MC

# Decay card for the signal process (EvtGen format).
# The neutral kaon produced in the D0 semileptonic decay (anti-K0) is observed as K_S0 -> pi+ pi-.
decay_card_signal = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D0 anti-D0 PHSP;
  Enddecay

  Decay D0
  1.0000 K_S0 pi- e+ nu_e PHSP;
  Enddecay

  Decay anti-D0
  1.0000 K+ pi- PHSP;
  Enddecay

  Decay K_S0
  1.0000 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# Exclusive MC for the signal process (100k events)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psipp_D0_KSpi_enue"
  config.related_dataset = data_3773                # associated real dataset
  config.events          = 100000                   # 100k events
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Tag-based event selection (BOSS) ###
alg_name = "D0TagKSpiEnu"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.773]})       # sqrt(s) = 3.773 GeV
   .set_alias({"std::vector<double>" => "Vdouble"})
   .with_decay_card(decay_card_signal)

# Tag side: the anti-D0 is reconstructed as a single tag in three hadronic modes
alg.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1                                         # pin the tagged side to the anti-D0
  t.window :deltaE, abs: 0.02                        # opt-in deltaE window (mode-dependent windows refined in ROOT)
end

# Signal side: D0 -> K_S0 pi- e+ nu_e with K_S0 -> pi+ pi-
alg.signal_side do |s|
  s.charged(pip: 1, pim: 2, ep: 1)                   # pi+ pi- (from K_S0) + pi- + e+
  s.require_charge 0                                 # net charge 0
  s.missing :nu_e                                    # massless missing neutrino
end

# Kinematic fit: 4C constraint plus D mass constraint, loose chi^2 cut
alg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:pip, :pim).between(0.4856, 0.5096)                            # |M(pipi) - m_K_S0| < 12 MeV
  f.invariant_mass_of(:pip, :pim, :pim, :ep, :nu_e).constrain_to_nominal_mass_of(:D0) # D mass constraint
  f.chi2_cut 200
end

# Inexpressible BOSS-side procedures
alg.note(:secondary_vertex_fit, "K_S0 reconstructed from pi+ pi- via a secondary-vertex fit; the decay-length "
                                "significance of the K_S0 vertex is required to exceed 2 sigma (no tag-DSL method)")
alg.note(:pid_selection_cuts,   "signal-side PID: pions identified with CL_pi > CL_K; the electron candidate "
                                "requires E/p > 0.8 and an EMC shower shape consistent with an electron "
                                "(signal-side PID thresholds are fixed v1 defaults, not DSL-tunable)")
alg.note(:background_veto,      "events consistent with D0 -> K_S0 pi+ pi- are vetoed against the "
                                "D0 -> K_S0 pi- e+ nu_e signal")
alg.note(:tag_deltaE_windows,   "the deltaE windows are tag-mode-dependent; a single representative window is "
                                "declared in the tag selection and the per-mode windows are refined in ROOT "
                                "(deltaE is stored unconditionally)")

alg.apply

# Execute on real data, inclusive MC, and the signal exclusive MC
root_files = alg.execute_on([data_3773, incMC_3773, exMC_signal])