### Dataset preparation ###
# 4.600 GeV real data (567 pb^-1) and the corresponding inclusive MC
data_4600  = DatasetManager.real_data.find("703_4600")       # select the 4.600 GeV real dataset
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")    # corresponding inclusive MC sample

# Decay card for e+e- -> Lambda_c+ anti-Lambda_c- produced in phase space
# (top mother psi(4260) is the BESIII KKMC convention when no intermediate resonance is given)
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Lambda_c+ anti-Lambda_c- PHSP;
    Enddecay
    End
DECAYCARD

# Exclusive MC: 500k events, phase-space Lambda_c+ anti-Lambda_c- (models both tag modes)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4600_LambdacLambdacbar"
  config.related_dataset = data_4600           # associate with the 4.600 GeV real dataset
  config.events          = 500000              # 500k events
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) — tag-based (double-tag) analysis ###
alg_name = "LambdacTagSemilep"
alg = TagAnalysis.new(alg_name)                # TagAnalysis (subclass of Algorithm); no Selection object
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 4.600]}) # declare the CMS energy (used by the 4C fit)
   .with_decay_card(decay_card_signal)

# Tag side: hadronic anti-Lambda_c- tags (charm -1 pins the anti-Lambda_c- side)
alg.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKsP, :LambdacPtoKPiP      # pbar K_S0 and pbar K+ pi-
  t.charm -1
end

# Signal side: inclusive semileptonic Lambda_c+ -> X e+ nu_e
# (everything the tag did NOT use; no explicit |cos(theta)|/|Vz|/|Vr| cuts are set)
alg.signal_side do |s|
  s.charged(ep: 1, at_least: true)             # at least one positron; other charged tracks from X allowed
  s.missing :nu_e                              # massless missing neutrino
end

# 4C kinematic fit: tag system + e+ + nu constrained to the lab-frame CMS four-momentum.
# The missing-neutrino declaration automatically stores the fitted missing mass
# (m_P4_miss_fit), the unfitted missing-mass observables (m_Umiss, m_Umiss2) and q^2 (m_q2).
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200                               # chi^2 < 200
end

# Positron PID criteria are stricter than the fixed SimplePIDSvc electron recipe used by
# the tag analysis, so they are preserved as a note.
alg.note(:pid_correction_method,
         "signal-side positron PID requires CL(e) > 0.001, " \
         "CL(e)/(CL(e)+CL(pi)+CL(K)+CL(p)) > 0.8, and E_e/p_e > 0.8; the tag-analysis DSL " \
         "applies only the fixed v1 SimplePIDSvc electron recipe, so these criteria and the " \
         "hadron-misidentification PID efficiency unfolding matrix are applied in the ROOT analysis")

alg.apply                                      # takes no Selection argument (tag-based surface)

# Execute the algorithm on data, inclusive MC and signal exclusive MC
root_files = alg.execute_on([data_4600, incMC_4600, exMC_signal])