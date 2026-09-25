### Dataset preparation ###
# Ten center-of-mass energy points, 4.600 - 4.840 GeV
data_sample_names = %w[703_4600 706_4610 706_4620 706_4640 706_4660
                       706_4680 706_4700 707_4740 707_4750 707_4840]
real_datasets   = data_sample_names.map { |n| DatasetManager.real_data.find(n) }     # real data at each point
inc_mc_datasets = data_sample_names.map { |n| DatasetManager.inclusive_mc.find(n) } # matching inclusive MC

# Decay card for e+e- -> Lambda_c+ Lambda_c-
#   tag side   : anti-Lambda_c- -> p- K+ pi-
#   signal side: Lambda_c+ -> p+ pi0 (pi0 -> gamma gamma)
decay_card_pi0 = <<~DECAYCARD
    Decay psi(4260)
    1.0 Lambda_c+ anti-Lambda_c- PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0 anti-p- K+ pi- PHSP;
    Enddecay

    Decay Lambda_c+
    1.0 p+ pi0 PHSP;
    Enddecay

    Decay pi0
    1.0 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Decay card with the eta signal channel: Lambda_c+ -> p+ eta (eta -> gamma gamma)
decay_card_eta = <<~DECAYCARD
    Decay psi(4260)
    1.0 Lambda_c+ anti-Lambda_c- PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0 anti-p- K+ pi- PHSP;
    Enddecay

    Decay Lambda_c+
    1.0 p+ eta PHSP;
    Enddecay

    Decay eta
    1.0 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC: same card for each energy point of the scan, one set per signal channel
exMC_pi0 = DatasetManager.create_exclusive_mc_for(real_datasets) do |config|
  config.sample_name   = "sig_LcLc_pLc_pi0"
  config.events        = 100000
  config.decay_card    = decay_card_pi0
  config.cross_section = :default
end

exMC_eta = DatasetManager.create_exclusive_mc_for(real_datasets) do |config|
  config.sample_name   = "sig_LcLc_pLc_eta"
  config.events        = 100000
  config.decay_card    = decay_card_eta
  config.cross_section = :default
end

### Event selection (BOSS): tag-based double-tag analysis ###

# ---------------- Signal channel I : Lambda_c+ -> p+ pi0 ----------------
alg_pi0 = TagAnalysis.new("LcDtagPi0")
alg_pi0.set_header(["LcDtagPi0Alg/LcDtagPi0.h"])
       .set_constant({"ECMS" => [:double, 4.600]})   # nominal; fit uses the per-run measured beam energy

# Tag side: anti-Lambda_c- reconstructed through the nine hadronic modes
alg_pi0.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKPiP, :LambdacPtoKsP, :LambdacPtoKPiPi0P, :LambdacPtoKsPi0P,
          :LambdacPtoKsPiPiP, :LambdacPtoLambdaPi, :LambdacPtoLambdaPiPi0,
          :LambdacPtoLambdaPiPiPi, :LambdacPtoSigmaPi
  t.charm -1                                        # pin the tagged side to anti-Lambda_c-
end

# Signal side: exactly one proton and two photons (photon pair -> pi0 candidate)
alg_pi0.signal_side do |s|
  s.charged(prp: 1)                                 # exactly one signal-side proton
  s.photons 2                                       # two signal-side photons
  s.min_photon_angle 30.0                           # photon-proton opening angle > 30 degrees
  s.min_photon_energy 0.025                         # photon energy > 25 MeV
end

# 4C kinematic fit (chi2 < 200); gamma-gamma mass window selects pi0
alg_pi0.fit do |f|
  f.constrain_four_momentum                         # four-momentum conservation
  f.chi2_cut 200                                    # chi2 < 200
  f.invariant_mass_of(:gamma, :gamma).between(0.115, 0.150)   # pi0 mass window
end

# BOSS-side procedures with no dedicated DSL construct
alg_pi0.note(:delta_e_cut,
             "signal-side DeltaE_{p2gamma} required in (-0.080, 0.035) GeV after the 4C fit")
       .note(:background_veto,
             "veto M(p_signal, pi-_tag) in [1.111, 1.121] GeV/c^2 (Lambda) and " \
             "M(pi0_signal, pi+pi-_tag) in [0.733, 0.833] GeV/c^2 (omega)")
       .note(:photon_shower_shape,
             "signal-side photons additionally require lateral moment in [0.05, 0.40] and " \
             "E3x3/E5x5 > 0.85; endcap energy threshold 50 MeV (barrel 25 MeV)")
       .note(:candidate_selection,
             "when several candidates remain, keep the signal candidate with the smallest " \
             "|DeltaE_{p2gamma}| and the tag candidate with the smallest |DeltaE|")

alg_pi0.apply
alg_pi0.execute_on(real_datasets + inc_mc_datasets + exMC_pi0)

# ---------------- Signal channel II : Lambda_c+ -> p+ eta ----------------
alg_eta = TagAnalysis.new("LcDtagEta")
alg_eta.set_header(["LcDtagEtaAlg/LcDtagEta.h"])
       .set_constant({"ECMS" => [:double, 4.600]})   # nominal; fit uses the per-run measured beam energy

# Tag side: identical nine hadronic anti-Lambda_c- modes
alg_eta.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKPiP, :LambdacPtoKsP, :LambdacPtoKPiPi0P, :LambdacPtoKsPi0P,
          :LambdacPtoKsPiPiP, :LambdacPtoLambdaPi, :LambdacPtoLambdaPiPi0,
          :LambdacPtoLambdaPiPiPi, :LambdacPtoSigmaPi
  t.charm -1
end

# Signal side: exactly one proton and two photons (photon pair -> eta candidate)
alg_eta.signal_side do |s|
  s.charged(prp: 1)
  s.photons 2
  s.min_photon_angle 30.0
  s.min_photon_energy 0.025
end

# 4C kinematic fit (chi2 < 200); gamma-gamma mass window selects eta
alg_eta.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
  f.invariant_mass_of(:gamma, :gamma).between(0.490, 0.583)   # eta mass window
end

alg_eta.note(:delta_e_cut,
             "signal-side DeltaE_{p2gamma} required in (-0.080, 0.035) GeV after the 4C fit")
       .note(:background_veto,
             "veto M(p_signal, pi-_tag) in [1.111, 1.121] GeV/c^2 (Lambda) and " \
             "M(pi0_signal, pi+pi-_tag) in [0.733, 0.833] GeV/c^2 (omega)")
       .note(:photon_shower_shape,
             "signal-side photons additionally require lateral moment in [0.05, 0.40] and " \
             "E3x3/E5x5 > 0.85; endcap energy threshold 50 MeV (barrel 25 MeV)")
       .note(:candidate_selection,
             "when several candidates remain, keep the signal candidate with the smallest " \
             "|DeltaE_{p2gamma}| and the tag candidate with the smallest |DeltaE|")

alg_eta.apply
alg_eta.execute_on(real_datasets + inc_mc_datasets + exMC_eta)