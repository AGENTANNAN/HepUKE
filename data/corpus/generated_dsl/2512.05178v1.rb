# frozen_string_literal: true

# ============================================================================
# Dataset preparation — Lambda_c+ scan at sqrt(s) = 4.600, 4.612, 4.628,
# 4.641, 4.661, 4.682, 4.699 GeV (tag analysis: e+e- -> Lambda_c+ Lambda_c-)
# ============================================================================
scan_points = {
  "703_4600" => 4.600,
  "706_4610" => 4.612,
  "706_4620" => 4.628,
  "706_4640" => 4.641,
  "706_4660" => 4.661,
  "706_4680" => 4.682,
  "706_4700" => 4.699,
}
data_points  = scan_points.keys.map { |s| DatasetManager.real_data.find(s) }     # real data at each energy point
incMC_points = scan_points.keys.map { |s| DatasetManager.inclusive_mc.find(s) }  # matching inclusive MC at each point

# ---------------------------------------------------------------------------
# Decay card — mode 1: Lambda_c+ -> Sigma+ pi- e+ nu_e, Sigma+ -> p pi0,
#                     generated tag side anti-Lambda_c- -> anti-p K_S0
# ---------------------------------------------------------------------------
decay_card_mode1 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay

  Decay Lambda_c+
  1.0000 Sigma+ pi- e+ nu_e PHSP;
  Enddecay

  Decay anti-Lambda_c-
  1.0000 anti-p- K_S0 PHSP;
  Enddecay

  Decay Sigma+
  1.0000 p+ pi0 PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# ---------------------------------------------------------------------------
# Decay card — modes 2+3: Lambda_c+ -> Sigma+/- pi∓ e+ nu_e,
#                         Sigma+/- -> n pi+/-,
#                         generated tag side anti-Lambda_c- -> anti-p K_S0
# ---------------------------------------------------------------------------
decay_card_mode23 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay

  Decay Lambda_c+
  0.5000 Sigma+ pi- e+ nu_e PHSP;
  0.5000 Sigma- pi+ e+ nu_e PHSP;
  Enddecay

  Decay anti-Lambda_c-
  1.0000 anti-p- K_S0 PHSP;
  Enddecay

  Decay Sigma+
  1.0000 n0 pi+ PHSP;
  Enddecay

  Decay Sigma-
  1.0000 n0 pi- PHSP;
  Enddecay

  End
DECAYCARD

# 200k-event exclusive signal MC at every one of the seven scan points
exMC_mode1 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_lambdac_scan_mode1_sigmaP"
  config.events        = 200_000
  config.decay_card    = decay_card_mode1
  config.cross_section = :default
end

exMC_mode23 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_lambdac_scan_mode23_sigmaN"
  config.events        = 200_000
  config.decay_card    = decay_card_mode23
  config.cross_section = :default
end

# ============================================================================
# Event selection — Mode 1: Lambda_c+ -> Sigma+ pi- e+ nu_e, Sigma+ -> p pi0
# ============================================================================
alg_name_mode1 = "LambdacSemilepSigmaP"
alg_mode1 = TagAnalysis.new(alg_name_mode1)
alg_mode1.set_header(["#{alg_name_mode1}Alg/#{alg_name_mode1}.h"])
         .set_constant({ "ECMS" => [:double, 4.600] })
         .with_decay_card(decay_card_mode1)

# Tag side: single tag anti-Lambda_c- through the 12 hadronic modes (both charges scanned)
alg_mode1.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKsP, :LambdacPtoKPiP, :LambdacPtoKsPi0P, :LambdacPtoKsPiPiP,
          :LambdacPtoKPiPi0P, :LambdacPtoLambdaPi, :LambdacPtoLambdaPiPi0,
          :LambdacPtoLambdaPiPiPi, :LambdacPtoSigma0Pi, :LambdacPtoSigmaPi0,
          :LambdacPtoSigmaPiPi, :LambdacPtoPiPiP
  # charm omitted -> both tag charges scanned
end

# Signal side: p, pi-, e+, and the two pi0 photons; missing nu_e
alg_mode1.signal_side do |s|
  s.photons 2                                        # pi0 -> gamma gamma
  s.charged(prp: 1, pim: 1, ep: 1, at_least: true)   # p, pi-, e+ (extra charged tracks allowed)
  s.missing :nu_e                                    # missing nu_e (massless)
  s.min_photon_angle 10.0                            # photon angle to nearest track > 10 degrees
  s.min_photon_energy 0.025                          # photon energy > 0.025 GeV
end

# 4C kinematic fit with chi2 < 200 and the M(p pi0) ~ M(Sigma+) mass constraint
alg_mode1.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:prp, :gamma, :gamma).between(1.176, 1.200)  # M(p pi0) in [1.176, 1.200] GeV/c^2
  f.chi2_cut 200
end

alg_mode1
  .note(:tag_candidate_ranking, "when more than one single-tag candidate is found, retain the candidate with minimum |deltaE|")
  .note(:event_net_charge, "event-level net charge required to be zero; the tagged anti-Lambda_c- contributes -1 while the signal side (p, pi-, e+) contributes +1")
  .note(:energy_scan, "seven sqrt(s) points 4.600-4.699 GeV share a single nominal ECMS constant; the tag kinematic fit reads the per-run measured beam energy from the conditions DB (dtag_reconstruction beam_energy :db)")
  .apply

alg_mode1.execute_on(data_points + incMC_points + exMC_mode1)

# ============================================================================
# Event selection — Modes 2+3: Lambda_c+ -> Sigma+/- pi∓ e+ nu_e,
#                               Sigma+/- -> n pi+/-
# ============================================================================
alg_name_mode23 = "LambdacSemilepSigmaN"
alg_mode23 = TagAnalysis.new(alg_name_mode23)
alg_mode23.set_header(["#{alg_name_mode23}Alg/#{alg_name_mode23}.h"])
          .set_constant({ "ECMS" => [:double, 4.600] })
          .with_decay_card(decay_card_mode23)

# Tag side: same 12 hadronic single-tag modes, both charges scanned
alg_mode23.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKsP, :LambdacPtoKPiP, :LambdacPtoKsPi0P, :LambdacPtoKsPiPiP,
          :LambdacPtoKPiPi0P, :LambdacPtoLambdaPi, :LambdacPtoLambdaPiPi0,
          :LambdacPtoLambdaPiPiPi, :LambdacPtoSigma0Pi, :LambdacPtoSigmaPi0,
          :LambdacPtoSigmaPiPi, :LambdacPtoPiPiP
end

# Signal side: zero photons, at least one pi- and one e+, missing nu_e
alg_mode23.signal_side do |s|
  s.photons 0
  s.charged(pim: 1, ep: 1, at_least: true)   # pi-, e+ (extra charged tracks allowed)
  s.missing :nu_e                            # missing nu_e (massless)
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
end

# 4C kinematic fit with chi2 < 200
alg_mode23.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg_mode23
  .note(:tag_candidate_ranking, "when more than one single-tag candidate is found, retain the candidate with minimum |deltaE|")
  .note(:neutron_detection, "the neutron from Sigma -> n pi is detected via the EMC neutron shower; it is not a charged track and its momentum is taken from the shower, so it is not part of the 4C fit")
  .note(:pi0_veto, "neutron-mode pi0 veto from EMC shower shape: E3x3/E5x5 > 0.9, second moment < 20 cm^2, M(gamma gamma) in [0.115, 0.150] GeV/c^2, chi2(1C) < 20")
  .note(:efficiency_curve, "neutron selection windows M(n pi+) in [1.15, 1.23] GeV/c^2 and M(n pi-) in [1.16, 1.24] GeV/c^2 applied on the reconstructed neutron")
  .note(:background_veto, "further suppression: M(p pi-) > 1.13 GeV/c^2, M(Sigma+ pi- pi(e)+) < 2.27 GeV/c^2, cos theta(pi, e) < 0.95")
  .note(:energy_scan, "seven sqrt(s) points 4.600-4.699 GeV share a single nominal ECMS constant; the tag kinematic fit reads the per-run measured beam energy from the conditions DB (dtag_reconstruction beam_energy :db)")
  .apply

alg_mode23.execute_on(data_points + incMC_points + exMC_mode23)