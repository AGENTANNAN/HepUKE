# DSL for paper 2406.05827v2: Measurement of the integrated luminosity of
# the data collected at 3.773 GeV by BESIII from 2021 to 2024
# Large-angle Bhabha scattering e+e- → e+e-; three data periods, total ~20.3 fb-1
# Babayaga@NLO generator; run-by-run calibrated E_cm

# Data at 3.773 GeV: ψ(3770) sample (BOSS 712)
psipp_data = DatasetManager.real_data.find("712_3773")
psipp_incMC = DatasetManager.inclusive_mc.find("712_3773")

# Decay card for Bhabha signal MC: e+e- → e+e- (large-angle)
# Using KKMC generator with psi(4260) top mother (BESIII convention)
decay_card_bhabha = <<~DECAYCARD
    Decay psi(4260)
    1.000 e+ e- PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for Bhabha events
exMC_bhabha = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "bhabha_luminosity"
  config.related_dataset = psipp_data
  config.events = 500_000
  config.decay_card = decay_card_bhabha
  config.cross_section = :default
end

alg = Algorithm.new("BhabhaLuminosity")
alg.set_header(["BhabhaLuminosityAlg/BhabhaLuminosity.h"])
   .set_constant({"ECMS" => [:double, 3.773]})
   .note(:luminosity_measurement, "Integrated luminosity measured via large-angle Bhabha scattering e+e- → e+e-; Babayaga@NLO generator used for theoretical cross section; three data periods: DATA I (4.995 fb-1), DATA II (8.157 fb-1), DATA III (4.191 fb-1); total 20.275 fb-1")
   .note(:babayaga_generator, "Babayaga@NLO event generator used for Bhabha cross-section calculation with 0.5% precision; run-by-run calibrated E_cm values used")
   .note(:bhabha_event_selection, "Large-angle Bhabha events selected with both e+ and e- in barrel EMC region |cosθ| < 0.8; back-to-back topology; E/p > 0.8 for electron identification; 2-track events with net charge zero")
   .note(:run_by_run_calibration, "Center-of-mass energy calibrated run-by-run using dimuon events e+e- → μ+μ-; E_cm variations accounted for in luminosity calculation")
   .note(:systematic_uncertainties, "Systematic uncertainties from event selection (0.3%), background subtraction (0.1%), tracking/EMC efficiency (0.3%), radiative correction (0.2%), beam energy (0.1%), and Babayaga theory (0.5%); total 0.7%")

sel = Selection.new

# Charged tracks: exactly 2, net charge zero, e+e- topology
sel.select_track do
  cos_theta  0.93
  Vz         10.0
  Vr         1.0
  nChrp      ">=1"
  nChrn      ">=1"
  nChrg      "==2"
  nNet       "==0"
end

# No photon requirement (vetoes radiative Bhabha events with hard photons)
sel.select_photon do
  tdc_emc_start     0
  tdc_emc_end       14
  energyThreshold_b 0.025
  energyThreshold_e 0.050
  angle_to_track    10.0
end

# Lepton identification: both tracks identified as electrons
sel.pid(method: :probability) do
  prob_cut 0.001
  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.5,
                                  treat_as_electron_if_energy_over_momentum_above: 0.8
  nlp ">=1"
  nlm ">=1"
end

# 4C kinematic fit: e+e- constrained to CMS energy
sel.kinematic_fit([:lp, :lm]) do
  nominal
  constrain_four_momentum
  chi2_cut 200  # loose; tight cut applied in ROOT
end

alg.with_decay_card(decay_card_bhabha).apply(sel)
alg.execute_on([psipp_data, psipp_incMC, exMC_bhabha])