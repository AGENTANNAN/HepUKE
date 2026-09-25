# =============================================================================
#  Lambda_c+ semileptonic decays with a single-tag technique
#    mode 1 : Lambda_c+ -> Lambda pi+ pi- e+ nu_e   (Lambda -> p pi-)
#    mode 2 : Lambda_c+ -> p K_S0 pi- e+ nu_e       (K_S0  -> pi+ pi-)
#  Tag : anti-Lambda_c- (charm flavour -1) reconstructed from nine ST modes
#  Data : 4.600, 4.612, 4.628, 4.641, 4.661, 4.682, 4.699 GeV
# =============================================================================

### Dataset preparation ###

# --- seven real-data samples (BOSS_version_CMSenergy) ---
data_4600 = DatasetManager.real_data.find("703_4600")   # 4599.53 MeV
data_4612 = DatasetManager.real_data.find("706_4610")   # 4611.86 MeV -> closest sample to 4.612 GeV
data_4628 = DatasetManager.real_data.find("706_4620")   # 4628.00 MeV
data_4641 = DatasetManager.real_data.find("706_4640")   # 4640.91 MeV
data_4661 = DatasetManager.real_data.find("706_4660")   # 4661.24 MeV
data_4682 = DatasetManager.real_data.find("706_4680")   # 4681.92 MeV
data_4699 = DatasetManager.real_data.find("706_4700")   # 4698.82 MeV

# --- inclusive MC at six of the seven points ---
incmc_4600 = DatasetManager.inclusive_mc.find("703_4600")
incmc_4628 = DatasetManager.inclusive_mc.find("706_4620")
incmc_4641 = DatasetManager.inclusive_mc.find("706_4640")
incmc_4661 = DatasetManager.inclusive_mc.find("706_4660")
incmc_4682 = DatasetManager.inclusive_mc.find("706_4680")
incmc_4699 = DatasetManager.inclusive_mc.find("706_4700")
# the 4.612 GeV inclusive MC is absent from the standard table -> nearest sample stands in
incmc_4612 = incmc_4600   # closest available inclusive MC used for the 4.612 GeV point

### Decay cards and exclusive MC ###

# mode 1 : Lambda_c+ -> Lambda pi+ pi- e+ nu_e ; tag side anti-Lambda_c- -> pbar K_S0
decay_card_mode1 = <<~DECAYCARD
  Decay psi(4260)
  1.0 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay

  Decay Lambda_c+
  1.0 Lambda pi+ pi- e+ nu_e PHSP;
  Enddecay

  Decay Lambda
  1.0 p+ pi- HypWK;
  Enddecay

  Decay anti-Lambda_c-
  1.0 anti-p- K_S0 PHSP;
  Enddecay

  Decay K_S0
  1.0 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# mode 2 : Lambda_c+ -> p K_S0 pi- e+ nu_e ; tag side anti-Lambda_c- -> pbar K_S0
decay_card_mode2 = <<~DECAYCARD
  Decay psi(4260)
  1.0 Lambda_c+ anti-Lambda_c- PHSP;
  Enddecay

  Decay Lambda_c+
  1.0 p+ K_S0 pi- e+ nu_e PHSP;
  Enddecay

  Decay K_S0
  1.0 pi+ pi- PHSP;
  Enddecay

  Decay anti-Lambda_c-
  1.0 anti-p- K_S0 PHSP;
  Enddecay

  End
DECAYCARD

# 200k-event exclusive MC for each of the two signal modes
exMC_mode1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_Lc_semilep_mode1"
  config.related_dataset = data_4600
  config.events         = 200_000
  config.decay_card     = decay_card_mode1
  config.cross_section  = :default
end

exMC_mode2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_Lc_semilep_mode2"
  config.related_dataset = data_4600
  config.events         = 200_000
  config.decay_card     = decay_card_mode2
  config.cross_section  = :default
end

# common dataset list processed by both signal-mode algorithms
datasets_common = [data_4600, data_4612, data_4628, data_4641, data_4661, data_4682, data_4699,
                   incmc_4600, incmc_4628, incmc_4641, incmc_4661, incmc_4682, incmc_4699]

# =============================================================================
#  Mode 1 : Lambda_c+ -> Lambda pi+ pi- e+ nu_e
# =============================================================================
alg_mode1 = TagAnalysis.new("LcSemilepMode1")
alg_mode1.set_header(["LcSemilepMode1Alg/LcSemilepMode1Alg.h"])
         .set_constant({"ECMS" => [:double, 4.600]})   # nominal beam energy (per-run energy taken from MeasuredEcmsSvc)
         .with_decay_card(decay_card_mode1)

# tag side : anti-Lambda_c- single tag from nine ST modes, charm flavour -1
alg_mode1.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKsP,          # pbar K_S0
          :LambdacPtoKPiP,         # K+ pi-  (pbar K+ pi-)
          :LambdacPtoKsPi0P,       # K_S0 pi0
          :LambdacPtoKsPiPiP,      # K_S0 pi+ pi-
          :LambdacPtoKPiPi0P,      # K+ pi- pi0
          :LambdacPtoLambdaPi,     # Lambdabar pi-
          :LambdacPtoLambdaPiPi0,  # Lambdabar pi- pi0
          :LambdacPtoLambdaPiPiPi, # Lambdabar pi- pi+ pi-
          :LambdacPtoPiPiPiP       # pi+ pi+ pi- type
  t.charm(-1)
end

# signal side : the five tracks the tag did not use + missing massless nu_e
alg_mode1.signal_side do |s|
  s.charged(prp: 1, pip: 1, pim: 2, ep: 1)   # 1 p, 1 pi+, 2 pi-, 1 e+  (no photons used)
  s.require_charge(1)                        # net charge +1
  s.missing :nu_e                            # one missing massless nu_e (semileptonic tag)
end

# 4C kinematic fit constraining the tag to the nominal Lambda_c mass
alg_mode1.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.chi2_cut 200
end

alg_mode1
  .note(:track_quality_selection,
        "Signal-side tracks required tight: |Vz| < 10 cm, Vr < 1 cm, |cos(theta)| < 0.93. "
        "Lambda daughter tracks are looser: |Vz| < 20 cm and no Vr cut.")
  .note(:secondary_vertex_fit,
        "Lambda reconstructed from a secondary-vertex fit (chi2 < 100, positive decay length); "
        "mass window M(p pi-) in [1.09, 1.14] GeV.")
  .note(:pid_correction_method,
        "PID by likelihood ratios: p requires L(p)>L(K), L(p)>L(pi), L(p)>0; pi requires "
        "L(pi)>L(K), L(pi)>0; e+ requires L(e)>0.001 and L(e)/(L(e)+L(pi)+L(K))>0.99. "
        "FSR showers inside a 5 degree cone around the positron are added back.")
  .note(:background_veto,
        "cos(theta_{e,pi}) < 0.88; M(Lambda pi+ pi- e) < 2.20 GeV; "
        "cos(theta_{p_miss,gamma}) < 0.82; U_miss = E_miss - |p_miss| in [-0.08, 0.08] GeV.")
  .note(:efficiency_curve,
        "Seven energy points (4.600-4.699 GeV) are combined; the 4C fit uses the per-run measured "
        "CMS four-vector (MeasuredEcmsSvc) and ECMS = 4.600 GeV is only the nominal constant.")

alg_mode1.apply                                   # tag analysis apply takes no Selection argument
alg_mode1.execute_on(datasets_common + [exMC_mode1])

# =============================================================================
#  Mode 2 : Lambda_c+ -> p K_S0 pi- e+ nu_e
# =============================================================================
alg_mode2 = TagAnalysis.new("LcSemilepMode2")
alg_mode2.set_header(["LcSemilepMode2Alg/LcSemilepMode2Alg.h"])
         .set_constant({"ECMS" => [:double, 4.600]})   # nominal beam energy (per-run energy taken from MeasuredEcmsSvc)
         .with_decay_card(decay_card_mode2)

# tag side : identical anti-Lambda_c- single-tag chain as mode 1
alg_mode2.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKsP,          # pbar K_S0
          :LambdacPtoKPiP,         # K+ pi-  (pbar K+ pi-)
          :LambdacPtoKsPi0P,       # K_S0 pi0
          :LambdacPtoKsPiPiP,      # K_S0 pi+ pi-
          :LambdacPtoKPiPi0P,      # K+ pi- pi0
          :LambdacPtoLambdaPi,     # Lambdabar pi-
          :LambdacPtoLambdaPiPi0,  # Lambdabar pi- pi0
          :LambdacPtoLambdaPiPiPi, # Lambdabar pi- pi+ pi-
          :LambdacPtoPiPiPiP       # pi+ pi+ pi- type
  t.charm(-1)
end

# signal side : identical five-track content (p, pi+, 2 pi-, e+) + missing nu_e
alg_mode2.signal_side do |s|
  s.charged(prp: 1, pip: 1, pim: 2, ep: 1)   # 1 p, 1 pi+, 2 pi-, 1 e+  (no photons used)
  s.require_charge(1)                        # net charge +1
  s.missing :nu_e                            # one missing massless nu_e
end

# same shared 4C kinematic fit constraining the tag to the nominal Lambda_c mass
alg_mode2.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.chi2_cut 200
end

alg_mode2
  .note(:track_quality_selection,
        "Signal-side tracks required tight: |Vz| < 10 cm, Vr < 1 cm, |cos(theta)| < 0.93. "
        "K_S0 daughter tracks are looser: |Vz| < 20 cm and no Vr cut.")
  .note(:secondary_vertex_fit,
        "K_S0 reconstructed from a secondary-vertex fit (chi2 < 100, positive decay length); "
        "mass window M(pi+ pi-) in [0.490, 0.504] GeV.")
  .note(:pid_correction_method,
        "PID by likelihood ratios: p requires L(p)>L(K), L(p)>L(pi), L(p)>0; pi requires "
        "L(pi)>L(K), L(pi)>0; e+ requires L(e)>0.001 and L(e)/(L(e)+L(pi)+L(K))>0.98 "
        "(loosened from 0.99 in mode 1). FSR showers inside a 5 degree cone around the positron "
        "are added back.")
  .note(:background_veto,
        "cos(theta_{e,pi}) < 0.92; M(p K_S0 pi- e) < 2.28 GeV; "
        "cos(theta_{p_miss,gamma}) < 0.90; U_miss = E_miss - |p_miss| in [-0.08, 0.08] GeV.")
  .note(:efficiency_curve,
        "Seven energy points (4.600-4.699 GeV) are combined; the 4C fit uses the per-run measured "
        "CMS four-vector (MeasuredEcmsSvc) and ECMS = 4.600 GeV is only the nominal constant.")

alg_mode2.apply                                   # tag analysis apply takes no Selection argument
alg_mode2.execute_on(datasets_common + [exMC_mode2])