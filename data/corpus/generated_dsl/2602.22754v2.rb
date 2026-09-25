# =====================================================================
# Λc+ → Σ0 K_S0 π+  and  Λc+ → Σ0 K_S0 K+   —   single-tag (ST) analysis
# 13 XYZ energy points, 4.600 – 4.950 GeV
# =====================================================================

### Datasets ###
# BESIII sample-name convention: [BOSS version]_[CMS energy in MeV]
energy_point_names = %w[
  703_4600
  706_4610 706_4620 706_4640 706_4660 706_4680 706_4700
  707_4740 707_4750 707_4780 707_4840 707_4914 707_4946
]
# Real data at each energy point
data_points = energy_point_names.map { |n| DatasetManager.real_data.find(n) }
# Matching inclusive MC at each energy point
incmc_points = energy_point_names.map { |n| DatasetManager.inclusive_mc.find(n) }

### Decay cards (EvtGen, phase-space) ###
# Λc+ → Σ0 K_S0 π+ ,  Σ0 → γΛ , Λ → pπ− , K_S0 → π+π−
decay_card_pi = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Lambda_c+ anti-Lambda_c-     PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000 Sigma0 K_S0 pi+              PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0000 anti-Sigma0 K_S0 pi-         PHSP;
    Enddecay

    Decay Sigma0
    1.0000 gamma Lambda0                PHSP;
    Enddecay

    Decay anti-Sigma0
    1.0000 gamma anti-Lambda0           PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-                       PHSP;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+                  PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi-                      PHSP;
    Enddecay

    End
DECAYCARD

# Λc+ → Σ0 K_S0 K+ ,  Σ0 → γΛ , Λ → pπ− , K_S0 → π+π−
decay_card_k = <<~DECAYCARD
    Decay psi(4260)
    1.0000 Lambda_c+ anti-Lambda_c-     PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000 Sigma0 K_S0 K+               PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0000 anti-Sigma0 K_S0 K-          PHSP;
    Enddecay

    Decay Sigma0
    1.0000 gamma Lambda0                PHSP;
    Enddecay

    Decay anti-Sigma0
    1.0000 gamma anti-Lambda0           PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-                       PHSP;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+                  PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi-                      PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive signal MC: 200k events per mode per energy point ###
exMC_pi = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_lctosigma0kspi"   # auto-suffixed per energy point
  config.events        = 200_000
  config.decay_card    = decay_card_pi
  config.cross_section = :default                # default production cross section
end

exMC_k = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_lctosigma0ksk"    # auto-suffixed per energy point
  config.events        = 200_000
  config.decay_card    = decay_card_k
  config.cross_section = :default                # default production cross section
end

# =====================================================================
# TagAnalysis — Λc+ reconstructed as a single tag (DTagTool);
# the recoil anti-Λc− is the missing particle on the signal side.
# =====================================================================

### Mode I : Λc+ → Σ0 K_S0 π+ ###
alg_name_pi = "LcSTSigma0KsPi"
alg_pi = TagAnalysis.new(alg_name_pi)
alg_pi.set_header(["#{alg_name_pi}Alg/#{alg_name_pi}.h"])
      .set_constant({ "ECMS" => [:double, 4.600] })   # 13-point scan: per-run measured beam energy is used by the fit
      .note(:tag_side_reconstruction,
            "all tag-side reconstruction is performed inside DTagAlg and is not re-expressed " \
            "by the tag DSL: charged tracks with |cosθ|<0.93, |Vz|<10 cm, Vr<1 cm, >=2 positive " \
            "and >=2 negative tracks, net charge 0; photons with E>25 MeV (barrel)/>50 MeV " \
            "(endcap), TDC in [0,700] ns and angle to the nearest track >10 deg; PID by the " \
            "probability method (prob>0.001) separating p, K and pi (>=1 p, >=1 pi+, >=1 pi-); " \
            "Λ→pπ- and K_S0→π+π- secondary-vertex fits minimising the mass difference with " \
            "M(pπ-) in [1.111,1.121] GeV/c2, M(π+π-) in [0.487,0.511] GeV/c2, daughter-track " \
            "|Vz|<20 cm and vertex chi2<100; Σ0 formed from γΛ by a Kalman fit constraining " \
            "M(γΛ) to the Σ0 nominal mass (chi2<200).")
      .with_decay_card(decay_card_pi)

# Tag side: single tag Λc+ → Σ0 K_S0 π+
alg_pi.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoSigma0KsPiP
  t.charm 1                      # tag the Λc+ (not the anti-Λc-)
end

# Signal side: only the recoil anti-Λc- is left over (missing, massive)
alg_pi.signal_side do |s|
  s.missing :"anti-Lambda_c-", mass: 2.28646
end

# Kinematic fit: 4C against √s + Λc+ mass constraint on the tag group;
# the missing anti-Λc- recoil is mass-constrained by the missing-particle treatment.
# chi2<200 at reconstruction — the optimised χ²_4C<29 is applied in ROOT,
# inside the M_BC signal window [2.282, 2.291] GeV/c² (mBC is stored, not cut here).
alg_pi.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.chi2_cut 200
end

alg_pi.apply
root_files_pi = alg_pi.execute_on(data_points + incmc_points + exMC_pi)

### Mode II : Λc+ → Σ0 K_S0 K+ ###
alg_name_k = "LcSTSigma0KsK"
alg_k = TagAnalysis.new(alg_name_k)
alg_k.set_header(["#{alg_name_k}Alg/#{alg_name_k}.h"])
     .set_constant({ "ECMS" => [:double, 4.600] })   # 13-point scan: per-run measured beam energy is used by the fit
     .note(:tag_side_reconstruction,
           "all tag-side reconstruction is performed inside DTagAlg and is not re-expressed " \
           "by the tag DSL: charged tracks with |cosθ|<0.93, |Vz|<10 cm, Vr<1 cm, >=2 positive " \
           "and >=2 negative tracks, net charge 0; photons with E>25 MeV (barrel)/>50 MeV " \
           "(endcap), TDC in [0,700] ns and angle to the nearest track >10 deg; PID by the " \
           "probability method (prob>0.001) separating p, K and pi (>=1 p, >=1 pi+, >=1 pi-, " \
           ">=1 K+); Λ→pπ- and K_S0→π+π- secondary-vertex fits minimising the mass difference " \
           "with M(pπ-) in [1.111,1.121] GeV/c2, M(π+π-) in [0.487,0.511] GeV/c2, daughter-track " \
           "|Vz|<20 cm and vertex chi2<100; Σ0 formed from γΛ by a Kalman fit constraining " \
           "M(γΛ) to the Σ0 nominal mass (chi2<200).")
     .with_decay_card(decay_card_k)

# Tag side: single tag Λc+ → Σ0 K_S0 K+
alg_k.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoSigma0KsKP
  t.charm 1                      # tag the Λc+ (not the anti-Λc-)
end

# Signal side: only the recoil anti-Λc- is left over (missing, massive)
alg_k.signal_side do |s|
  s.missing :"anti-Lambda_c-", mass: 2.28646
end

# Kinematic fit: 4C against √s + Λc+ mass constraint on the tag group;
# chi2<200 at reconstruction — the optimised χ²_4C<171 is applied in ROOT,
# inside the M_BC signal window [2.282, 2.291] GeV/c² (mBC is stored, not cut here).
alg_k.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.chi2_cut 200
end

alg_k.apply
root_files_k = alg_k.execute_on(data_points + incmc_points + exMC_k)