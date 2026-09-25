### Dataset preparation ###
# Real data at the eight e+e- centre-of-mass energies of the scan (7.33 fb^-1 total)
data_4128 = DatasetManager.real_data.find("705_4130")   # 4128.48 MeV
data_4157 = DatasetManager.real_data.find("705_4160")   # 4157.83 MeV
data_4178 = DatasetManager.real_data.find("703_4180")   # 4178 MeV
data_4190 = DatasetManager.real_data.find("703_4190")   # 4188.8 MeV
data_4200 = DatasetManager.real_data.find("703_4200")   # 4198.9 MeV
data_4210 = DatasetManager.real_data.find("703_4210")   # 4209.2 MeV
data_4220 = DatasetManager.real_data.find("703_4220")   # 4218.7 MeV
data_4230 = DatasetManager.real_data.find("703_4230")   # 4226.26 MeV

data_points = [data_4128, data_4157, data_4178, data_4190,
               data_4200, data_4210, data_4220, data_4230]

# Matching inclusive MC samples at the same energy points
incMC_4128 = DatasetManager.inclusive_mc.find("705_4130")
incMC_4157 = DatasetManager.inclusive_mc.find("705_4160")
incMC_4178 = DatasetManager.inclusive_mc.find("703_4180")
incMC_4190 = DatasetManager.inclusive_mc.find("703_4190")
incMC_4200 = DatasetManager.inclusive_mc.find("703_4200")
incMC_4210 = DatasetManager.inclusive_mc.find("703_4210")
incMC_4220 = DatasetManager.inclusive_mc.find("703_4220")
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")

incMC_points = [incMC_4128, incMC_4157, incMC_4178, incMC_4190,
                incMC_4200, incMC_4210, incMC_4220, incMC_4230]

# Decay card for e+e- -> D_s*+ D_s-, D_s*+ -> gamma D_s+, D_s+ -> gamma K*+, K*+ -> K+ pi0
# Tag (recoil) decay D_s- -> anti-K0 K-, with pi0 -> gamma gamma and K_S0 -> pi+ pi-
decay_card_kpi0 = <<~DECAYCARD
    Decay psi(4260)
    1.0000 D_s*+ D_s-              PHSP;
    Enddecay

    Decay D_s*+
    1.0000 gamma D_s+              PHSP;
    Enddecay

    Decay D_s+
    1.0000 gamma K*+               PHSP;
    Enddecay

    Decay K*+
    1.0000 K+ pi0                  PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma             PHSP;
    Enddecay

    Decay D_s-
    1.0000 anti-K0 K-              PHSP;
    Enddecay

    Decay anti-K0
    1.0000 K_S0                    PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi-                 PHSP;
    Enddecay

    End
DECAYCARD

# Decay card for the K_S0 pi+ sub-mode of K*+
decay_card_kspi = <<~DECAYCARD
    Decay psi(4260)
    1.0000 D_s*+ D_s-              PHSP;
    Enddecay

    Decay D_s*+
    1.0000 gamma D_s+              PHSP;
    Enddecay

    Decay D_s+
    1.0000 gamma K*+               PHSP;
    Enddecay

    Decay K*+
    1.0000 K_S0 pi+                PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi-                 PHSP;
    Enddecay

    Decay D_s-
    1.0000 anti-K0 K-              PHSP;
    Enddecay

    Decay anti-K0
    1.0000 K_S0                    PHSP;
    Enddecay

    End
DECAYCARD

# Signal exclusive MC: 150k events per K*+ sub-mode, generated at every scan point
# (same signal MC running over the distinct energy points of the scan)
exMCs_kpi0 = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_dsstar_ds_gamma_kpi0"
  config.events        = 150_000
  config.decay_card    = decay_card_kpi0
  config.cross_section = :default
end

exMCs_kspi = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_dsstar_ds_gamma_kspi"
  config.events        = 150_000
  config.decay_card    = decay_card_kspi
  config.cross_section = :default
end

### Event selection (BOSS) — tag-based (D_s tag / recoil D_s- side) ###

## ---------------------------------------------------------------------
## Channel I: D_s*+ -> gamma D_s+, D_s+ -> gamma K*+, K*+ -> K+ pi0
## ---------------------------------------------------------------------
alg_kpi0 = TagAnalysis.new("DsStarDsGamKPi0")
alg_kpi0.set_header(["DsStarDsGamKPi0Alg/DsStarDsGamKPi0.h"])
        .set_constant({"ECMS" => [:double, 4.178]})
        .set_alias({"std::vector<double>" => "Vdouble"})

# Tag side: the recoiling D_s- (charm -1) in the hadronic modes K_S0 K, K K pi, K_S0 K pi pi
alg_kpi0.tag_side(:Ds) do |t|
  t.modes :DstoKsK, :DstoKKPi, :DstoKsKPiPi
  t.charm -1
end

# Signal side: three good showers (radiative gamma and the pi0 -> gamma gamma pair),
# exactly one K+ from K*+ -> K+ pi0
alg_kpi0.signal_side do |s|
  s.photons 3                    # at least three photons
  s.min_photon_angle 10.0        # min angle to charged tracks in degrees
  s.min_photon_energy 0.025      # E_gamma > 25 MeV
  s.charged(kp: 1)               # exactly one K+
  s.require_charge 1
end

# 4C kinematic fit: total four-momentum to the measured CMS four-vector,
# with M(gamma gamma) constrained to the nominal pi0 mass; chi2 < 200
alg_kpi0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.chi2_cut 200
end

# BOSS-side procedures that have no dedicated DSL primitive
alg_kpi0
  .note(:radiative_photon_assignment,
        "the most energetic signal photon (E_gamma > 0.55 GeV) is taken as the radiative
         gamma from D_s*+ -> gamma D_s+; the requirement M(K+ pi0) in 0.83-0.94 GeV/c^2 is
         applied downstream on the stored four-momenta")
  .note(:background_veto,
        "extra photons not used in the 4C fit are vetoed when any gamma-gamma combination
         falls in the pi0 (0.115-0.150 GeV/c^2) or eta (0.50-0.57 GeV/c^2) mass window,
         and M(gamma gamma_h) > 0.62 GeV/c^2 is required; no tag-layer primitive for
         extra-shower combination vetoes")
  .note(:tag_recoil_window,
        "the single-tag recoil mass against the D_s- tag is stored for each energy point and
         windowed per E_cm (2.040-2.220 GeV/c^2) in the ROOT analysis; no tag-side cut is
         applied in BOSS (store-not-cut)")

alg_kpi0.with_decay_card(decay_card_kpi0).apply     # tag layer: apply takes no Selection
alg_kpi0.execute_on(data_points + incMC_points + exMCs_kpi0)

## ---------------------------------------------------------------------
## Channel II: D_s*+ -> gamma D_s+, D_s+ -> gamma K*+, K*+ -> K_S0 pi+
## ---------------------------------------------------------------------
alg_kspi = TagAnalysis.new("DsStarDsGamKsPi")
alg_kspi.set_header(["DsStarDsGamKsPiAlg/DsStarDsGamKsPi.h"])
        .set_constant({"ECMS" => [:double, 4.178]})
        .set_alias({"std::vector<double>" => "Vdouble"})

# Same tag side
alg_kspi.tag_side(:Ds) do |t|
  t.modes :DstoKsK, :DstoKKPi, :DstoKsKPiPi
  t.charm -1
end

# Signal side: one good shower (radiative gamma); two pi+ and one pi- from K_S0 -> pi+ pi-
# and from K*+ -> K_S0 pi+
alg_kspi.signal_side do |s|
  s.photons 1                    # at least one photon
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025      # E_gamma > 25 MeV
  s.charged(pip: 2, pim: 1)      # exactly two pi+ and one pi-
  s.require_charge 1
end

# Same 4C fit, with M(pi+ pi-) constrained to the nominal K_S0 mass; chi2 < 200
alg_kspi.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:pip, :pim).constrain_to_nominal_mass_of(:K_S0)
  f.chi2_cut 200
end

# BOSS-side procedures that have no dedicated DSL primitive
alg_kspi
  .note(:ks0_reconstruction,
        "the signal-side K_S0 is reconstructed from two oppositely charged tracks assigned as
         pi+ pi- without PID, requiring |Vz| < 20 cm, secondary-vertex chi2 < 100,
         0.487 < M(pi+ pi-) < 0.511 GeV/c^2 and a decay length > 2 sigma from the IP; the tag
         signal-side layer carries no secondary-vertex-fit primitive")
  .note(:radiative_photon_assignment,
        "the most energetic signal photon (E_gamma > 0.55 GeV) is taken as the radiative gamma
         from D_s*+ -> gamma D_s+; M(K_S0 pi+) in 0.83-0.94 GeV/c^2 is required downstream")
  .note(:background_veto,
        "extra photons not used in the 4C fit are vetoed when any gamma-gamma combination falls
         in the pi0 (0.115-0.150 GeV/c^2) or eta (0.50-0.57 GeV/c^2) mass window")
  .note(:tag_recoil_window,
        "the single-tag recoil mass against the D_s- tag is stored for each energy point and
         windowed per E_cm (2.040-2.220 GeV/c^2) in the ROOT analysis; no tag-side cut is
         applied in BOSS (store-not-cut)")

alg_kspi.with_decay_card(decay_card_kspi).apply
alg_kspi.execute_on(data_points + incMC_points + exMCs_kspi)