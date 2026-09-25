# Paper: 2406.18083v2
# KS0-KL0 asymmetries in Λc+ → pKL0, pKL0π+π-, pKL0π0
# DT tag-based: 12 ST tag modes for anti-Λc- → signal Λc+ decays
# Energy range: 4.600-4.699 GeV (BOSS 706), 4.5 fb⁻¹

# Datasets: BOSS 706, 7 energy points
ds_4610 = DatasetManager.real_data.find("706_4610")
ds_4620 = DatasetManager.real_data.find("706_4620")
ds_4640 = DatasetManager.real_data.find("706_4640")
ds_4660 = DatasetManager.real_data.find("706_4660")
ds_4680 = DatasetManager.real_data.find("706_4680")
ds_4700 = DatasetManager.real_data.find("706_4700")
# Note: 4.612 GeV corresponds to dataset entry 706_4610 (4611.86 MeV);
# 4.682 GeV corresponds to dataset entry 706_4680 (4681.92 MeV);
# the paper groups: 4600→4610, 4612→4610, 4628→4620, 4641→4640,
# 4661→4660, 4682→4680, 4698→4700

data_points = [ds_4610, ds_4620, ds_4640, ds_4660, ds_4680, ds_4700]

incMC_4610 = DatasetManager.inclusive_mc.find("706_4610")
incMC_4620 = DatasetManager.inclusive_mc.find("706_4620")
incMC_4640 = DatasetManager.inclusive_mc.find("706_4640")
incMC_4660 = DatasetManager.inclusive_mc.find("706_4660")
incMC_4680 = DatasetManager.inclusive_mc.find("706_4680")
incMC_4700 = DatasetManager.inclusive_mc.find("706_4700")

incMC_points = [incMC_4610, incMC_4620, incMC_4640, incMC_4660, incMC_4680, incMC_4700]

# Tag mode list (anti-Λc- reconstructed via 12 CF modes)
# DTagAlg channel names for Lambda_c+ (charge conjugate of the anti-Λc- tag):
# 1. p K_S0 → :LambdacPtoKsP
# 2. p K- π+ → :LambdacPtoKPiP
# 3. p K_S0 π0 → :LambdacPtoKsPi0P
# 4. p K_S0 π+ π- → :LambdacPtoKsPiPiP
# 5. p K- π+ π0 → :LambdacPtoKPiPi0P
# 6. p π+ π- → :LambdacPtoPiPiP
# 7. Λ π+ → :LambdacPtoLambdaPi
# 8. Λ π+ π0 → :LambdacPtoLambdaPiPi0
# 9. Λ π+ π- π+ → :LambdacPtoLambdaPiPiPi
# 10. Σ0 π+ → :LambdacPtoSigma0Pi
# 11. Σ+ π0 → :LambdacPtoSigmaPPi0
# 12. Σ+ π+ π- → :LambdacPtoSigmaPPiPi
tag_modes = [:LambdacPtoKsP, :LambdacPtoKPiP, :LambdacPtoKsPi0P,
             :LambdacPtoKsPiPiP, :LambdacPtoKPiPi0P, :LambdacPtoPiPiP,
             :LambdacPtoLambdaPi, :LambdacPtoLambdaPiPi0, :LambdacPtoLambdaPiPiPi,
             :LambdacPtoSigma0Pi, :LambdacPtoSigmaPPi0, :LambdacPtoSigmaPPiPi]

# Decay card: e+e- → Λc+ anti-Λc- at threshold
# Signal side: Λc+ → pKL0 (or pKL0π+π- or pKL0π0)
# KL0 treated as missing particle
decay_card = <<~DECAYCARD
    Decay psi(4260)
    1.0000 anti-Lambda_c- Lambda_c+ PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000 p K_L0 PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0000 anti-p- K_S0 PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

exMCs = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "Lc_to_pKL0"
  config.events = 200000
  config.decay_card = decay_card
  config.cross_section = :default
end

# ============================================================
# Analysis A: Λc+ → p KL0
# ============================================================
alg_pKL0 = TagAnalysis.new("LcpKL0_DT")
alg_pKL0.set_header(["LcpKL0Alg/LcpKL0.h"])
         .set_constant({ "ECMS" => [:double, 4.641] })
         .note(:ecms_per_dataset,
           "ECMS per dataset: 4611.86 (4610), 4628.00 (4620), 4640.91 (4640),
            4661.24 (4660), 4681.92 (4680), 4698.82 (4700) MeV;
            per-run MeasuredEcmsSvc handles the actual value")

alg_pKL0.tag_side(:Lambdac) do |t|
  t.modes(*tag_modes)
  t.charm -1
end

alg_pKL0.signal_side do |s|
  s.charged(prp: 1)
  s.require_charge 1
  s.missing :K_L0
end

# 5C fit: 4-momentum + tag Lambda_c mass
# (Signal Lambda_c mass constraint would make it 6C but is handled by
#  constraining the tag + missing + signal system in ROOT)
alg_pKL0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.chi2_cut 200
end

alg_pKL0
  .note(:signal_lc_mass_constraint,
    "6C kinematic fit: 4-momentum + tag Λc mass + signal Λc mass.
     The signal Λc mass constraint (M(prp + KL0) = m_Λc) is not expressible
     in the current DSL; the readback of the 5C (4C + tag mass) chi2 is
     stored, and the full 6C is performed in the ROOT analysis stage.")
  .note(:chi2_optimized,
    "Optimized χ² < 60 for pKL0, < 25 for pKL0π+π-, < 20 for pKL0π0
     (from s/√(s+b) optimisation in ROOT)")
  .note(:background_vetoes,
    "Λ→pπ- veto: 1.11-1.12 GeV/c²; KS0→π+π- veto: 0.48-0.52 GeV/c²;
     M_recoil(p) > 1.0 GeV/c² for pKL0π+π- mode")
  .note(:mc_angular,
    "Λc+→pKL0 angular distribution modeled with decay asymmetry parameters from Ref. [42];
     signal models for pKL0π+π- and pKL0π0 tuned from data;
     inclusive MC ~40x data luminosity")
  .note(:dt_yield_formula,
    "DT yields extracted via simultaneous unbinned ML fit to M²_miss distributions
     at all 7 c.m. energies; KS0X reference mode used for BF determination")
  .with_decay_card(decay_card)
  .apply

alg_pKL0.execute_on(data_points + incMC_points + exMCs)

# ============================================================
# Analysis B: Λc+ → p KL0 π+ π-
# ============================================================
alg_pKL0pipi = TagAnalysis.new("LcpKL0PiPi_DT")
alg_pKL0pipi.set_header(["LcpKL0PiPiAlg/LcpKL0PiPi.h"])
            .set_constant({ "ECMS" => [:double, 4.641] })
            .note(:ecms_per_dataset,
              "ECMS per dataset: 4611.86 (4610), 4628.00 (4620), 4640.91 (4640),
               4661.24 (4660), 4681.92 (4680), 4698.82 (4700) MeV;
               per-run MeasuredEcmsSvc handles the actual value")

alg_pKL0pipi.tag_side(:Lambdac) do |t|
  t.modes(*tag_modes)
  t.charm -1
end

alg_pKL0pipi.signal_side do |s|
  s.charged(prp: 1, pip: 1, pim: 1)
  s.require_charge 1   # +1 (p) + 1 (π+) - 1 (π-) = +1
  s.missing :K_L0
end

alg_pKL0pipi.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.chi2_cut 200
end

alg_pKL0pipi
  .note(:signal_lc_mass_constraint,
    "6C kinematic fit: 4-momentum + tag Λc mass + signal Λc mass.
     Signal Λc = prp + π+ + π- + KL0(missing)")
  .note(:chi2_optimized,
    "Optimized χ² < 25 for this mode (from s/√(s+b) in ROOT)")
  .note(:background_vetoes,
    "Λ→pπ- veto: 1.11-1.12 GeV/c²; KS0→π+π- veto: 0.48-0.52 GeV/c²;
     M_recoil(p) > 1.0 GeV/c²")
  .note(:amplitude_model,
    "Amplitude model includes Σ*, Δ*, N*, K* and ρ intermediate resonances;
     signal models tuned from data")
  .with_decay_card(decay_card)
  .apply

alg_pKL0pipi.execute_on(data_points + incMC_points + exMCs)

# ============================================================
# Analysis C: Λc+ → p KL0 π0
# ============================================================
alg_pKL0pi0 = TagAnalysis.new("LcpKL0Pi0_DT")
alg_pKL0pi0.set_header(["LcpKL0Pi0Alg/LcpKL0Pi0.h"])
           .set_constant({ "ECMS" => [:double, 4.641] })
           .note(:ecms_per_dataset,
             "ECMS per dataset: 4611.86 (4610), 4628.00 (4620), 4640.91 (4640),
              4661.24 (4660), 4681.92 (4680), 4698.82 (4700) MeV;
              per-run MeasuredEcmsSvc handles the actual value")

alg_pKL0pi0.tag_side(:Lambdac) do |t|
  t.modes(*tag_modes)
  t.charm -1
end

alg_pKL0pi0.signal_side do |s|
  s.photons 2
  s.min_photon_angle 10.0
  s.charged(prp: 1)
  s.require_charge 1
  s.missing :K_L0
end

# 6C fit: 4-momentum + tag Λc mass + signal Λc mass + π0 mass
# π0 mass constraint via invariant mass of signal photons
alg_pKL0pi0.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.invariant_mass_of(:tag1).constrain_to_nominal_mass_of(:"Lambda_c+")
  f.chi2_cut 200
end

alg_pKL0pi0
  .note(:signal_lc_mass_constraint,
    "Signal Λc mass constraint (M(prp + π0 + KL0) = m_Λc) handled in ROOT;
     Provides additional 1C beyond the 6C declared here")
  .note(:chi2_optimized,
    "Optimized χ² < 20 for this mode (from s/√(s+b) in ROOT)")
  .note(:background_vetoes,
    "Λ→pπ- veto: 1.11-1.12 GeV/c²; KS0→π+π- veto: 0.48-0.52 GeV/c²;
     M_recoil(p) > 0.65 GeV/c²; Σ+→pπ0 veto: 1.17-1.20 GeV/c²;
     pKS0(→π0π0) background suppressed via recoil mass")
  .note(:photon_isolation,
    "Photons isolated from proton (>20°); π0: M(γγ) ∈ [0.115,0.150] GeV/c²,
     1C fit χ² < 20 (applied in ROOT stage)")
  .note(:truth_matching,
    "Truth-matching method for pure signal shape: opening angle
     between truth and reconstructed photons < 10°")
  .with_decay_card(decay_card)
  .apply

alg_pKL0pi0.execute_on(data_points + incMC_points + exMCs)