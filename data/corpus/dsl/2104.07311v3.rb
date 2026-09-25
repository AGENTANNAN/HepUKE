# DSL for: Inclusive semielectronic decay branching fraction of D_s+
# Paper: 2104.07311v3
# Classification: TagAnalysis (ST + missing, semileptonic tag)
# Data: √s = 4.178–4.230 GeV, three data-set groups
# Process: e+e- → D_s*± D_s∓, tag D_s- → K+K-π-, signal D_s+ → X e+ ν_e

# Three data-set groups (analyzed separately due to detector conditions)
data_4178 = DatasetManager.real_data.find("703_4180")
data_4190 = DatasetManager.real_data.find("703_4190")
data_4200 = DatasetManager.real_data.find("703_4200")
data_4210 = DatasetManager.real_data.find("703_4210")
data_4220 = DatasetManager.real_data.find("703_4220")
data_4230 = DatasetManager.real_data.find("703_4230")

# Group 1: 4.178 GeV
datasets_group1 = [data_4178]

# Group 2: 4.189–4.219 GeV
datasets_group2 = [data_4190, data_4200, data_4210, data_4220]

# Group 3: 4.225–4.230 GeV
datasets_group3 = [data_4230]

all_datasets = datasets_group1 + datasets_group2 + datasets_group3

incMC_4178 = DatasetManager.inclusive_mc.find("703_4180")
incMC_4190 = DatasetManager.inclusive_mc.find("703_4190")
incMC_4200 = DatasetManager.inclusive_mc.find("703_4200")
incMC_4210 = DatasetManager.inclusive_mc.find("703_4210")
incMC_4220 = DatasetManager.inclusive_mc.find("703_4220")
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")

all_incMC = [incMC_4178, incMC_4190, incMC_4200, incMC_4210, incMC_4220, incMC_4230]

# Signal MC: D_s+ → X e+ ν_e inclusive
sig_mc = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_Ds_to_Xepnue_inclusive"
  config.related_dataset = data_4178
  config.events          = 1_000_000
  config.decay_card      = <<~DECAYCARD
    Decay D_s+
    1.0000 e+ nu_e phi PHSP;
    Enddecay
    Decay phi
    1.0000 K+ K- PHSP;
    Enddecay
    End
  DECAYCARD
  config.cross_section   = :default
end

# ===========================================================================
# TagAnalysis: ST + missing (semileptonic tag)
# Single tag D_s- → K+K-π- (only tag mode used)
# Signal side: e+ + missing ν_e
# ===========================================================================

alg = TagAnalysis.new("DsToInclusiveSemielectronic")
alg.set_header(["DsToInclusiveSemielectronicAlg/DsToInclusiveSemielectronic.h"])
  .set_constant({ "ECMS" => [:double, 4.178] })

# Tag side: D_s- reconstructed via K+K-π- (single mode)
alg.tag_side(:Ds) do |t|
  t.modes :DstoKKPi
  t.charm -1
end

# Signal side: positron from D_s+ → X e+ ν_e, missing neutrino
# Positron p > 200 MeV/c (applied in ROOT via momentum binning)
alg.signal_side do |s|
  s.photons 0
  s.charged(ep: 1)
  s.missing :nu_e
  s.require_charge 1
end

# No kinematic fit on signal side — semileptonic with missing neutrino
# Tag-side fit: recoil mass window used for D_s* D_s event selection
# Recoil mass against tag D_s- computed from 4-momentum conservation
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

alg.note(:recoil_mass_windows,
  "Recoil mass M_Rec windows per data set: [2057, 2177] MeV/c² for 4.178 GeV; [2145, 2190] MeV/c² for 4.189–4.219 GeV; [2150, 2200] MeV/c² for 4.225–4.230 GeV. Applied in ROOT analysis.")

alg.note(:tag_invariant_mass_fit,
  "Single-tag yield n_ST determined via unbinned fit to D_s- → K+K-π- invariant mass. Signal shape: MC convolved with Gaussian (width/mean free). Background: second-order Chebyshev polynomial. Fit range: ±40 MeV/c² around D_s mass (4.178 GeV), ±35 MeV/c² (other data sets).")

alg.note(:multiple_tag_candidates,
  "Multiple tag candidates per event allowed to increase efficiency and minimize tag bias b_tag. About 20% of events have >1 candidate after recoil-mass selection. Background from extra candidates subtracted in fitting.")

alg.note(:positron_momentum_cut,
  "Recoil-side positron candidates required p > 200 MeV/c. Sorted into 18 momentum bins of 50 MeV/c from 200 to 1100 MeV/c. Momentum spectrum extrapolated below 200 MeV/c via fit to assumed spectral shape from 6 observed exclusive modes.")

alg.note(:pid_unfolding,
  "Complex PID unfolding methodology: 3×3 matrix inversion (e, π, K) to correct for PID efficiencies and misidentification rates. Right-sign (RS) minus wrong-sign (WS) subtraction to remove charge-symmetric backgrounds. Tracking unfolding matrix corrects for resolution and FSR effects. Fully implemented in ROOT analysis — inexpressible in BOSS DSL.")

alg.note(:tau_subtraction,
  "Contribution of D_s+ → τ+ ν_τ → e+ ν_e ν_τ ν̄_τ subtracted from positron spectrum using known branching fraction B(D_s+ → τ+ ν_τ) and MC-predicted momentum spectrum with tag bias correction b_tag,τ.")

alg.note(:tag_bias,
  "Tag bias b_tag = ε_ST'/ε_ST accounts for difference in tag reconstruction efficiency given semileptonic decay on recoil side. Determined from MC for each observed exclusive semielectronic mode, weighted by branching fractions. b_tag ≈ 1.004–1.007 depending on data set.")

alg.note(:spectrum_extrapolation,
  "Positron yield for p_e ≤ 200 MeV/c extrapolated by fitting MC-based spectral shape to measured yields with p_e > 200 MeV/c. Shape constructed from 6 observed exclusive modes (φ e+ν, η e+ν, η' e+ν, K0 e+ν, K*(892)0 e+ν, f0(980) e+ν) with form-factor parameterizations. Systematic uncertainty includes unobserved modes (h1(1415), f1(1510), γ e+ν).")

alg.note(:branching_fraction_formula,
  "B(D_s+ → X e+ ν_e) = N_DT / (n_ST · b_tag). N_DT = efficiency-corrected double-tag yield summed over p_e bins with extrapolation to p_e ≤ 200 MeV/c. Three data sets analyzed independently and combined. Final result: (6.30 ± 0.13 ± 0.10)%.")

alg.apply
alg.execute_on([data_4178, incMC_4178, sig_mc])