# Measurement of Born cross section and effective form factor
# for e+e- -> n nbar at 18 center-of-mass energies between 2.0 and 3.08 GeV
# BESIII: 647.9 pb-1 total, R-scan data

# --- ConExc decay card for continuum n nbar production ---
# Mode 79 = n nbar (PDG codes 2112 -2112)
decay_card_nnbar = <<~DECAYCARD
    Decay vpho
    1.000 ConExc 79;
    Enddecay
End
DECAYCARD

# --- Datasets (R-scan at 18 energy points) ---
data_2000 = DatasetManager.real_data.find("713_Rscan_2000")
data_2050 = DatasetManager.real_data.find("713_Rscan_2050")
data_2100 = DatasetManager.real_data.find("713_Rscan_2100")
data_2125 = DatasetManager.real_data.find("713_Rscan_2125")
data_2150 = DatasetManager.real_data.find("713_Rscan_2150")
data_2175 = DatasetManager.real_data.find("713_Rscan_2175")
data_2200 = DatasetManager.real_data.find("713_Rscan_2200")
data_2232 = DatasetManager.real_data.find("713_Rscan_2232")
data_2309 = DatasetManager.real_data.find("713_Rscan_2309")
data_2386 = DatasetManager.real_data.find("713_Rscan_2386")
data_2396 = DatasetManager.real_data.find("713_Rscan_2396")
data_2645 = DatasetManager.real_data.find("713_Rscan_2645")
data_2900 = DatasetManager.real_data.find("713_Rscan_2900")
data_2950 = DatasetManager.real_data.find("713_Rscan_2950")
data_2981 = DatasetManager.real_data.find("713_Rscan_2981")
data_3000 = DatasetManager.real_data.find("713_Rscan_3000")
data_3020 = DatasetManager.real_data.find("713_Rscan_3020")
data_3080 = DatasetManager.real_data.find("713_Rscan_3080")

all_scan_data = [
  data_2000, data_2050, data_2100, data_2125, data_2150, data_2175,
  data_2200, data_2232, data_2309, data_2386, data_2396, data_2645,
  data_2900, data_2950, data_2981, data_3000, data_3020, data_3080
]

# Inclusive MC for each energy point
all_scan_incMC = all_scan_data.map { |d| DatasetManager.inclusive_mc.find(d.sample_name) }

# Exclusive signal MC across the energy scan
exMC_nnbar = DatasetManager.create_exclusive_mc_for(all_scan_data) do |config|
  config.sample_name = "nnbar_signal"
  config.events = 100000
  config.decay_card = decay_card_nnbar
  config.cross_section = :straight_line
end

# --- Algorithm ---
alg = Algorithm.new("NNbarXS")
alg.set_header(["NNbarXSAlg/NNbarXS.h"])
    .set_constant({ "ECMS" => [:double, 2.0] })

# Neutral-only event selection (no charged tracks)
sel = Selection.new
sel.select_track {
      nChrp "==0"
      nChrn "==0"
      nTot "==0"
    }
    .select_photon {
      tdc_emc_start 0
      tdc_emc_end 14
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      nGam ">=2"
    }

alg.with_decay_card(decay_card_nnbar).apply(sel)

# Most selection criteria are highly detector-specific and inexpressible
# in standard BOSS DSL; captured via note() for downstream systematic treatment
alg.note(:signal_categories, "Three statistically independent event categories:
  A: TOF + EMC response for both n and nbar (knock-off proton in TOF plastic
     scintillator + hadronic showers in EMC). N_50hits in [30,140].
  B: EMC showers for both, TOF only for nbar (|DeltaT_nbar| > 0.5 ns).
     BDT discriminator > 0.1 for background rejection.
  C: EMC showers only for both (no TOF). Most energetic shower identified as
     nbar (annihilation signature: second moment > 20 cm^2).
  Combined via inverse-variance weighting of statistically independent results.")
  .note(:antineutron_reconstruction, "Anti-neutron identified via annihilation in EMC:
  most energetic shower, |cos(theta)| < 0.7-0.8, E > 0.5 GeV.
  EMC shower second moment and N_50hits for hadronic vs EM shower discrimination.
  Neutron: EMC energy 0.06-0.70 GeV (category-dependent).
  Back-to-back topology: opening angle > 150-170 deg between n and nbar.")
  .note(:neutral_tof_algorithm, "Neutral-particle TOF time reconstruction:
  flight time measured for n and nbar candidates relative to photon expectation.
  DeltaT_n and DeltaT_nbar used for signal/background discrimination.
  Efficiency verified with e+e- -> gamma gamma and J/psi -> n nbar control samples
  (agreement with world average within 1%).")
  .note(:muc_cosmic_veto, "Muon Counter (MUC) used for cosmic ray rejection:
  last 3 layers required to have no hit response. Efficiency validated with
  J/psi -> p nbar pi- and J/psi -> pbar n pi+ control channels.
  Cosmic background studied with non-collision data samples.")
  .note(:bdt_category_b, "Boosted Decision Tree for category B using EMC and TOF
  observables; discriminator > 0.1 for signal-like events.
  Training and validation on signal MC vs multi-hadronic/di-gamma backgrounds.")
  .note(:efficiency_correction, "Data-driven efficiency corrections using control samples:
  - J/psi -> p nbar pi- and J/psi -> pbar n pi+ for n/nbar EMC and TOF efficiency
    correction matrix M_jk = eps_data(p,cos_theta) / eps_MC(p,cos_theta)
  - e+e- -> p pbar for E_extra cut efficiency and EMC neutral trigger efficiency.
  Assumption: hadronic shower in EMC from (anti-)neutron similar to (anti-)proton.")
  .note(:trigger_correction, "EMC neutral trigger efficiency modelled as:
  Trg(E) = 0.5 + 0.5 Erf((E - a)/b), parameters a = 0.758 +/- 0.005,
  b = 0.334 +/- 0.009 from e+e- -> p pbar data. Reweighted with corrected
  total energy deposition spectrum from signal process.")
  .note(:isr_vp_correction, "Radiative correction factor (1+delta) with ISR and
  vacuum polarization from ConExc generator (NLO). VP by Jegerlehner.
  Systematic from form factor parametrization variation and sampling within
  input model uncertainty band.")
  .note(:signal_extraction, "Unbinned maximum likelihood fit to:
  DeltaT_n distribution (category A) or opening angle between n and nbar
  (categories B, C). Global fit across 18 energy points minimizing NLL.
  Signal PDF from ConExc MC; background PDFs from multi-hadronic (lund),
  di-gamma (babayaga NNLO), beam-associated (non-collision data), cosmic rays.
  MINOS error analysis for parameter uncertainty.")
  .note(:oscillation_analysis, "Effective form factor |G| shows oscillatory
  behavior around dipole law G_D. Periodic structure F_osc = A exp(-B p)
  cos(C p + D) fitted simultaneously with proton data (common momentum
  frequency C = 6.5 +/- 0.1 GeV^-1). Phase difference |D_p - D_n| = 123 +/- 12 deg.
  Possible explanations: final-state re-scattering interference or resonant structure.")
  .note(:born_cross_section, "sigma_B = N_s / (L_int * eps * (1+delta)).
  sigma_B(n nbar) compared with sigma_B(p pbar) from same data:
  ratio R_np = sigma_B(nnbar)/sigma_B(ppbar) in [0.25, 1.0],
  contradicting FENICE result (R_np = 1.69 +/- 0.49 > 1).
  Photon-neutron interaction weaker than photon-proton, as theoretically expected.")
  .note(:category_consistency, "Results from categories A, B, C consistent within
  1 sigma at all energies. Inverse-variance weighted combination using
  generalized least squares with correlation matrix from shared systematics
  (luminosity, ISR, angular model).")

alg.execute_on(all_scan_data + all_scan_incMC + exMC_nnbar)