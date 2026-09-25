# Single Inclusive π± and K± Production in e+e- Annihilation at √s 2.000-3.671 GeV
# Paper: 2502.16084v2
# This is a continuum cross-section measurement across 8 energy points with generic
# hadronic event selection + PID. No specific decay channel is reconstructed; the
# analysis measures differential cross sections dσ/(σ_had dp) for π± and K±.
# Much of the analysis logic (hadronic event selection with Bhabha/γγ rejection,
# prong counting, PID efficiency matrix inversion, cross-section extraction,
# ISR/beam-energy-spread corrections) is done in ROOT and inexpressible in DSL.

DatasetManager.load_real_data('assets/BES3_dataset.md')
DatasetManager.load_inclusive_mc('assets/BES3_incMC.md')

# Five energy points with exact matches in the 713 Rscan dataset table:
data_2000 = DatasetManager.real_data.find("713_Rscan_2000")
data_2200 = DatasetManager.real_data.find("713_Rscan_2200")
data_2396 = DatasetManager.real_data.find("713_Rscan_2396")
data_2644 = DatasetManager.real_data.find("713_Rscan_2644")
data_2900 = DatasetManager.real_data.find("713_Rscan_2900")

incMC_2000 = DatasetManager.inclusive_mc.find("713_Rscan_2000")
incMC_2200 = DatasetManager.inclusive_mc.find("713_Rscan_2200")
incMC_2396 = DatasetManager.inclusive_mc.find("713_Rscan_2396")
incMC_2644 = DatasetManager.inclusive_mc.find("713_Rscan_2644")
incMC_2900 = DatasetManager.inclusive_mc.find("713_Rscan_2900")

scan_data  = [data_2000, data_2200, data_2396, data_2644, data_2900]
scan_incMC = [incMC_2000, incMC_2200, incMC_2396, incMC_2644, incMC_2900]

# ConExc R-value mode for continuum MC generation
decay_card_continuum = <<~DECAYCARD
  Decay vpho
  1 ConExc 74110;
  Enddecay
  End
DECAYCARD

exMC_continuum = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name   = "sig_continuum_hadronic"
  config.events        = 1_000_000
  config.decay_card    = decay_card_continuum
  config.cross_section = :default
end

algorithm = Algorithm.new("InclusiveHadronAnalysis")
algorithm
  .set_header(["InclusiveHadronAnalysisAlg/InclusiveHadronAnalysis.h"])

event_selection = Selection.new
  .select_track do
    cos_theta   0.93
    Vz          10.0
    Vr          1.0
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    identify :kaon, against: [:pion, :proton]
  end

algorithm
  .note(:analysis_type,
    "This is an inclusive cross-section measurement, not a specific decay channel reconstruction. " \
    "There is no kinematic fit — π± and K± are identified from all hadronic events via PID. " \
    "The DSL scope covers only track/shower selection and PID; the rest is handled in ROOT.")
  .note(:hadronic_event_selection,
    "Inclusive hadronic events selected per Ref.[13]: (1) Bhabha (e+e-→e+e-) rejected by dedicated EMC " \
    "shower requirements on cluster timing, energy, and topology; (2) γγ events rejected by EMC shower " \
    "balance and acollinearity cuts; (3) events with 0 or 1 prong removed to suppress QED backgrounds; " \
    "(4) 2/3-prong events pass additional QED rejection (acoplanarity, momentum balance, shower energy); " \
    "(5) >=4 prong events treated as hadronic directly. Residual QED background from MC (BABAYAGA, " \
    "DIAG36, EKHARA, GALUGA). Beam-associated background via Vz_evt sideband method.")
  .note(:pid_method,
    "PID via dE/dx + TOF combined probability P(h) for h=π,K,p. Track identified as species with " \
    "highest probability. PID efficiency matrix (ε_g→h for g,h=π,K,p) determined from control samples: " \
    "π± via J/ψ→π+π-π0, K± via J/ψ→KS K± π∓, p via J/ψ/ψ(2S)→p pbar π+π-. " \
    "Matrix inversion applied per (p, cosθ) bin to correct for misidentification. " \
    "Electron/muon contamination suppressed by hadronic event selection; opposite-charge mis-ID negligible.")
  .note(:cross_section_formula,
    "dσ/(σ_had dp) = (N_obs_πK / N_had_obs) × (1/Δp) × f_πK, where f_πK is the correction factor " \
    "accounting for detection efficiency and ISR effects. N_had_obs = N_had_tot - N_bkg, with N_bkg " \
    "estimated from QED MC and beam-associated background sideband.")
  .note(:energy_points,
    "Table I lists 8 energy points: 2.0000, 2.2000, 2.3960, 2.6444, 2.9000, 3.0500, 3.5000, 3.6710 GeV. " \
    "Five points (2.000, 2.200, 2.396, 2.644, 2.900) match 713 Rscan entries. " \
    "Three points (3.050, 3.500, 3.671) have no standard 713 Rscan sample name; they require " \
    "dedicated dataset lookups from the corresponding BOSS run periods.")
  .note(:mc_generation,
    "ConExc R-value mode 74110 (light hadrons, 0.36-5.00 GeV) used for continuum MC generation. " \
    "Multi-energy scan: ECMS injected per-point by execute_on from each dataset's CMS energy.")
  .with_decay_card(decay_card_continuum)
  .apply(event_selection)
  .execute_on(scan_data + scan_incMC + exMC_continuum)