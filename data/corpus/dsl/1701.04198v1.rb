# =============================================================================
# 1701.04198v1
# Cross section measurements of e+e- -> p pbar pi0
# at center-of-mass energies between 4.008 and 4.600 GeV (BESIII)
#
# BOSS-side spec: dataset preparation + event selection up to the final 4C
# kinematic fit. The Born cross-section extraction, the PWA and the Y(4260)
# upper-limit fit are ROOT-level and are not part of this spec.
# =============================================================================

### ---------------------------------------------------------------------------
### Datasets — 13 energy points between sqrt(s) = 4.008 and 4.600 GeV
### Sample names chosen so that the tabulated luminosities match Table 1 of
### the paper (the narrow 4.189/4.208/4.217/4.242/4.308/4.387 points come from
### the 4260/4360 scans).
### ---------------------------------------------------------------------------
data_4008 = DatasetManager.real_data.find("703_4009")        # 482.0 pb^-1
data_4085 = DatasetManager.real_data.find("703_4090")        #  52.9 pb^-1
data_4189 = DatasetManager.real_data.find("703_4190scan")    #  43.3 pb^-1
data_4208 = DatasetManager.real_data.find("703_4210scan")    #  55.0 pb^-1
data_4217 = DatasetManager.real_data.find("703_4220scan")    #  54.6 pb^-1
data_4226 = DatasetManager.real_data.find("703_4230")        # 1056.4 pb^-1
data_4242 = DatasetManager.real_data.find("703_4245")        #  55.9 pb^-1
data_4258 = DatasetManager.real_data.find("703_4260")        # 828.4 pb^-1
data_4308 = DatasetManager.real_data.find("703_4310")        #  45.1 pb^-1
data_4358 = DatasetManager.real_data.find("703_4360")        # 543.9 pb^-1
data_4387 = DatasetManager.real_data.find("703_4390")        #  55.6 pb^-1
data_4416 = DatasetManager.real_data.find("703_4420")        # 1043.9 pb^-1
data_4600 = DatasetManager.real_data.find("703_4600")        # 586.9 pb^-1

data_points = [data_4008, data_4085, data_4189, data_4208, data_4217, data_4226,
               data_4242, data_4258, data_4308, data_4358, data_4387, data_4416,
               data_4600]

# Inclusive MC is only available for a subset of the scan points.
incMC_points = [DatasetManager.inclusive_mc.find("703_4009"),
                DatasetManager.inclusive_mc.find("703_4190scan"),
                DatasetManager.inclusive_mc.find("703_4210scan"),
                DatasetManager.inclusive_mc.find("703_4220scan"),
                DatasetManager.inclusive_mc.find("703_4230"),
                DatasetManager.inclusive_mc.find("703_4260"),
                DatasetManager.inclusive_mc.find("703_4360"),
                DatasetManager.inclusive_mc.find("703_4420"),
                DatasetManager.inclusive_mc.find("703_4600")]

### ---------------------------------------------------------------------------
### Decay cards
### ---------------------------------------------------------------------------
# Signal: e+e- -> p pbar pi0. This is a continuum Born-cross-section
# measurement (ISR and vacuum-polarisation corrections enter the published
# result), so a ConExc card is used: mode 50 = p pbar pi0.
# 'Particle vpho' is intentionally OMITTED so the DSL injects the correct
# per-point sqrt(s) into every scan job.
decay_card_signal = <<~DECAYCARD
  Decay vpho
  1 ConExc 50;
  Enddecay

  Decay vhdr
  1 p+ anti-p- pi0 PHSP;
  Enddecay

  Decay pi0
  1 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Background: inclusive Y(4260) decays at sqrt(s) = 4.26 GeV (equiv. to
# 825.6 pb^-1), generated with the KKMC + psi(4260) convention.
decay_card_y4260_inclusive = <<~DECAYCARD
  Decay psi(4260)
  1.0000 generic generic generic GENERIC;
  Enddecay

  End
DECAYCARD

### ---------------------------------------------------------------------------
### Exclusive MC — signal p pbar pi0 at all 13 energy points (200k events each)
### ---------------------------------------------------------------------------
exMC_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "ppbarpi0_signal_exclusive_mc"   # auto-suffixed per energy point
  config.events        = 200_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default                          # required field; inert for ConExc
end

# Inclusive Y(4260) background MC at sqrt(s) = 4.26 GeV
exMC_y4260_inclusive = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "y4260_inclusive_mc"
  config.related_dataset = data_4258
  config.events          = 500_000
  config.decay_card      = decay_card_y4260_inclusive
  config.cross_section   = :default
end

### ---------------------------------------------------------------------------
### Event selection (BOSS)
### ---------------------------------------------------------------------------
alg_name = "PpbarPi0"
algorithm = Algorithm.new(alg_name)
algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
# Multi-energy cross-section scan: ECMS is injected per job at run time, so it
# is deliberately NOT set via set_constant here.

event_selection = Selection.new
# --- Charged tracks: exactly two tracks with opposite charge -----------------
event_selection
  .select_track do
    cos_theta 0.93     # |cos(theta)| < 0.93
    Vz        10.0     # |Vz| < 10 cm along the beam direction
    Vr        1.0      # Vr < 1 cm in the plane perpendicular to the beam
    nChrp     "==1"
    nChrn     "==1"
    nNet      "==0"
  end
  # --- Photon candidates: >= 2 photons ----------------------------------------
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025    # 25 MeV in the barrel  (|cos(theta)| < 0.8)
    energyThreshold_e 0.050    # 50 MeV in the endcap (0.86 < |cos(theta)| < 0.92)
    nGam              ">=2"
  end
  # --- PID: one track identified as proton, the other as anti-proton ---------
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]   # p+ and anti-p- in one line
    nprp     "==1"
    nprm     "==1"
  end
  # --- Photon isolation: > 10 deg from the proton, > 30 deg from the
  #     anti-proton (suppresses photons from anti-proton annihilation) --------
  .select_isolated_photon do
    angle_to_prp_track 10.0
    angle_to_prm_track 30.0
    nGam               ">=2"
  end
  # --- 4C kinematic fit on p pbar gamma gamma --------------------------------
  .kinematic_fit([:prp, :prm, :gamma, :gamma]) do
    nominal
    constrain_four_momentum   # total 4-momentum constrained to the e+e- CMS
    chi2_cut 200              # loose BOSS cut; the paper's chi2_4C < 30 is applied in ROOT
  end

algorithm
  .note(:kinematic_fit_chi2, "the paper keeps the combination with the smallest chi2_4C among all p pbar gamma gamma pairings and requires chi2_4C < 30; BOSS applies only the loose default chi2_cut 200 and the tight cut is applied in ROOT")
  .note(:pi0_mass_window, "pi0 candidates are selected with |M(gamma gamma) - m_pi0| < 15 MeV/c2; this window uses post-kinematic-fit quantities and is therefore applied in the ROOT analysis")
  .note(:pwa_efficiency, "detection efficiency is derived from MC generated according to the partial wave analysis result (N(1440)/N(1520), rho(2150), rho3(1990) and 1-- PHSP components); the PWA fit itself is not expressible in the BOSS DSL")
  .note(:pi0_sideband_background, "non-pi0 background is estimated from the pi0 sidebands 0.07 < M(gamma gamma) < 0.10 GeV/c2 and 0.17 < M(gamma gamma) < 0.20 GeV/c2; the contamination is 0.3% at sqrt(s) = 4.258 GeV and is neglected")
  .note(:conexc_mode_range, "ConExc mode 50 (p pbar pi0) is tabulated for sqrt(s) = 4.009-4.20 GeV; the higher scan points of this analysis lie beyond that range and the ISR line shape should be cross-checked against the measured cross section used iteratively in the paper")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

### ---------------------------------------------------------------------------
### Execution
### ---------------------------------------------------------------------------
root_files = algorithm.execute_on(data_points + incMC_points + [exMC_signal, exMC_y4260_inclusive].flatten)
