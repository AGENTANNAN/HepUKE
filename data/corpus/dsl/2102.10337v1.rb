# =============================================================================
# BESIII paper arXiv:2102.10337v1
# Measurement of proton electromagnetic form factors in the time-like region
# using initial state radiation (ISR) at BESIII.
#
# Process: e+e- -> p pbar gamma_ISR using 7.5 fb-1 at seven c.m. energy points
# from 3.773 to 4.600 GeV.
#
# The ISR photon is from initial-state radiation and identified as the
# highest-energy photon in the event. The signal is modelled with the KKMC
# generator (psi(4260) as the top mother), not ConExc.
# =============================================================================

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

### Dataset preparation — seven energy points covering 3.773-4.600 GeV ###
data_points = [
  DatasetManager.real_data.find("712_3773"),   # sqrt(s) = 3.773 GeV, ~2.9 fb-1
  DatasetManager.real_data.find("703_4180"),   # sqrt(s) = 4.180 GeV
  DatasetManager.real_data.find("703_4230"),   # sqrt(s) = 4.230 GeV
  DatasetManager.real_data.find("703_4260"),   # sqrt(s) = 4.260 GeV
  DatasetManager.real_data.find("703_4360"),   # sqrt(s) = 4.360 GeV
  DatasetManager.real_data.find("703_4420"),   # sqrt(s) = 4.420 GeV
  DatasetManager.real_data.find("703_4600"),   # sqrt(s) = 4.600 GeV
]

incMC_points = data_points.map { |dp| DatasetManager.inclusive_mc.find(dp.sample_name) }

# Signal decay card: e+e- -> p pbar gamma (KKMC + psi(4260) top mother,
# gamma is the ISR photon).
decay_card_signal = <<~DECAYCARD
  Decay psi(4260)
  1 p+ anti-p- gamma PHSP;
  Enddecay
  End
DECAYCARD

# Exclusive signal MC: one sample per c.m. energy point, so the detection
# efficiency is evaluated at the same sqrt(s) as the data.
exMC_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "sig_ppbar_gamma_isr"
  config.events        = 200_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

# =============================================================================
# Event selection
# =============================================================================
# Multi-energy measurement: ECMS is not set as a constant here; each job
# takes its c.m. energy from the dataset at run time.
alg_name = "PpbarGammaISR"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])

sel = Selection.new

# --- Charged track selection ---------------------------------------------------
sel.select_track do
      cos_theta 0.93     # |cos(theta)| < 0.93
      Vz        10.0     # |Vz| < 10 cm along the beam direction
      Vr        1.0      # Vr < 1 cm in the plane perpendicular to the beam
      nChrp     "==1"    # exactly one positively charged track (proton candidate)
      nChrn     "==1"    # exactly one negatively charged track (antiproton candidate)
      nNet      "==0"    # net charge zero
    end

# --- Photon selection ----------------------------------------------------------
sel.select_photon do
      tdc_emc_start     0
      tdc_emc_end       14
      energyThreshold_b 0.025    # E > 25 MeV in barrel EMC (|cos(theta)| < 0.80)
      energyThreshold_e 0.050    # E > 50 MeV in endcap EMC (0.86 < |cos(theta)| < 0.92)
      angle_to_track    10.0     # > 10 deg from the nearest charged track
      nGam              ">=1"    # at least one photon (the ISR photon)
    end

# --- Particle identification ---------------------------------------------------
# Combined dE/dx and TOF confidence levels: the proton hypothesis must have
# higher probability than the kaon or pion hypotheses for both tracks.
sel.pid(method: :probability) do
      prob_cut 0.001
      identify :proton, against: [:kaon, :pion]
      nprp "==1"
      nprm "==1"
    end

# --- E/p requirement for proton candidates (Bhabha suppression) ----------------
# E/p < 0.5 is required for both proton and antiproton, where E is the energy
# deposited in the EMC and p is the momentum measured in the MDC. Candidates
# failing this cut (E/p >= 0.5) are removed.
sel.for_each(:prp) do
      define(:e_over_p) { eraw / p }
      where { e_over_p >= 0.5 }
      remove
    end
   .for_each(:prm) do
      define(:e_over_p) { eraw / p }
      where { e_over_p >= 0.5 }
      remove
    end

# --- ISR photon isolation -----------------------------------------------------
# Photons from antiproton/proton interactions in the detector material are
# suppressed by requiring the photon to be well separated from both proton tracks.
sel.select_isolated_photon do
      angle_to_prp_track 20.0     # > 20 deg from the proton track
      angle_to_prm_track 20.0     # > 20 deg from the antiproton track
      nGam               ">=1"
    end

# --- High-energy ISR photon requirement ----------------------------------------
# The ISR photon is the highest-energy photon in the event and must have
# E > 0.4 GeV: photons below this threshold are removed, and the
# highest-energy remaining photon is tagged as the ISR candidate.
sel.for_each(:gamma) do
      where { energy < 0.4 }
      remove
    end
   .for_each(:gamma) do
      define(:egamma) { energy }
      best { maximize { "egamma" } }
      store(:isr_gamma_index)
    end

# --- 4C kinematic fit ----------------------------------------------------------
# Four-momentum conservation constrains the p pbar gamma system to the initial
# e+e- four-momentum (four constraints).
sel.kinematic_fit([:prp, :prm, :gamma]) do
      nominal
      constrain_four_momentum
      chi2_cut 50
    end

# --- Notes — inexpressible BOSS-side procedures -------------------------------
alg.note(:isr_technique,
         "The ISR technique is used: the ISR photon is radiated from the initial " \
         "e+e- state, and its energy determines the effective c.m. energy of the " \
         "p pbar system as s' = s - 2*sqrt(s)*E_gamma. This allows a continuous " \
         "scan of the p pbar invariant mass from threshold up to the full c.m. energy.")
   .note(:generator,
         "Signal is modelled with the KKMC generator using psi(4260) as the top " \
         "mother, not ConExc. The decay p+ anti-p- gamma is generated with a " \
         "phase-space (PHSP) model.")
   .note(:ep_cut,
         "E/p < 0.5 for protons and antiprotons suppresses Bhabha (e+e-) background " \
         "where electrons deposit more energy in the EMC per unit momentum than protons.")
   .note(:isr_photon_tag,
         "The ISR photon is identified as the highest-energy photon with " \
         "E > 0.4 GeV in the event. Photons below 0.4 GeV are removed before the " \
         "kinematic fit; the highest-energy remaining photon is stored.")
   .note(:cos_theta_p,
         "The proton polar angle cos(theta_p) in the p pbar rest frame is used to " \
         "extract the form factor ratio |G_E|/|G_M|. This analysis is performed at " \
         "the ROOT level on the kinematic-fit-corrected four-momenta.")
   .note(:form_factor_extraction,
         "The proton electromagnetic form factors |G_E| and |G_M| and their ratio " \
         "R = |G_E|/|G_M| are extracted from the p pbar invariant mass distribution " \
         "and the proton angular distribution cos(theta_p). The efficiency correction " \
         "depends on the form factor ratio and is iterated until convergence. " \
         "Performed at the ROOT level.")
   .note(:cross_section,
         "The Born cross section sigma(e+e- -> p pbar) is measured as a function of " \
         "the p pbar invariant mass M(p pbar) from threshold up to ~4.6 GeV/c^2, " \
         "combining data from all seven energy points via the ISR technique. The " \
         "radiative correction factor (1+delta) is applied to convert the observed " \
         "cross section to the Born-level cross section.")
   .note(:energy_points,
         "Seven c.m. energy points: 3.773 GeV (712_3773, ~2.9 fb-1), " \
         "4.180 GeV (703_4180), 4.230 GeV (703_4230), 4.260 GeV (703_4260), " \
         "4.360 GeV (703_4360), 4.420 GeV (703_4420), 4.600 GeV (703_4600). " \
         "Total integrated luminosity ~7.5 fb-1.")
   .note(:background_study,
         "Main backgrounds: e+e- -> pi+ pi- (misidentified as p pbar), " \
         "e+e- -> e+ e- (Bhabha, suppressed by E/p), e+e- -> mu+ mu-, and " \
         "e+e- -> p pbar pi0 (with a missing pi0). Background yields are estimated " \
         "from inclusive MC and sideband studies at the ROOT level.")
   .note(:systematic_uncertainties,
         "Systematic uncertainties include tracking efficiency, PID efficiency, " \
         "photon detection efficiency, ISR photon selection, kinematic fit, " \
         "luminosity, radiative corrections, generator model dependence, and the " \
         "form-factor-dependent efficiency correction. Estimated at the ROOT level.")
   .note(:bhabha_suppression,
         "In addition to the E/p < 0.5 cut, Bhabha events are further suppressed " \
         "by requiring cos(theta_p) < 0.8 at c.m. energies above threshold where " \
         "protons are forward-peaked. Applied at the ROOT level.")

# --- Execution ----------------------------------------------------------------
alg.with_decay_card(decay_card_signal).apply(sel)
alg.execute_on(data_points + incMC_points + exMC_signal)