# Paper: arXiv:2401.09468v2
# Measurement of Born cross section of e+e- → Σ+ anti-Σ- at √s between 3.510 and 4.951 GeV
# 41 energy points, 24.1 fb⁻¹ total integrated luminosity
# Process: e+e- → Σ+ anti-Σ-, Σ+ → p π0, anti-Σ- → anti-p π0, π0 → γγ
# Ordinary analysis with 6C kinematic fit (4C + 2 π0 mass constraints)
# Signal MC: KKMC + PHSP model, 100k events per energy point

### Dataset preparation ###

DatasetManager.load_real_data('assets/BES3_dataset.md')
DatasetManager.load_inclusive_mc('assets/BES3_incMC.md')

# Representative datasets from the 41 energy points (3.510–4.951 GeV)
# Data taken with BOSS 703, 704, 705, 706, 707, 712 across the R-scan
# Full 41 energy points (approximate list based on BESIII R-scan configuration):
#   3510, 3582, 3650, 3670, 3773, 3808, 3867, 3871, 3896, 4009, 4085,
#   4128, 4157, 4178, 4189, 4199, 4209, 4219, 4226, 4236, 4244, 4258,
#   4267, 4278, 4288, 4308, 4315, 4340, 4358, 4377, 4396, 4416, 4436,
#   4467, 4527, 4575, 4600, 4612, 4628, 4641, 4661
# Representative subset for DSL expression:
data_3773 = DatasetManager.real_data.find("712_3773")   # 3.773 GeV (representative)
data_4180 = DatasetManager.real_data.find("703_4180")   # 4.178 GeV
data_4260 = DatasetManager.real_data.find("703_4260")   # 4.260 GeV

scan_data = [data_3773, data_4180, data_4260]

incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")
incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")

# Decay card: e+e- → Σ+ anti-Σ- via continuum (vpho)
# Σ+ → p π0, anti-Σ- → anti-p π0, π0 → γγ
# EvtGen names: Sigma+ (Σ+), anti-Sigma- (anti-Σ-), p+, anti-p-, pi0, gamma
# PHSP model for signal generation; KKMC provides ISR simulation
decay_card = <<~DECAYCARD
    Decay vpho
    1.0000  Sigma+  anti-Sigma-   PHSP;
    Enddecay

    Decay Sigma+
    1.0000  p+  pi0               PHSP;
    Enddecay

    Decay anti-Sigma-
    1.0000  anti-p-  pi0          PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma  gamma          PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC at 3.773 GeV (representative)
# 100k events per energy point; Born cross section from ConExc internal table
sig_mc = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_SigmaSigma_3773"
  config.related_dataset = data_3773
  config.events          = 100_000
  config.decay_card      = decay_card
  config.cross_section   = :default
end

### Event selection ###

algorithm = Algorithm.new("SigmaSigmaScan")
algorithm.set_header(["SigmaSigmaScanAlg/SigmaSigmaScan.h"])
  .set_constant({ "ECMS" => [:double, 3.773] })

selection = Selection.new
  .select_track do
    cos_theta 0.93       # |cosθ| < 0.93
    Vz        10.0       # |Vz| < 10 cm
    Vr        1.0        # |Vr| < 1 cm
    nChrp     ">=1"      # at least 1 positive track (proton from Σ+)
    nChrn     ">=1"      # at least 1 negative track (anti-proton from Σ-)
    nNet      "==0"      # net charge = 0
  end
  .select_photon do
    tdc_emc_start     0         # EMC timing > 0 ns
    tdc_emc_end       14        # EMC timing < 700 ns (14 × 50 ns)
    energyThreshold_b 0.025     # barrel: E > 25 MeV
    energyThreshold_e 0.050     # endcap: E > 50 MeV
    angle_to_track    10.0      # > 10° from any charged track
    nGam              ">=4"     # at least 4 photons (2 π0 → 4γ)
  end
  .pid(method: :probability) do
    prob_cut 0.001
    # Proton/anti-proton identification: probability-based PID as DSL placeholder.
    # The actual analysis uses momentum-based PID: p > 0.5 GeV/c.
    # See algorithm.note(:proton_momentum_pid) for details.
    identify :proton, against: [:kaon, :pion]
    nprp ">=1"                  # at least 1 proton from Σ+ decay
    nprm ">=1"                  # at least 1 anti-proton from Σ- decay
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    # 1C Kalman fit per π0 candidate: constrain γγ invariant mass to π0 nominal mass
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=2"                  # at least 2 π0 candidates from 4 photons
  end
  # 6C kinematic fit: 4C energy-momentum + 2C from pi0 mass constraints (already applied via kalman_kinematic_fit)
  .kinematic_fit([:prp, :prm, :pi0, :pi0]) do
    nominal
    constrain_four_momentum     # 4C: total energy-momentum conservation
    chi2_cut 200                # loose BOSS cut; tighter selection applied in ROOT
  end

# ============================================================================
# Algorithm Notes
# ============================================================================

algorithm.note(:multi_energy,
  "Multi-energy measurement at 41 CM energy points from 3.510 to 4.951 GeV, " \
  "covering a total integrated luminosity of 24.1 fb⁻¹. ECMS varies per run; " \
  "3.773 GeV is used as a representative value in the algorithm constant. " \
  "Per-run beam energy is injected at runtime via MeasuredEcmsSvc. " \
  "Representative BOSS datasets shown: 712_3773 (3.773 GeV), 703_4180 (4.178 GeV), " \
  "703_4260 (4.260 GeV). The complete 41-point dataset list should be loaded " \
  "from a YAML configuration file; this representative subset suffices for DSL " \
  "expression of the selection logic.")

algorithm.note(:born_cross_section_formula,
  "Born cross section: σ_B = N_sig / [L_int · ε · B · (1+δ_ISR) · 1/|1-Π|²], " \
  "where N_sig is the signal yield extracted from fits to the M(pπ⁰) and " \
  "M(p̄π⁰) distributions, ε is the detection efficiency, " \
  "B = B(Σ⁺→pπ⁰) · B(Σ̄⁻→p̄π⁰) · B(π⁰→γγ)² is the product of branching " \
  "fractions, (1+δ_ISR) is the ISR correction factor, and 1/|1-Π|² is the " \
  "vacuum polarization correction. The dressed cross section is obtained as " \
  "σ_D = σ_B · |1-Π|². All cross section calculations are performed at the " \
  "ROOT analysis level.")

algorithm.note(:isr_vp_corrections,
  "ISR correction factor (1+δ_ISR) obtained iteratively using the KKMC event " \
  "generator with a PHSP signal model for e⁺e⁻ → Σ⁺Σ̄⁻. Vacuum polarization " \
  "factor 1/|1-Π|² obtained from ConExc generator internal tables (ConExc " \
  "mode 44 for charged Σ pair production, valid for 2.308–5.000 GeV). Both " \
  "correction factors are applied at the ROOT level during cross section " \
  "extraction. The ISR correction accounts for the radiative tail and " \
  "energy-dependent efficiency variations.")

algorithm.note(:sideband_background,
  "Sideband background subtraction method: signal and sideband regions are " \
  "defined in the M(pπ⁰) and M(p̄π⁰) invariant mass distributions. The Σ⁺ " \
  "signal region is |M(pπ⁰) − m_Σ⁺| < 2.5σ (∼2–3 MeV/c² around " \
  "m_Σ⁺ = 1.18937 GeV/c²) and the Σ̄⁻ signal region is analogously defined " \
  "around m_Σ̄⁻. Sideband regions at 5σ–7.5σ on both sides of the Σ peak " \
  "are used to estimate the combinatorial background through interpolation. " \
  "Background subtraction is performed at the ROOT analysis level.")

algorithm.note(:iterative_efficiency,
  "Iterative efficiency correction method: the detection efficiency depends " \
  "on the Born cross section line shape, which is unknown a priori. An " \
  "iterative procedure is employed: (1) initial efficiency obtained from " \
  "PHSP signal MC assuming a flat cross section; (2) Born cross section " \
  "extracted using this efficiency; (3) signal MC reweighted according to " \
  "the extracted cross section line shape; (4) updated efficiency computed " \
  "from reweighted MC; (5) steps 2–4 repeated until convergence (typically " \
  "2–3 iterations). The iteration is performed at the ROOT analysis level.")

algorithm.note(:proton_momentum_pid,
  "Proton/anti-proton identification by momentum threshold: charged tracks " \
  "with momentum p > 0.5 GeV/c in the laboratory frame are assigned as proton " \
  "candidates (positive charge) or anti-proton candidates (negative charge). " \
  "This momentum-based criterion replaces the standard probability-based PID " \
  "and is motivated by the fact that protons from Σ decay have momenta " \
  "well above 0.5 GeV/c at BESIII scan energies, while pions and kaons from " \
  "background processes typically have lower momenta. The DSL probability " \
  "PID block is a placeholder; the actual momentum threshold must be " \
  "implemented in the C++ analysis algorithm code.")

algorithm.note(:sigma_mass_window,
  "Σ mass window definition for signal and sideband regions. Signal region: " \
  "|M(pπ⁰) − m_Σ⁺| < 2.5σ and |M(p̄π⁰) − m_Σ̄⁻| < 2.5σ, where σ is the " \
  "mass resolution obtained from signal MC (∼1 MeV/c², dominated by the " \
  "π⁰ mass resolution). Sideband regions: 5σ < |M − m_Σ| < 7.5σ on both " \
  "sides of the Σ peak. Events in the signal region are used for yield " \
  "extraction; events in the sideband regions are used for background " \
  "estimation via linear interpolation. Applied at ROOT level.")

algorithm.note(:bhabha_suppression,
  "Bhabha background suppression: requirement E/p < 0.8 on all charged " \
  "tracks, where E is the energy deposited in the EMC and p is the track " \
  "momentum measured in the MDC. This cut rejects e⁺e⁻ → e⁺e⁻ (Bhabha " \
  "scattering) events where electrons can be misidentified as protons or " \
  "anti-protons due to shower fluctuations in the EMC. Implemented as a " \
  "track-level cut in the C++ analysis code; not expressible at the DSL " \
  "PID level.")

algorithm.note(:best_sigma_pair,
  "Best Σ⁺Σ̄⁻ pair selection: when multiple proton/anti-proton candidates " \
  "and π⁰ combinations survive the kinematic fit, a two-step selection is " \
  "applied. (1) The combination with the minimum χ² from the 6C kinematic " \
  "fit is retained. (2) If multiple Σ⁺Σ̄⁻ pairs remain after step (1), the " \
  "pair with the minimum mass difference |M(pπ⁰) − M(p̄π⁰)| is selected. " \
  "This two-step procedure minimizes both the kinematic fit quality and the " \
  "Σ⁺/Σ̄⁻ mass asymmetry, improving the signal purity. Implemented in the " \
  "C++ analysis algorithm code.")

algorithm.note(:signal_mc,
  "Signal MC generated with the KKMC event generator using a PHSP (phase " \
  "space) model for e⁺e⁻ → Σ⁺Σ̄⁻ at each of the 41 energy points. 100,000 " \
  "events per energy point. The PHSP model assumes an isotropic angular " \
  "distribution for Σ⁺Σ̄⁻ production, which is a reasonable approximation " \
  "for the near-threshold energy region. ConExc internal Born cross section " \
  "tables (mode 44, charged Σ pair production) provide the normalization " \
  "for the ISR calculation. The exclusive MC block shown is a representative " \
  "sample at 3.773 GeV; in practice one sample per energy point is needed.")

algorithm.note(:systematics,
  "Systematic uncertainties evaluated per energy point and combined in " \
  "quadrature. Sources include: tracking efficiency (1.0% per track), " \
  "proton/anti-proton PID efficiency (momentum-dependent, evaluated using " \
  "control samples), photon detection efficiency (1.0% per photon, ∼4% " \
  "for 4 photons), π⁰ reconstruction efficiency, 6C kinematic fit " \
  "efficiency, signal/sideband mass window definition, ISR correction " \
  "factor uncertainty, vacuum polarization correction uncertainty, " \
  "integrated luminosity (1.0%), branching fraction uncertainties from " \
  "PDG, MC statistics, and fit range dependence. All systematic studies " \
  "are performed at the ROOT analysis level.")

algorithm.with_decay_card(decay_card).apply(selection)
algorithm.execute_on(scan_data + [incMC_3773, incMC_4180, sig_mc])