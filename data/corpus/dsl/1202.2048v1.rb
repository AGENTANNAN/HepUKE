# =============================================================================
# BESIII: Precision measurement of the branching fractions of
#         J/psi -> pi+ pi- pi0  and  psi' -> pi+ pi- pi0
# arXiv:1202.2048v1
# J/psi sample: 2.25 x 10^8 ; psi' sample: 1.06 x 10^8 (2009 data)
# =============================================================================

### Dataset preparation ###
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# ---- Decay cards -------------------------------------------------------------
# J/psi -> pi+ pi- pi0 : the decay is dominated by the intermediate rho(770) pi state,
# which gives a good description of the data.
decay_card_jpsi = <<~DECAYCARD
  Decay J/psi
  1.000 rho0 pi0   VSS;
  Enddecay

  Decay rho0
  1.000 pi+ pi-    VSS;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# psi' -> pi+ pi- pi0 : a mixture of rho(770) pi and P-wave phase space is used,
# modelled here with the phase-space generator.
decay_card_psip = <<~DECAYCARD
  Decay psi(2S)
  1.000 pi+ pi- pi0 PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# ---- Exclusive MC samples ----------------------------------------------------
exMC_jpsi = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "jpsi_pip_pim_pi0"
  c.related_dataset = jpsi_data
  c.events          = 1000000
  c.decay_card      = decay_card_jpsi
  c.cross_section   = :default
end

exMC_psip = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "psip_pip_pim_pi0"
  c.related_dataset = psip_data
  c.events          = 1000000
  c.decay_card      = decay_card_psip
  c.cross_section   = :default
end

# =============================================================================
# ALGORITHM 1: J/psi -> pi+ pi- pi0   (2.25 x 10^8 J/psi events)
# =============================================================================
alg_jpsi = Algorithm.new("JpsiToPiPiPi0")
alg_jpsi.set_header(["JpsiToPiPiPi0Alg/JpsiToPiPiPi0.h"])
        .set_constant({"ECMS" => [:double, 3.097]})   # J/psi center-of-mass energy (GeV)

sel_jpsi = Selection.new
sel_jpsi.select_track {
           cos_theta  0.93   # |cos(theta)| < 0.93
           Vz         10.0   # within +-10 cm of the IP along the beam direction
           Vr         1.0    # within 1 cm of the beam line in the transverse plane
           nChrp      "==1"  # exactly one positively charged track (pi+)
           nChrn      "==1"  # exactly one negatively charged track (pi-)
           nNet       "==0"
         }
         .select_photon {
           # EMC showers: > 25 MeV in the barrel (|cos|<0.8), > 50 MeV in the endcaps
           # (0.86<|cos|<0.92); the barrel/endcap transition region is excluded and
           # showers within 10 degrees of a charged track are rejected
           tdc_emc_start     0
           tdc_emc_end       14
           angle_to_track    10.0
           energyThreshold_b 0.025
           energyThreshold_e 0.050
           nGam              ">=2"   # at least two photon candidates for the pi0
         }
         .assign({chrgp: :pip, chrgn: :pim})   # pi+ / pi- (no hadron PID required)
         .kalman_kinematic_fit([:gamma, :gamma]) {
           # pi0 mass constraint; the photon pair with the smallest chi2 is kept as the
           # pi0 candidate if chi2 < 50, otherwise the event is rejected
           invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
           chi2_cut 50
           npi0 ">=1"
         }
         # Full-event kinematic fit with the initial J/psi four-momentum as constraint
         # (4C) together with the pi0 mass constraint (1C) -> 5C fit, chi2 < 50
         .kinematic_fit([:pip, :pim, :pi0]) {
           nominal
           constrain_four_momentum
           chi2_cut 50
         }
         # Competing hypothesis: the same event fitted with the charged tracks assumed to
         # be kaons. If the kaon hypothesis gives a smaller chi2 the event is rejected
         # (veto applied in the ROOT analysis from the stored chi2 value).
         .kinematic_fit([:kp, :km, :pi0]) {
           use_track_index_from_nominal_kmfit
           constrain_four_momentum
         }

alg_jpsi
  .note(:kaon_hypothesis_veto,
        "the same event is refitted under the K+ K- pi0 hypothesis (kinematic fit repeated with the " \
        "charged particles assumed to be kaons); if this leads to a smaller chi2 the event is rejected. " \
        "The competing chi2 is stored in the NTuple and the veto is applied in the ROOT analysis.")
  .note(:pi0_mass_window,
        "the invariant mass of the two photon candidates is required to be compatible with the pi0 mass, " \
        "0.11 < m(gamma gamma) < 0.15 GeV/c^2; this uses the kinematic-fit-corrected four-momenta and is " \
        "applied at the ROOT level.")
  .note(:track_efficiency_correction,
        "the tracking efficiency in simulation is corrected as a function of polar angle and track " \
        "momentum to match a tagged J/psi -> 3pi control sample (on average about 2% lower in data).")
  .note(:pi0_efficiency_correction,
        "the pi0 reconstruction efficiency in simulation is corrected as a function of pi0 momentum " \
        "using a control sample of two tracks plus two photons (differences of order 0.5%).")
  .note(:reweighting,
        "the signal MC events are reweighted in the Dalitz variables to the observed data distribution " \
        "to account for differences between the generated and observed distributions.")
  .note(:continuum_background,
        "background from the continuum is estimated from off-resonance data samples taken at 3.08 GeV " \
        "(282 nb^-1) and 3.650 GeV, normalised by luminosity and cross-section energy dependence.")
  .note(:resonance_background,
        "background from other resonance processes is estimated using the inclusive MC sample; the " \
        "inclusive simulation is assigned a 100% uncertainty based on checks with pi+pi-gamma and " \
        "pi+pi-pi0 gamma final states.")
  .note(:trigger_efficiency,
        "the trigger efficiency is taken as (100 - 0.2)% for hadronic events containing charged particles.")

alg_jpsi.with_decay_card(decay_card_jpsi).apply(sel_jpsi)
root_files_jpsi = alg_jpsi.execute_on([jpsi_data, jpsi_incMC, exMC_jpsi])

# =============================================================================
# ALGORITHM 2: psi' -> pi+ pi- pi0   (1.06 x 10^8 psi' events)
# =============================================================================
alg_psip = Algorithm.new("PsipToPiPiPi0")
alg_psip.set_header(["PsipToPiPiPi0Alg/PsipToPiPiPi0.h"])
        .set_constant({"ECMS" => [:double, 3.686]})   # psi(2S) center-of-mass energy (GeV)

sel_psip = Selection.new
sel_psip.select_track {
           cos_theta  0.93
           Vz         10.0
           Vr         1.0
           nChrp      "==1"
           nChrn      "==1"
           nNet       "==0"
         }
         .select_photon {
           tdc_emc_start     0
           tdc_emc_end       14
           angle_to_track    10.0
           energyThreshold_b 0.025
           energyThreshold_e 0.050
           nGam              ">=2"
         }
         # Additional suppression of radiative psi' -> gamma e+ e- / mu+ mu- and of
         # psi' -> pi+ pi- J/psi: the energy deposit associated with a track must be
         # less than 0.8 GeV (removes electrons)
         .for_each(:charged) {
           where { eraw > 0.8 }
           remove
         }
         .assign({chrgp: :pip, chrgn: :pim})
         .kalman_kinematic_fit([:gamma, :gamma]) {
           invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
           chi2_cut 50
           npi0 ">=1"
         }
         .kinematic_fit([:pip, :pim, :pi0]) {
           nominal
           constrain_four_momentum
           chi2_cut 50
         }
         .kinematic_fit([:kp, :km, :pi0]) {
           use_track_index_from_nominal_kmfit
           constrain_four_momentum
         }

alg_psip
  .note(:muon_veto,
        "the penetration depth into the muon system (MUC) is required to be less than 40 cm for each " \
        "charged track, suppressing psi' -> mu+ mu- and psi' -> pi+ pi- J/psi (J/psi -> l+ l-) backgrounds.")
  .note(:pipi_mass_veto,
        "the invariant mass of the two charged pion candidates is required to be less than 3 GeV/c^2 to " \
        "suppress backgrounds from radiative decays to J/psi and chi_c states.")
  .note(:kaon_hypothesis_veto,
        "the event is refitted with the charged tracks assumed to be kaons; if the kaon hypothesis yields " \
        "a smaller chi2 the event is rejected (chi2 stored in the NTuple, veto applied at the ROOT level).")
  .note(:pi0_mass_window,
        "0.11 < m(gamma gamma) < 0.15 GeV/c^2 is required for the pi0 candidate from the kinematic-fit " \
        "corrected four-momenta (applied at the ROOT level).")
  .note(:track_efficiency_correction,
        "tracking efficiency in simulation corrected as a function of polar angle and momentum from a " \
        "tagged J/psi -> 3pi control sample.")
  .note(:pi0_efficiency_correction,
        "pi0 reconstruction efficiency in simulation corrected as a function of pi0 momentum.")
  .note(:reweighting,
        "MC events are reweighted to the data Dalitz distribution; for the psi' sample a comparison with " \
        "a sample generated from amplitudes extracted from a phenomenological fit is used to assign the " \
        "model systematic uncertainty.")
  .note(:continuum_background,
        "continuum background estimated from off-resonance data at 3.08 GeV and 3.650 GeV; 820 +- 55 " \
        "events are expected for the psi' sample.")
  .note(:resonance_background,
        "resonance background estimated from the inclusive MC sample, assigned a 100% uncertainty.")

alg_psip.with_decay_card(decay_card_psip).apply(sel_psip)
root_files_psip = alg_psip.execute_on([psip_data, psip_incMC, exMC_psip])
