# ============================================================================
# arXiv:1505.00539v2 (BESIII)
# Search for the isospin violating decay Y(4260) -> J/psi eta pi0,
# with J/psi -> e+ e- / mu+ mu-, eta -> gamma gamma, pi0 -> gamma gamma.
#
# Data: sqrt(s) = 4.009, 4.226, 4.257, 4.358, 4.416 and 4.599 GeV.
# ISR is simulated with KKMC assuming a Y(4260) Breit-Wigner line shape for
# the Born cross section of e+e- -> J/psi eta pi0; FSR with PHOTOS.
# No signal is observed; only upper limits on the Born cross section are set.
#
# BOSS part only: decay cards, exclusive MC, event selection up to the 4C
# kinematic fit.
# ============================================================================

### Dataset description ###
# Real data at the six centre-of-mass energies used in the search
data_4009 = DatasetManager.real_data.find("703_4009")   # sqrt(s) = 4.009 GeV
data_4226 = DatasetManager.real_data.find("703_4230")   # sqrt(s) = 4.226 GeV
data_4257 = DatasetManager.real_data.find("703_4260")   # sqrt(s) = 4.257 GeV
data_4358 = DatasetManager.real_data.find("703_4360")   # sqrt(s) = 4.358 GeV
data_4416 = DatasetManager.real_data.find("703_4420")   # sqrt(s) = 4.416 GeV
data_4599 = DatasetManager.real_data.find("703_4600")   # sqrt(s) = 4.599 GeV

# Corresponding inclusive MC samples
incMC_4009 = DatasetManager.inclusive_mc.find("703_4009")
incMC_4226 = DatasetManager.inclusive_mc.find("703_4230")
incMC_4257 = DatasetManager.inclusive_mc.find("703_4260")
incMC_4358 = DatasetManager.inclusive_mc.find("703_4360")
incMC_4416 = DatasetManager.inclusive_mc.find("703_4420")
incMC_4599 = DatasetManager.inclusive_mc.find("703_4600")

data_points = [data_4009, data_4226, data_4257, data_4358, data_4416, data_4599]
incMCs      = [incMC_4009, incMC_4226, incMC_4257, incMC_4358, incMC_4416, incMC_4599]

# ----------------------------------------------------------------------------
# Decay cards: one per J/psi lepton channel (EvtGen syntax / EvtGen names).
# The top mother is psi(4260) (BESIII KKMC convention for e+e- annihilation)
# because the Born cross section is assumed to follow the Y(4260) line shape.
# ----------------------------------------------------------------------------
# e+e- -> J/psi eta pi0 with J/psi -> e+ e-
decay_card_ee = <<~DECAYCARD
    Decay psi(4260)
    1.0000 J/psi eta pi0 PHSP;
    Enddecay

    Decay J/psi
    1.0000 e+ e- PHOTOS VLL;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# e+e- -> J/psi eta pi0 with J/psi -> mu+ mu-
decay_card_mumu = <<~DECAYCARD
    Decay psi(4260)
    1.0000 J/psi eta pi0 PHSP;
    Enddecay

    Decay J/psi
    1.0000 mu+ mu- PHOTOS VLL;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# ----------------------------------------------------------------------------
# Exclusive MC: a large signal sample generated uniformly in phase space at
# each centre-of-mass energy, for each J/psi lepton channel.
# ----------------------------------------------------------------------------
exMC_ee = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "JpsiEtaPi0_ee"    # auto-suffixed with the dataset key
  config.events        = 100_000
  config.decay_card    = decay_card_ee
  config.cross_section = :default
end

exMC_mumu = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "JpsiEtaPi0_mumu"
  config.events        = 100_000
  config.decay_card    = decay_card_mumu
  config.cross_section = :default
end

### Event selection (BOSS) ###
# --- J/psi -> e+ e- channel -------------------------------------------------
alg_name_ee = "JpsiEtaPi0EE"
alg_ee = Algorithm.new(alg_name_ee)
alg_ee.set_header(["#{alg_name_ee}Alg/#{alg_name_ee}.h"])
      .set_alias({"std::vector<double>" => "Vdouble"})
      # Multi-energy scan: ECMS is injected per dataset at run time, so no
      # ECMS constant is declared here.

# --- J/psi -> mu+ mu- channel -----------------------------------------------
alg_name_mumu = "JpsiEtaPi0MuMu"
alg_mumu = Algorithm.new(alg_name_mumu)
alg_mumu.set_header(["#{alg_name_mumu}Alg/#{alg_name_mumu}.h"])
        .set_alias({"std::vector<double>" => "Vdouble"})

# Common selection: two charged tracks with zero net charge, at least four
# good photons, lepton identification for J/psi -> l+ l-, eta/pi0
# reconstruction from photon pairs, then the 4C kinematic fit.
event_selection = Selection.new
event_selection.select_track {
                  cos_theta   0.93   # |cos(theta)| < 0.93 in the MDC
                  Vz          10.0   # |Vz| < 10 cm along the beam direction
                  Vr          1.0    # Vr < 1 cm in the plane perpendicular to the beam
                  nChrp       "==1"  # exactly one positively charged track
                  nChrn       "==1"  # exactly one negatively charged track
                  nNet        "==0"  # net charge zero
                }
               .select_photon {
                  tdc_emc_start     0     # EMC timing window
                  tdc_emc_end       14
                  angle_to_track    5.0   # photon at least 5 deg from any charged track
                  energyThreshold_b 0.025 # barrel minimum EMC energy (GeV)
                  energyThreshold_e 0.050 # endcap minimum EMC energy (GeV)
                  nGam              ">=4" # at least four good photon candidates
                }
               .pid(method: :probability) {
                  prob_cut 0.001
                  # Leptons from J/psi decay have momentum above 1 GeV/c;
                  # the EMC energy deposit separates electrons from muons
                  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                                 treat_as_electron_if_energy_above: 0.6
                  nlp "==1"
                  nlm "==1"
                }
               # eta -> gamma gamma and pi0 -> gamma gamma, each reconstructed
               # with a 1C (mass) Kalman fit on a photon pair
               .kalman_kinematic_fit([:gamma, :gamma]) {
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
                  chi2_cut 25
                  neta "==1"
                }
               .kalman_kinematic_fit([:gamma, :gamma]) {
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                  chi2_cut 25
                  npi0 "==1"
                }
               # 4C kinematic fit under the hypothesis e+e- -> gamma gamma gamma gamma l+ l-
               # (the eta and pi0 four-momenta carry the four photons), selecting the
               # photon-pair combination that best matches the eta pi0 hypothesis
               .kinematic_fit([:eta, :pi0, :lp, :lm]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 200   # loose BOSS-pass cut; the paper's chi2 < 40 is applied in ROOT
                }

# Dedicated selection objects for the two lepton channels
selection_ee   = event_selection.dup
selection_mumu = event_selection.dup

# Notes for BOSS-side criteria without a dedicated DSL construct
alg_ee
  .note(:chi2_selection, "The published analysis uses a four-constraint (4C) fit " \
        "under the hypothesis e+e- -> gamma gamma gamma gamma l+ l- with " \
        "chi2_4C < 40; for events with more than four photons the four photons " \
        "returning the smallest chi2 are assigned to the eta and pi0. The nominal " \
        "BOSS pass uses the loose default chi2_cut 200 and the tight value is " \
        "applied at the ROOT stage.")
  .note(:background_veto, "The e+e- -> pi0 pi0 J/psi background is removed by " \
        "rejecting events in which any combination of photon pairs falls in the " \
        "pi0 pi0 region, i.e. both pairs satisfy |M(gamma gamma) - m(pi0)| < " \
        "10 MeV/c^2. This veto is applied on the photon pairs of the 4C-fitted " \
        "event and is therefore implemented at the ROOT stage.")
  .note(:combination_selection, "Among the photon-pair combinations the one " \
        "closest to the eta pi0 signal region is chosen by minimising " \
        "sqrt(|(M(g1 g2) - m_eta)/sigma_eta|^2 + |(M(g3 g4) - m_pi0)/sigma_pi0|^2), " \
        "with sigma_eta and sigma_pi0 the resolutions from signal MC. In the DSL " \
        "this combination choice is delegated to the Kalman fits and to the " \
        "smallest-chi2 combination selection of the kinematic fit.")
  .note(:pid_correction_method, "Electron and muon candidates are separated by E/p " \
        "(E = EMC energy deposit, p = MDC momentum): E/p > 0.7 for electrons and " \
        "E/p < 0.3 for muons. In addition at least one of the two muons is required " \
        "to have at least five layers with valid hits in the MUC, which suppresses " \
        "pion tracks in the final state. These requirements are not expressible in " \
        "the DSL lepton-identification block.")
  .note(:background_veto, "The eta pi0 sidebands are used to estimate the residual " \
        "background: 0.3978 < M(g1 g2) < 0.4578 and 0.6378 < M(g1 g2) < 0.6978 GeV/c^2 " \
        "for the eta, 0.0849 < M(g3 g4) < 0.1049 and 0.1649 < M(g3 g4) < 0.1849 GeV/c^2 " \
        "for the pi0, normalised to the two-dimensional sideband regions.")
  .note(:helix_correction, "The helix parameters of the charged tracks are corrected " \
        "according to the psi(2S) -> gamma chi_cJ control-sample method before the " \
        "kinematic fit; the efficiency difference with and without the correction is " \
        "taken as a systematic uncertainty.")

alg_mumu
  .note(:chi2_selection, "The published analysis uses a four-constraint (4C) fit " \
        "under the hypothesis e+e- -> gamma gamma gamma gamma l+ l- with " \
        "chi2_4C < 40; for events with more than four photons the four photons " \
        "returning the smallest chi2 are assigned to the eta and pi0. The nominal " \
        "BOSS pass uses the loose default chi2_cut 200 and the tight value is " \
        "applied at the ROOT stage.")
  .note(:background_veto, "The e+e- -> pi0 pi0 J/psi background is removed by " \
        "rejecting events in which any combination of photon pairs falls in the " \
        "pi0 pi0 region, i.e. both pairs satisfy |M(gamma gamma) - m(pi0)| < " \
        "10 MeV/c^2. This veto is applied on the photon pairs of the 4C-fitted " \
        "event and is therefore implemented at the ROOT stage.")
  .note(:combination_selection, "Among the photon-pair combinations the one " \
        "closest to the eta pi0 signal region is chosen by minimising " \
        "sqrt(|(M(g1 g2) - m_eta)/sigma_eta|^2 + |(M(g3 g4) - m_pi0)/sigma_pi0|^2), " \
        "with sigma_eta and sigma_pi0 the resolutions from signal MC. In the DSL " \
        "this combination choice is delegated to the Kalman fits and to the " \
        "smallest-chi2 combination selection of the kinematic fit.")
  .note(:pid_correction_method, "Electron and muon candidates are separated by E/p " \
        "(E = EMC energy deposit, p = MDC momentum): E/p > 0.7 for electrons and " \
        "E/p < 0.3 for muons. In addition at least one of the two muons is required " \
        "to have at least five layers with valid hits in the MUC, which suppresses " \
        "pion tracks in the final state. These requirements are not expressible in " \
        "the DSL lepton-identification block.")
  .note(:background_veto, "The eta pi0 sidebands are used to estimate the residual " \
        "background: 0.3978 < M(g1 g2) < 0.4578 and 0.6378 < M(g1 g2) < 0.6978 GeV/c^2 " \
        "for the eta, 0.0849 < M(g3 g4) < 0.1049 and 0.1649 < M(g3 g4) < 0.1849 GeV/c^2 " \
        "for the pi0, normalised to the two-dimensional sideband regions.")
  .note(:helix_correction, "The helix parameters of the charged tracks are corrected " \
        "according to the psi(2S) -> gamma chi_cJ control-sample method before the " \
        "kinematic fit; the efficiency difference with and without the correction is " \
        "taken as a systematic uncertainty.")

# Generate the BOSS algorithm packages (the J/psi, eta and pi0 mass windows and
# the Born cross-section upper limits are evaluated in ROOT)
alg_ee.with_decay_card(decay_card_ee).apply(selection_ee)
alg_mumu.with_decay_card(decay_card_mumu).apply(selection_mumu)

# Execute on real data, inclusive MC and the channel-specific exclusive MC
root_files_ee   = alg_ee.execute_on(data_points + incMCs + exMC_ee)
root_files_mumu = alg_mumu.execute_on(data_points + incMCs + exMC_mumu)
