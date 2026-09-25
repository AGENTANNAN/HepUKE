# ============================================================================
# arXiv:1503.08203v2 (BESIII)
# Observation of X(3823) in e+e- -> pi+ pi- X(3823),
# with X(3823) -> gamma chi_c1, gamma chi_c2, chi_c1,c2 -> gamma J/psi,
# J/psi -> l+ l- (l = e, mu).
#
# Data: 4.67 fb^-1 at sqrt(s) = 4.230, 4.260, 4.360, 4.420, 4.600 GeV.
# ISR is simulated with KKMC, FSR with PHOTOS; the Born cross section of
# e+e- -> pi+ pi- X(3823) is assumed to follow the e+e- -> pi+ pi- psi(2S)
# lineshape, hence the psi(4260) top mother (BESIII convention) is used.
#
# BOSS part only: decay cards, exclusive MC, event selection up to the 4C
# kinematic fit.
# ============================================================================

### Dataset description ###
# Real data at the five centre-of-mass energies used in the analysis
data_4230 = DatasetManager.real_data.find("703_4230")
data_4260 = DatasetManager.real_data.find("703_4260")
data_4360 = DatasetManager.real_data.find("703_4360")
data_4420 = DatasetManager.real_data.find("703_4420")
data_4600 = DatasetManager.real_data.find("703_4600")

# Corresponding inclusive MC samples
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")
incMC_4360 = DatasetManager.inclusive_mc.find("703_4360")
incMC_4420 = DatasetManager.inclusive_mc.find("703_4420")
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")

data_points = [data_4230, data_4260, data_4360, data_4420, data_4600]
incMCs      = [incMC_4230, incMC_4260, incMC_4360, incMC_4420, incMC_4600]

# ----------------------------------------------------------------------------
# Decay cards: one card per reconstructed decay mode of X(3823) and per
# lepton channel of J/psi -> l+ l- (EvtGen particle names, EvtGen syntax).
# ----------------------------------------------------------------------------
# X(3823) -> gamma chi_c1, chi_c1 -> gamma J/psi, J/psi -> e+ e-
decay_card_c1_e = <<~DECAYCARD
    Decay psi(4260)
    1.0000 pi+ pi- X(3823) PHSP;
    Enddecay

    Decay X(3823)
    1.0000 gamma chi_c1 PHSP;
    Enddecay

    Decay chi_c1
    1.0000 gamma J/psi PHSP;
    Enddecay

    Decay J/psi
    1.0000 e+ e- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# X(3823) -> gamma chi_c1, chi_c1 -> gamma J/psi, J/psi -> mu+ mu-
decay_card_c1_mu = <<~DECAYCARD
    Decay psi(4260)
    1.0000 pi+ pi- X(3823) PHSP;
    Enddecay

    Decay X(3823)
    1.0000 gamma chi_c1 PHSP;
    Enddecay

    Decay chi_c1
    1.0000 gamma J/psi PHSP;
    Enddecay

    Decay J/psi
    1.0000 mu+ mu- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# X(3823) -> gamma chi_c2, chi_c2 -> gamma J/psi, J/psi -> e+ e-
decay_card_c2_e = <<~DECAYCARD
    Decay psi(4260)
    1.0000 pi+ pi- X(3823) PHSP;
    Enddecay

    Decay X(3823)
    1.0000 gamma chi_c2 PHSP;
    Enddecay

    Decay chi_c2
    1.0000 gamma J/psi PHSP;
    Enddecay

    Decay J/psi
    1.0000 e+ e- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# X(3823) -> gamma chi_c2, chi_c2 -> gamma J/psi, J/psi -> mu+ mu-
decay_card_c2_mu = <<~DECAYCARD
    Decay psi(4260)
    1.0000 pi+ pi- X(3823) PHSP;
    Enddecay

    Decay X(3823)
    1.0000 gamma chi_c2 PHSP;
    Enddecay

    Decay chi_c2
    1.0000 gamma J/psi PHSP;
    Enddecay

    Decay J/psi
    1.0000 mu+ mu- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# ----------------------------------------------------------------------------
# Exclusive MC: 40,000 signal events generated at each centre-of-mass energy
# (one sample per energy point and per decay mode/channel).
# ----------------------------------------------------------------------------
exMC_c1_e = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "X3823_gammac1_ee"   # auto-suffixed with the dataset key
  config.events        = 40_000
  config.decay_card    = decay_card_c1_e
  config.cross_section = :default
end

exMC_c1_mu = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "X3823_gammac1_mumu"
  config.events        = 40_000
  config.decay_card    = decay_card_c1_mu
  config.cross_section = :default
end

exMC_c2_e = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "X3823_gammac2_ee"
  config.events        = 40_000
  config.decay_card    = decay_card_c2_e
  config.cross_section = :default
end

exMC_c2_mu = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "X3823_gammac2_mumu"
  config.events        = 40_000
  config.decay_card    = decay_card_c2_mu
  config.cross_section = :default
end

# Signal MC for the chi_c1 mode (both lepton channels) and for the chi_c2 mode
exMC_chi_c1 = exMC_c1_e + exMC_c1_mu
exMC_chi_c2 = exMC_c2_e + exMC_c2_mu

### Event selection (BOSS) ###
# The gamma chi_c1 and gamma chi_c2 signals share the final state
# pi+ pi- gamma gamma l+ l- and identical selection criteria; the two
# modes are separated only by the M(gamma_H J/psi) window applied later
# in the ROOT analysis. One Algorithm per mode is kept for bookkeeping,
# with a common selection chain.

# --- gamma chi_c1 mode ------------------------------------------------------
alg_name_c1 = "X3823ToGammaChiC1"
alg_c1 = Algorithm.new(alg_name_c1)
alg_c1.set_header(["#{alg_name_c1}Alg/#{alg_name_c1}.h"])
       .set_alias({"std::vector<double>" => "Vdouble"})
       # Multi-energy scan: ECMS is injected per dataset at run time, so no
       # ECMS constant is declared here.

# --- gamma chi_c2 mode ------------------------------------------------------
alg_name_c2 = "X3823ToGammaChiC2"
alg_c2 = Algorithm.new(alg_name_c2)
alg_c2.set_header(["#{alg_name_c2}Alg/#{alg_name_c2}.h"])
       .set_alias({"std::vector<double>" => "Vdouble"})

# Common selection: four charged tracks with zero net charge, at least two
# good photons, e/mu identification for J/psi -> l+ l-, then the 4C fit to
# e+e- -> pi+ pi- gamma gamma l+ l-.
event_selection = Selection.new
event_selection.select_track {
                  cos_theta   0.93   # |cos(theta)| < 0.93 for charged tracks
                  Vz          10.0   # |Vz| < 10 cm
                  Vr          1.0    # Vr < 1 cm
                  nChrp       "==2"  # exactly two positively charged tracks
                  nChrn       "==2"  # exactly two negatively charged tracks
                  nNet        "==0"  # zero net charge
                }
               .select_photon {
                  tdc_emc_start     0     # EMC timing window
                  tdc_emc_end       14
                  angle_to_track    10.0  # min angle to the nearest charged track (deg)
                  energyThreshold_b 0.025 # barrel energy threshold (GeV)
                  energyThreshold_e 0.050 # endcap energy threshold (GeV)
                  nGam              ">=2" # at least two good photon candidates
                }
               .pid(method: :probability) {
                  prob_cut 0.001
                  # J/psi -> l+ l- (l = e or mu): high momentum leptons
                  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                                 treat_as_electron_if_energy_above: 0.6
                  identify :pion, against: [:kaon, :proton]  # the two low momentum pions
                  npip "==1"
                  npim "==1"
                  nlp  "==1"
                  nlm  "==1"
                }
               # 4C kinematic fit: constrain the total four-momentum of the
               # detected particles to the initial four-momentum of the beams
               # (hypothesis e+e- -> pi+ pi- gamma gamma l+ l-).
               .kinematic_fit([:pip, :pim, :lp, :lm, :gamma, :gamma]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 200   # loose BOSS-pass cut; the paper's chi2 < 80 is applied in ROOT
                }

# Dedicated selection objects for the two modes
selection_c1 = event_selection.dup
selection_c2 = event_selection.dup

# Notes for BOSS-side criteria that have no dedicated DSL construct (all of
# them are applied downstream of the 4C fit, on fit-corrected quantities)
alg_c1
  .note(:chi2_selection, "The published selection requires chi2_4C < 80 for the " \
        "hypothesis e+e- -> pi+ pi- gamma gamma l+ l- (about 95% efficient for " \
        "signal). The nominal BOSS pass uses the loose default chi2_cut 200; the " \
        "tight value is applied at the ROOT stage.")
  .note(:background_veto, "Radiative Bhabha / radiative dimuon (gamma e+e- / " \
        "gamma mu+ mu-) background from photon conversion is rejected by requiring " \
        "cos(opening angle of the pi+ pi- pair) < 0.98 (efficiency loss < 1% for signal).")
  .note(:background_veto_etajpsi, "e+e- -> eta J/psi with eta -> pi+ pi- pi0 / " \
        "gamma pi+ pi- is rejected by M(gamma gamma pi+ pi-) > 0.57 GeV/c^2; " \
        "e+e- -> gamma_ISR psi(2S), e+e- -> eta psi(2S) and e+e- -> gamma gamma " \
        "psi(2S) are rejected by |M(pi+ pi- J/psi) - m(psi(2S))| > 6 MeV/c^2.")
  .note(:jpsi_mass_window, "The J/psi candidate is selected with " \
        "3.08 < M(l+ l-) < 3.13 GeV/c^2 (resolution 9 MeV/c^2 from MC), with " \
        "sidebands 3.01-3.06 and 3.15-3.20 GeV/c^2 used to model the non-J/psi " \
        "background. The two photons are ranked by energy: the higher-energy " \
        "photon (gamma_H) forms chi_c1,c2 via 3.490 < M(gamma_H J/psi) < 3.530 GeV/c^2 " \
        "(chi_c1) and 3.536 < M(gamma_H J/psi) < 3.576 GeV/c^2 (chi_c2); the " \
        "lower-energy photon is assigned to the X(3823) decay.")
  .note(:efficiency_curve, "The Born cross section of e+e- -> pi+ pi- X(3823) is " \
        "assumed to follow the e+e- -> pi+ pi- psi(2S) lineshape between 4.1 and " \
        "4.6 GeV; the maximum ISR photon energy corresponds to the 4.1 GeV/c^2 " \
        "production threshold of the pi+ pi- X(3823) system. The resulting " \
        "radiative correction factor (1+delta), the vacuum polarisation factor " \
        "1/|1-Pi|^2 (0.5% from QED) and the detection efficiency enter the " \
        "cross-section measurement.")

alg_c2
  .note(:chi2_selection, "The published selection requires chi2_4C < 80 for the " \
        "hypothesis e+e- -> pi+ pi- gamma gamma l+ l- (about 95% efficient for " \
        "signal). The nominal BOSS pass uses the loose default chi2_cut 200; the " \
        "tight value is applied at the ROOT stage.")
  .note(:background_veto, "Radiative Bhabha / radiative dimuon (gamma e+e- / " \
        "gamma mu+ mu-) background from photon conversion is rejected by requiring " \
        "cos(opening angle of the pi+ pi- pair) < 0.98 (efficiency loss < 1% for signal).")
  .note(:background_veto_etajpsi, "e+e- -> eta J/psi with eta -> pi+ pi- pi0 / " \
        "gamma pi+ pi- is rejected by M(gamma gamma pi+ pi-) > 0.57 GeV/c^2; " \
        "e+e- -> gamma_ISR psi(2S), e+e- -> eta psi(2S) and e+e- -> gamma gamma " \
        "psi(2S) are rejected by |M(pi+ pi- J/psi) - m(psi(2S))| > 6 MeV/c^2.")
  .note(:jpsi_mass_window, "The J/psi candidate is selected with " \
        "3.08 < M(l+ l-) < 3.13 GeV/c^2 (resolution 9 MeV/c^2 from MC), with " \
        "sidebands 3.01-3.06 and 3.15-3.20 GeV/c^2 used to model the non-J/psi " \
        "background. The two photons are ranked by energy: the higher-energy " \
        "photon (gamma_H) forms chi_c1,c2 via 3.490 < M(gamma_H J/psi) < 3.530 GeV/c^2 " \
        "(chi_c1) and 3.536 < M(gamma_H J/psi) < 3.576 GeV/c^2 (chi_c2); the " \
        "lower-energy photon is assigned to the X(3823) decay.")
  .note(:efficiency_curve, "No X(3823) signal is observed in the gamma chi_c2 mode; " \
        "only an upper limit on its production rate is reported. The same efficiency " \
        "and radiative-correction treatment as for the gamma chi_c1 mode applies.")

# Generate the BOSS algorithm packages (the X(3823) recoil mass and the
# M(gamma_H J/psi) regions are analysed on the corrected four-momenta in ROOT)
alg_c1.with_decay_card(decay_card_c1_e).apply(selection_c1)
alg_c2.with_decay_card(decay_card_c2_e).apply(selection_c2)

# Execute on real data, inclusive MC and the mode-specific exclusive MC
root_files_chi_c1 = alg_c1.execute_on(data_points + incMCs + exMC_chi_c1)
root_files_chi_c2 = alg_c2.execute_on(data_points + incMCs + exMC_chi_c2)
