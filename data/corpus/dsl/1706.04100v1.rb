# =============================================================================
# BESIII: determination of the spin and parity of Z_c(3900)+
# arXiv:1706.04100v1
#
# PWA of e+e- -> pi+ pi- J/psi with J/psi -> l+ l- (l = mu, e) at
# sqrt(s) = 4.23 GeV (1092 pb^-1) and 4.26 GeV (827 pb^-1). The Z_c(3900)+
# spin-parity is determined to be J^P = 1+.
#
# The event selection follows Ref. [1,17]; the amplitude analysis (helicity-
# covariant amplitudes, Flatte-like Z_c lineshape, pi+pi- resonances, unbinned
# maximum-likelihood fit) is a ROOT-level procedure and is captured with notes.
# The J/psi -> mu+mu- and J/psi -> e+e- channels are independent decay modes and
# are handled by separate Algorithm objects (Rule T1).
# =============================================================================

### Dataset description ###
data_4230  = DatasetManager.real_data.find("703_4230")     # sqrt(s) = 4226.26 MeV
data_4260  = DatasetManager.real_data.find("703_4260")     # sqrt(s) = 4257.97 MeV
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")

energy_points = [data_4230, data_4260]

### Decay cards — e+e- -> pi+ pi- J/psi, J/psi -> l+ l- ###
decay_card_mumu = <<~DECAYCARD
    Decay psi(4260)
    1.000  pi+  pi-  J/psi   PHSP;
    Enddecay
    Decay J/psi
    1.000  mu+  mu-   PHSP;
    Enddecay
    End
DECAYCARD

decay_card_ee = <<~DECAYCARD
    Decay psi(4260)
    1.000  pi+  pi-  J/psi   PHSP;
    Enddecay
    Decay J/psi
    1.000  e+  e-   PHSP;
    Enddecay
    End
DECAYCARD

### Exclusive signal MC — both channels generated at each of the two energy points ###
exMC_mumu = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "exmc_pipijpsi_mumu"
  config.events        = 100000
  config.decay_card    = decay_card_mumu
  config.cross_section = :default
end

exMC_ee = DatasetManager.create_exclusive_mc_for(energy_points) do |config|
  config.sample_name   = "exmc_pipijpsi_ee"
  config.events        = 100000
  config.decay_card    = decay_card_ee
  config.cross_section = :default
end

### Event selection (BOSS) ###

# ---------------------------------------------------------------------------
# Channel I: e+e- -> pi+ pi- J/psi, J/psi -> mu+ mu-
# ---------------------------------------------------------------------------
alg_mumu = Algorithm.new("Zc3900PiPiJpsiMuMu")
alg_mumu.set_header(["Zc3900PiPiJpsiMuMuAlg/Zc3900PiPiJpsiMuMu.h"])
        # ECMS is not pinned: the same algorithm runs at the two c.m. energies,
        # the per-run energy being read from the database / injected by jobOptions.
sel_mumu = Selection.new
sel_mumu.select_track do
          cos_theta 0.93    # |cos(theta)| < 0.93
          Vz        10.0    # |Vz| < 10 cm
          Vr        1.0     # Vr < 1 cm
          nChrp     "==2"   # pi+ pi- mu+ mu-
          nChrn     "==2"
          nNet      "==0"
        end
        .pid(method: :probability) do
          # muon identification: small EMC energy deposit relative to momentum, MUC depth
          identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                          treat_as_electron_if_energy_above: 0.6
          nlp "==2"
          nlm "==2"
        end
        .kinematic_fit([:mup, :mum, :pip, :pim]) do
          nominal
          constrain_four_momentum                                # 4C: e+e- -> pi+ pi- J/psi
          invariant_mass_of(:mup, :mum).constrain_to_nominal_mass_of(:jpsi)  # 1C (5C total)
          chi2_cut 200
        end

alg_mumu.with_decay_card(decay_card_mumu).apply(sel_mumu)

# ---------------------------------------------------------------------------
# Channel II: e+e- -> pi+ pi- J/psi, J/psi -> e+ e-
# ---------------------------------------------------------------------------
alg_ee = Algorithm.new("Zc3900PiPiJpsiEE")
alg_ee.set_header(["Zc3900PiPiJpsiEEAlg/Zc3900PiPiJpsiEE.h"])
sel_ee = Selection.new
sel_ee.select_track do
        cos_theta 0.93
        Vz        10.0
        Vr        1.0
        nChrp     "==2"
        nChrn     "==2"
        nNet      "==0"
      end
      .pid(method: :probability) do
        # electron identification: E/p > 0.8 (large EMC energy deposit)
        identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                        treat_as_electron_if_energy_above: 0.6
        nlp "==2"
        nlm "==2"
      end
      .kinematic_fit([:ep, :em, :pip, :pim]) do
        nominal
        constrain_four_momentum                                # 4C: e+e- -> pi+ pi- J/psi
        invariant_mass_of(:ep, :em).constrain_to_nominal_mass_of(:jpsi)   # 1C (5C total)
        chi2_cut 200
      end

alg_ee.with_decay_card(decay_card_ee).apply(sel_ee)

### BOSS-side procedures that have no formal DSL construct ###
[alg_mumu, alg_ee].each do |alg|
  alg.note(:event_selection,
           "the e+e- -> pi+ pi- J/psi selection follows Ref. [1,17]: exactly four " \
           "charged tracks with |cos(theta)| < 0.93, |Vz| < 10 cm and Vr < 1 cm; " \
           "lepton identification with the joint EMC/dE/dx/TOF likelihoods " \
           "(muons: small EMC deposit relative to momentum and MUC depth; " \
           "electrons: E/p > 0.8); the J/psi is required in the " \
           "0.1 < M(l+ l-) < 3.2 GeV/c^2 window and a 4C kinematic fit " \
           "(or the J/psi mass-constrained fit) is applied. Selected candidates: " \
           "4154 (2447) events at 4.23 (4.26) GeV, with 365 (272) background " \
           "events estimated from the J/psi mass sidebands.")
  alg.note(:partial_wave_analysis,
           "the amplitude analysis is performed at the ROOT level: helicity-" \
           "covariant amplitudes with the Z_c(3900)+ and non-Z_c (R J/psi, " \
           "R -> pi+ pi-) components added coherently; pi+ pi- parametrised with " \
           "sigma, f0(980), f2(1270) and f0(1370); the Z_c lineshape is a " \
           "Flatte-like formula with the J/psi pi+/- and (D D*)+/- channels; the " \
           "relative magnitudes and phases of the coupling constants are " \
           "determined by an unbinned maximum-likelihood fit (MINUIT), with the " \
           "backgrounds subtracted from the likelihood; a simultaneous fit is " \
           "performed to the two data sets with the common Z_c parameters.")
  alg.note(:cross_section,
           "sigma = N_Zc / [ L (1+delta) eps B ], with the radiative correction " \
           "factor (1+delta) = 0.818; the Born cross section is " \
           "(22.0 +/- 1.0) pb at 4.23 GeV and (11.0 +/- 1.2) pb at 4.26 GeV. The " \
           "detection efficiency eps is obtained from a MC simulation generated " \
           "with the amplitude parameters determined in the PWA.")
end

### Execution ###
datasets = [data_4230, incMC_4230, data_4260, incMC_4260]

root_files_mumu = alg_mumu.execute_on(datasets + exMC_mumu)
root_files_ee   = alg_ee.execute_on(datasets + exMC_ee)
