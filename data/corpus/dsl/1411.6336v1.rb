# BESIII arXiv:1411.6336v1
# Search for e+e- -> gamma chi_cJ (J=0,1,2), chi_cJ -> gamma J/psi, J/psi -> mu+mu-
# Data at sqrt(s) = 4.009, 4.230, 4.260, 4.360 GeV
# Final state: gamma gamma mu+ mu-  (J/psi -> e+e- deliberately not used due to Bhabha background)

### Datasets — four CME points ###
energies = ["4009", "4230", "4260", "4360"]
data_samples  = energies.map { |e| DatasetManager.real_data.find("703_#{e}") }
incMC_samples = energies.map { |e| DatasetManager.inclusive_mc.find("703_#{e}") }

### Decay cards — one per chi_cJ state (identical final state, different card) ###
decay_card_chic0 = <<~DECAYCARD
  Decay psi(4260)
  1.000 gamma chi_c0                              PHSP;
  Enddecay

  Decay chi_c0
  1.000 gamma J/psi                               PHSP;
  Enddecay

  Decay J/psi
  1.000 mu+ mu-                                   PHOTOS VLL;
  Enddecay

  End
DECAYCARD

decay_card_chic1 = <<~DECAYCARD
  Decay psi(4260)
  1.000 gamma chi_c1                              PHSP;
  Enddecay

  Decay chi_c1
  1.000 gamma J/psi                               PHSP;
  Enddecay

  Decay J/psi
  1.000 mu+ mu-                                   PHOTOS VLL;
  Enddecay

  End
DECAYCARD

decay_card_chic2 = <<~DECAYCARD
  Decay psi(4260)
  1.000 gamma chi_c2                              PHSP;
  Enddecay

  Decay chi_c2
  1.000 gamma J/psi                               PHSP;
  Enddecay

  Decay J/psi
  1.000 mu+ mu-                                   PHOTOS VLL;
  Enddecay

  End
DECAYCARD

### Exclusive signal MC — one sample per chi_cJ per CME point ###
# Generator assumption: E1 transition for psi(4260) -> gamma chi_cJ; ISR effects via KKMC
# with the Y(4260) described by a PDG Breit-Wigner (top mother = psi(4260)).
exMC_chic0 = DatasetManager.create_exclusive_mc_for(data_samples) do |config|
  config.sample_name   = "gammachic0"
  config.events        = 200000
  config.decay_card    = decay_card_chic0
  config.cross_section = :default
end

exMC_chic1 = DatasetManager.create_exclusive_mc_for(data_samples) do |config|
  config.sample_name   = "gammachic1"
  config.events        = 200000
  config.decay_card    = decay_card_chic1
  config.cross_section = :default
end

exMC_chic2 = DatasetManager.create_exclusive_mc_for(data_samples) do |config|
  config.sample_name   = "gammachic2"
  config.events        = 200000
  config.decay_card    = decay_card_chic2
  config.cross_section = :default
end

### Common event selection (shared by chi_c0, chi_c1 and chi_c2 — identical final state) ###
common_selection = Selection.new
common_selection.select_track {
    cos_theta 0.93        # |cos(theta)| < 0.93 for each good charged track
    Vz        10.0        # |Vz| < 10 cm (beam direction)
    Vr        1.0         # |Vr| < 1 cm (transverse plane)
    nChrp     "==2"       # exactly two positively charged tracks
    nChrn     "==2"       # exactly two negatively charged tracks
    nNet      "==0"       # zero net charge
  }
  .select_photon {
    angle_to_track    20.0   # photon at least 20 degrees away from any charged track
    energyThreshold_b 0.025  # > 25 MeV in the EMC barrel (|cos(theta)| < 0.80)
    energyThreshold_e 0.050  # > 50 MeV in the EMC end-cap (0.86 < |cos(theta)| < 0.92)
    tdc_emc_start     0      # EMC time window 0-700 ns with the collision
    tdc_emc_end       14
    nGam              ">=2"  # at least two photon candidates
  }
  .pid(method: :probability) {
    # Both charged tracks are muons from J/psi -> mu+ mu- (high-momentum leptons)
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                    treat_as_electron_if_energy_above: 0.6
    nlp "==1"
    nlm "==1"
  }

### Algorithm 1: e+e- -> gamma chi_c0, chi_c0 -> gamma J/psi, J/psi -> mu+ mu- ###
alg_chic0 = Algorithm.new("GammaChic0")
alg_chic0.set_header(["GammaChic0Alg/GammaChic0.h"])
         .set_constant({"ECMS" => [:double, 4.260]})  # same selection run at 4.009/4.230/4.260/4.360 GeV

sel_chic0 = common_selection.dup
  # 5C kinematic fit: e+e- -> gamma gamma mu+ mu- with 4-momentum conservation
  # and M(mu+mu-) constrained to the nominal J/psi mass
  .kinematic_fit([:gamma, :gamma, :lp, :lm]) {
    nominal
    constrain_four_momentum
    invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)
    chi2_cut 200   # loose BOSS cut; published chi2_5C < 40 applied in ROOT
  }

alg_chic0.note(:muon_pid,
  "Muons are identified by E/p < 0.35 and p > 1.0 GeV/c, where E is the EMC energy " \
  "deposit and p the MDC momentum. The DSL lepton declaration exposes only the " \
  "momentum / EMC-energy thresholds (1.0 GeV/c, 0.6 GeV); the E/p < 0.35 requirement " \
  "must be applied as an additional track-level criterion.")
         .note(:background_veto,
  "Backgrounds from e+e- -> P J/psi with P -> gamma gamma (P = pi0, eta, eta') are " \
  "rejected by requiring |M(gamma gamma) - M(pi0)| > 0.025 GeV/c^2, " \
  "|M(gamma gamma) - M(eta)| > 0.03 GeV/c^2 and |M(gamma gamma) - M(eta')| > 0.02 GeV/c^2, " \
  "where M(gamma gamma) is the invariant mass of the two selected photons.")
         .note(:photon_energy_selection,
  "At sqrt(s) = 4.009 GeV the photon energy spectra from e+e- -> gamma chi_c1,2 and from " \
  "chi_c1,2 -> gamma J/psi overlap; the energy of the photon from the chi_c1,2 decay is " \
  "required to be below 0.403 GeV to separate them.")
         .note(:multi_candidate,
  "If more than one candidate survives in an event, the candidate with the smallest " \
  "chi2_5C is retained; this is handled by the combination ranking of the kinematic fit.")
         .note(:jpsi_sideband,
  "Control samples from J/psi sidebands 2.917 < M(mu+mu-) < 3.057 GeV/c^2 and " \
  "3.137 < M(mu+mu-) < 3.277 GeV/c^2 are studied by constraining M(mu+mu-) to 3.047 or " \
  "3.147 GeV/c^2 in the 5C fit; no chi_cJ signal is expected there.")

alg_chic0.with_decay_card(decay_card_chic0).apply(sel_chic0)
alg_chic0.execute_on(data_samples + incMC_samples + exMC_chic0)

### Algorithm 2: e+e- -> gamma chi_c1, chi_c1 -> gamma J/psi, J/psi -> mu+ mu- ###
alg_chic1 = Algorithm.new("GammaChic1")
alg_chic1.set_header(["GammaChic1Alg/GammaChic1.h"])
         .set_constant({"ECMS" => [:double, 4.260]})

sel_chic1 = common_selection.dup
  .kinematic_fit([:gamma, :gamma, :lp, :lm]) {
    nominal
    constrain_four_momentum
    invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)
    chi2_cut 200
  }

alg_chic1.note(:muon_pid,
  "Same muon identification as in the gamma chi_c0 channel (E/p < 0.35 and p > 1.0 GeV/c).")
         .note(:background_veto,
  "Same pi0 / eta / eta' vetoes on M(gamma gamma) as in the gamma chi_c0 channel.")
         .note(:photon_energy_selection,
  "At sqrt(s) = 4.009 GeV the photon from chi_c1 -> gamma J/psi is required to have " \
  "energy below 0.403 GeV to resolve the overlap with the prompt photon.")
         .note(:multi_candidate,
  "Smallest chi2_5C candidate selected per event.")

alg_chic1.with_decay_card(decay_card_chic1).apply(sel_chic1)
alg_chic1.execute_on(data_samples + incMC_samples + exMC_chic1)

### Algorithm 3: e+e- -> gamma chi_c2, chi_c2 -> gamma J/psi, J/psi -> mu+ mu- ###
alg_chic2 = Algorithm.new("GammaChic2")
alg_chic2.set_header(["GammaChic2Alg/GammaChic2.h"])
         .set_constant({"ECMS" => [:double, 4.260]})

sel_chic2 = common_selection.dup
  .kinematic_fit([:gamma, :gamma, :lp, :lm]) {
    nominal
    constrain_four_momentum
    invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)
    chi2_cut 200
  }

alg_chic2.note(:muon_pid,
  "Same muon identification as in the gamma chi_c0 channel (E/p < 0.35 and p > 1.0 GeV/c).")
         .note(:background_veto,
  "Same pi0 / eta / eta' vetoes on M(gamma gamma) as in the gamma chi_c0 channel.")
         .note(:photon_energy_selection,
  "At sqrt(s) = 4.009 GeV the photon from chi_c2 -> gamma J/psi is required to have " \
  "energy below 0.403 GeV to resolve the overlap with the prompt photon.")
         .note(:multi_candidate,
  "Smallest chi2_5C candidate selected per event.")

alg_chic2.with_decay_card(decay_card_chic2).apply(sel_chic2)
alg_chic2.execute_on(data_samples + incMC_samples + exMC_chic2)
