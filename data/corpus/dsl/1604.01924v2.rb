# arXiv:1604.01924v2 — Measurement of the leptonic decay width of J/psi using ISR
# e+e- -> J/psi gamma -> mu+ mu- gamma (untagged ISR photon), sqrt(s) = 3.773 GeV

### Dataset description ###
psi3770_data  = DatasetManager.real_data.find("712_3773")     # psi(3770) data at sqrt(s) = 3.773 GeV
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")  # Corresponding inclusive MC sample

# Decay card for the signal process e+e- -> mu+ mu- gamma (ISR), full theta_gamma range.
# No intermediate charmonium is imposed at generator level: the di-muon invariant-mass
# spectrum around the J/psi is produced by PHOKHARA through ISR radiation off the beams.
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000 mu+ mu- gamma    PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for the ISR signal process
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_ismumu_gamma"
  config.related_dataset = psi3770_data
  config.events          = 100000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "IsrMuMuGamma"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.773]})
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
  # At least two charged tracks with net charge zero; the pair closest to the IP is used
  # in case of three-track events (DocaVr/DocaVz matching the 1 cm / +-10 cm IP cylinder).
  .select_track do
    cos_theta 0.93        # fiducial volume of the MDC, 0.4 rad < theta < pi - 0.4 rad
    Vz        10.0        # |Vz| < 10 cm along the beam axis
    Vr        1.0         # Vr < 1 cm in the transverse plane
    nChrp     ">=1"       # at least two charged tracks ...
    nChrn     ">=1"
    nTot      ">=2"
    nNet      "==0"       # ... with net charge zero
  end
  # Lepton identification: P(mu) > P(e) for both charged tracks (electron suppression > 96%)
  .pid(method: :probability) do
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    nlp "==1"             # one mu+ candidate
    nlm "==1"             # one mu- candidate
  end
  # 1C kinematic fit: e+e- -> mu+ mu- gamma with a missing massless photon.
  # The constraint is a missing massless particle, imposing overall energy and
  # momentum balance on the two selected charged tracks.
  .kinematic_fit([:lp, :lm]) do
    nominal
    miss_track_of :gamma
    constrain_four_momentum
    chi2_cut 10           # chi2_1C < 10 (published selection)
  end
  # Competing-hypothesis fit (e+e- -> mu+ mu- with no ISR photon), stores chi2 for a
  # ROOT-level consistency check of the 1C hypothesis.
  .kinematic_fit([:lp, :lm]) do
    constrain_four_momentum
  end

my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Procedures that cannot be expressed in the DSL surface are recorded as notes:
my_algorithm
  .note(:pid_correction_method, "Electron suppression via P(mu) > P(e) computed from MDC, TOF and EMC
    information for both charged tracks; hadronic background further rejected with a TMVA Artificial
    Neural Network (yANN < 0.3 for both tracks).")
  .note(:theta_gamma_selection, "The predicted missing ISR photon angle w.r.t. the beam axis must satisfy
    theta_gamma < 0.3 rad or theta_gamma > pi - 0.3 rad in the lab frame (untagged ISR topology); this
    cut is applied on the direction returned by the 1C kinematic fit, i.e. on fit-corrected information.")
  .note(:transverse_momentum_cut, "Transverse momentum p_t > 300 MeV/c required for each of the two
    selected charged tracks.")
  .note(:efficiency_correction, "The true-MC J/psi sample covering the full theta_gamma range is used
    for the efficiency; per-track corrections for muon tracking, electron-PID and ANN efficiency are
    applied. Selection efficiency found to be (32.04 +- 0.09)%.")
  .note(:bin_by_bin_background, "Non-muon backgrounds (pi+pi-gamma, pi+pi-pi0gamma, K+K-gamma, continuum,
    psi(3770) -> non DD, J/psi -> non mu mu ...) are subtracted bin by bin in the 150 m_2mu mass bins
    between 2.8 and 3.4 GeV/c^2.")

root_files = my_algorithm.execute_on([psi3770_data, psi3770_incMC, exMC_signal])
