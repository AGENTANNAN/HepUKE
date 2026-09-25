# =============================================================================
# Search for the isospin-violating decay Y(4260) -> J/psi eta pi0 at BESIII
#   J/psi -> e+ e-  /  mu+ mu-   (PHOTOS)
#   eta   -> gamma gamma
#   pi0   -> gamma gamma
# BOSS part: dataset preparation + event selection up to the final 4C fit.
# Six XYZ scan points: 4.009, 4.226, 4.257, 4.358, 4.416, 4.599 GeV.
# =============================================================================

### --------------------------- Dataset preparation --------------------------- ###

# Scan points: BOSS sample-name tag -> centre-of-mass energy [GeV]
scan_points = [
  { tag: "4009", ecms: 4.009 },
  { tag: "4230", ecms: 4.226 },
  { tag: "4260", ecms: 4.257 },
  { tag: "4360", ecms: 4.358 },
  { tag: "4420", ecms: 4.416 },
  { tag: "4600", ecms: 4.599 }
]

# Real data and inclusive MC at each of the six energy points
data_points  = scan_points.map { |p| DatasetManager.real_data.find("703_#{p[:tag]}") }
incMC_points = scan_points.map { |p| DatasetManager.inclusive_mc.find("703_#{p[:tag]}") }

# Decay cards (EvtGen syntax) - one per J/psi lepton mode.
# psi(4260) is used as the KKMC top mother (BESIII convention for this process).
decay_card_ee = <<~DECAYCARD
  Decay psi(4260)
  1.0000  J/psi  eta  pi0   PHSP;
  Enddecay

  Decay J/psi
  1.0000  e+  e-   PHOTOS VLL;
  Enddecay

  Decay eta
  1.0000  gamma  gamma   PHSP;
  Enddecay

  Decay pi0
  1.0000  gamma  gamma   PHSP;
  Enddecay

  End
DECAYCARD

decay_card_mumu = <<~DECAYCARD
  Decay psi(4260)
  1.0000  J/psi  eta  pi0   PHSP;
  Enddecay

  Decay J/psi
  1.0000  mu+  mu-   PHOTOS VLL;
  Enddecay

  Decay eta
  1.0000  gamma  gamma   PHSP;
  Enddecay

  Decay pi0
  1.0000  gamma  gamma   PHSP;
  Enddecay

  End
DECAYCARD

# 100k-event exclusive phase-space MC for each J/psi lepton mode at every energy
# point (same signal MC across several distinct energy points -> create_exclusive_mc_for).
exMC_ee = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_y4260_jpsi_eta_pi0_ee"
  config.events        = 100_000
  config.decay_card    = decay_card_ee
  config.cross_section = :default
end

exMC_mumu = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_y4260_jpsi_eta_pi0_mumu"
  config.events        = 100_000
  config.decay_card    = decay_card_mumu
  config.cross_section = :default
end

### ----------------------- Event selection (BOSS) ----------------------- ###
# A single selection chain, shared by both J/psi lepton channels.
event_selection = Selection.new
  .select_track {               # |cosθ| < 0.93, |Vz| < 10 cm, Vr < 1 cm,
    cos_theta 0.93              #   1 positive + 1 negative track, net charge 0
    Vz        10.0
    Vr        1.0
    nChrp     "==1"
    nChrn     "==1"
    nNet      "==0"
  }
  .select_photon {              # >=4 photons, EMC timing 0-14,
    tdc_emc_start     0         # >25 MeV (barrel) / >50 MeV (endcap),
    tdc_emc_end       14        # and >=5 degrees from any charged track
    angle_to_track    5.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=4"
  }
  .pid(method: :probability) {  # probability-method PID, prob > 0.001
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6 # E>0.6 -> e, else mu
    nlp "==1"                   # exactly one l+ and one l-
    nlm "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # 1C fit: gamma gamma -> eta
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta "==1"
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # 1C fit: gamma gamma -> pi0
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 "==1"
  }
  .kinematic_fit([:eta, :pi0, :lp, :lm]) {    # 4C kinematic fit of eta pi0 l+ l-
    nominal                                    # gamma-pair combination with smallest chi2 chosen
    constrain_four_momentum
    chi2_cut 200                               # loose BOSS pass; chi2<40 applied later in ROOT
  }

### ---------------------- Algorithm per energy point ---------------------- ###
# ECMS must match each scan point, so one algorithm is created per point;
# both lepton channels run through the same (shared) selection chain.
scan_points.each_with_index do |p, i|
  alg_name = "JpsiEtaPi0_#{p[:tag]}"
  alg = Algorithm.new(alg_name)
  alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
     .set_constant({"ECMS" => [:double, p[:ecms]]})
  alg.with_decay_card(decay_card_ee).apply(event_selection.dup)
  alg.execute_on([data_points[i], incMC_points[i], exMC_ee[i], exMC_mumu[i]])
end