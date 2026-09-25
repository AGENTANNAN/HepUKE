# =====================================================================
# BESIII: Observation of e+e- -> gamma X(3872)
#         X(3872) -> pi+pi- J/psi,  J/psi -> l+l-  (l = e, mu)
# Data: sqrt(s) = 4.009, 4.229, 4.260, 4.360 GeV
# =====================================================================

### Dataset description ###
# Real data at the four center-of-mass energies used in the paper
data_4009 = DatasetManager.real_data.find("703_4009")   # 4.009 GeV
data_4229 = DatasetManager.real_data.find("703_4230")   # 4.229 GeV
data_4260 = DatasetManager.real_data.find("703_4260")   # 4.260 GeV
data_4360 = DatasetManager.real_data.find("703_4360")   # 4.360 GeV

# Corresponding inclusive MC samples (used for background study)
incMC_4009 = DatasetManager.inclusive_mc.find("703_4009")
incMC_4229 = DatasetManager.inclusive_mc.find("703_4230")
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")
incMC_4360 = DatasetManager.inclusive_mc.find("703_4360")

data_points  = [data_4009, data_4229, data_4260, data_4360]
incMC_points = [incMC_4009, incMC_4229, incMC_4260, incMC_4360]

# ---------------------------------------------------------------------
# Signal decay card: e+e- -> gamma X(3872), X(3872) -> rho0 J/psi,
# rho0 -> pi+pi-.  The radiative transition is generated with the
# psi(4260) top mother (BESIII KKMC/ISR convention); the Born cross
# section follows the e+e- -> pi+pi- J/psi line shape and the maximum
# ISR photon energy corresponds to the 3.9 GeV gamma X(3872) threshold.
# rho0 and J/psi are in a relative S-wave.
# ---------------------------------------------------------------------
# Mode I: J/psi -> e+ e-
decay_card_signal_ee = <<~DECAYCARD
    Decay psi(4260)
    1.0000 gamma X(3872) PHSP;
    Enddecay

    Decay X(3872)
    1.0000 rho0 J/psi PHSP;
    Enddecay

    Decay rho0
    1.0000 pi+ pi- PHSP;
    Enddecay

    Decay J/psi
    1.0000 e+ e- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# Mode II: J/psi -> mu+ mu-
decay_card_signal_mumu = <<~DECAYCARD
    Decay psi(4260)
    1.0000 gamma X(3872) PHSP;
    Enddecay

    Decay X(3872)
    1.0000 rho0 J/psi PHSP;
    Enddecay

    Decay rho0
    1.0000 pi+ pi- PHSP;
    Enddecay

    Decay J/psi
    1.0000 mu+ mu- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# ---------------------------------------------------------------------
# Background decay card explicitly simulated in the paper:
# e+e- -> (gamma_ISR) pi+pi- J/psi, which is the dominant peaking
# background under the X(3872) signal (the ISR photon is faked by the
# radiative photon).
# ---------------------------------------------------------------------
decay_card_bkg_pipiJpsi = <<~DECAYCARD
    Decay psi(4260)
    1.0000 pi+ pi- J/psi PHSP;
    Enddecay

    Decay J/psi
    1.0000 e+ e- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# Signal exclusive MC, one sample per center-of-mass energy point
exMCs_signal_ee = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_gammaX3872_pipijpsi_ee"    # auto-suffixed per energy point
  config.events        = 100_000
  config.decay_card    = decay_card_signal_ee
  config.cross_section = :default
end

exMCs_signal_mumu = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_gammaX3872_pipijpsi_mumu"  # auto-suffixed per energy point
  config.events        = 100_000
  config.decay_card    = decay_card_signal_mumu
  config.cross_section = :default
end

exMCs_bkg_pipiJpsi = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_pipijpsi_isr_bkg"          # auto-suffixed per energy point
  config.events        = 100_000
  config.decay_card    = decay_card_bkg_pipiJpsi
  config.cross_section = :default
end

# ---------------------------------------------------------------------
### Event selection (BOSS) ###
# ---------------------------------------------------------------------
alg_name = "gammaX3872"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
# Four good charged tracks with net charge zero
event_selection.select_track do
                  cos_theta 0.93      # |cos(theta)| < 0.93
                  Vz        10.0      # |Vz| < 10 cm
                  Vr        1.0       # Vr < 1 cm
                  nChrp     "==2"     # exactly two positive tracks (pi+/l+)
                  nChrn     "==2"     # exactly two negative tracks (pi-/l-)
                  nNet      "==0"     # net charge zero
                end
               # Good photon candidates (fiducial, shower quality, timing)
               .select_photon do
                  tdc_emc_start     0
                  tdc_emc_end       14
                  angle_to_track    10.0
                  energyThreshold_b 0.025
                  energyThreshold_e 0.050
                  nGam              ">=1"   # at least one radiative photon candidate
                end
               # Two pions and two leptons among the four charged tracks
               .pid(method: :probability) do
                  prob_cut 0.001
                  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                                 treat_as_electron_if_energy_above: 0.6
                  identify :pion, against: [:kaon, :proton]
                  npip "==1"
                  npim "==1"
                  nlp  "==1"
                  nlm  "==1"
                end
               # 4C kinematic fit to the hypothesis e+e- -> gamma pi+pi- l+l-,
               # i.e. total four-momentum of the measured particles equals the
               # initial four-momentum of the colliding beams.
               .kinematic_fit([:gamma, :pip, :pim, :lp, :lm]) do
                  nominal
                  constrain_four_momentum
                  chi2_cut 60   # paper: chi2 of the 4C fit < 60
                  # Reject e+e- -> eta J/psi (eta -> gamma pi+pi- / pi+pi-pi0)
                  invariant_mass_of(:gamma, :pip, :pim).within(0.6, 20.0)
                  # J/psi mass window for signal (3.08 < M(l+l-) < 3.12 GeV/c^2);
                  # the J/psi sidebands are defined in the ROOT analysis
                  invariant_mass_of(:lp, :lm).within(3.08, 3.12)
                end
               # The remaining selection criteria cannot be expressed as BOSS
               # constructs and are recorded for downstream handling.
               .note(:radiative_photon_selection,
                     "the largest-energy photon candidate is taken as the radiative photon; the 4C fit iterates over all photon candidates and keeps the smallest chi2 combination")
               .note(:background_veto,
                     "|cos(opening angle between the two pion candidates)| < 0.98, applied to reject radiative Bhabha (gamma e+e-) and radiative dimuon (gamma mu+mu-) backgrounds with photon conversion; efficiency loss for signal < 1%")

# Generate the algorithm for the signal decay card
my_algorithm.with_decay_card(decay_card_signal_ee).apply(event_selection)

# Execute on real data, inclusive MC, and all exclusive MC samples
root_files = my_algorithm.execute_on(data_points + incMC_points +
                                     exMCs_signal_ee + exMCs_signal_mumu +
                                     exMCs_bkg_pipiJpsi)
