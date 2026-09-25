# ============================================================================
# Dataset preparation
# ============================================================================
# 14 c.m. energy points of the continuum scan, 4.189 - 4.600 GeV (BOSS 7.0.3)
data_points = %w[
  703_4190 703_4200 703_4210 703_4220 703_4230 703_4237 703_4246
  703_4260 703_4270 703_4280 703_4310 703_4360 703_4420 703_4600
].map { |name| DatasetManager.real_data.find(name) }

# Inclusive MC samples at 4.258, 4.416 and 4.600 GeV
incMC_samples = %w[703_4260 703_4420 703_4600].map { |name| DatasetManager.inclusive_mc.find(name) }

# ---------------------------------------------------------------------------
# Decay cards (EvtGen). Continuum e+e- -> eta' J/psi; psi(4260) is used as the
# KKMC top mother, following the BESIII convention.
# ---------------------------------------------------------------------------
# Mode I : eta' -> gamma pi+ pi-
decay_card_gpipi_ee = <<~DECAYCARD
    Decay psi(4260)
    1.000 eta' J/psi PHSP;
    Enddecay

    Decay eta'
    1.000 gamma pi+ pi- PHSP;
    Enddecay

    Decay J/psi
    1.000 e+ e- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

decay_card_gpipi_mumu = <<~DECAYCARD
    Decay psi(4260)
    1.000 eta' J/psi PHSP;
    Enddecay

    Decay eta'
    1.000 gamma pi+ pi- PHSP;
    Enddecay

    Decay J/psi
    1.000 mu+ mu- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# Mode II : eta' -> eta pi+ pi-, eta -> gamma gamma
decay_card_etapipi_ee = <<~DECAYCARD
    Decay psi(4260)
    1.000 eta' J/psi PHSP;
    Enddecay

    Decay eta'
    1.000 eta pi+ pi- PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    Decay J/psi
    1.000 e+ e- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

decay_card_etapipi_mumu = <<~DECAYCARD
    Decay psi(4260)
    1.000 eta' J/psi PHSP;
    Enddecay

    Decay eta'
    1.000 eta pi+ pi- PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    Decay J/psi
    1.000 mu+ mu- PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# ---------------------------------------------------------------------------
# Exclusive MC: 200k events at every energy point for each of the four final
# states (gamma pi+ pi- e+e-, gamma pi+ pi- mu+mu-, gamma gamma pi+ pi- e+e-,
# gamma gamma pi+ pi- mu+mu-)
# ---------------------------------------------------------------------------
exMC_gpipi_ee = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_etaprime_gpipi_ee"
  config.events        = 200_000
  config.decay_card    = decay_card_gpipi_ee
  config.cross_section = :default
end

exMC_gpipi_mumu = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_etaprime_gpipi_mumu"
  config.events        = 200_000
  config.decay_card    = decay_card_gpipi_mumu
  config.cross_section = :default
end

exMC_etapipi_ee = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_etaprime_etapipi_ee"
  config.events        = 200_000
  config.decay_card    = decay_card_etapipi_ee
  config.cross_section = :default
end

exMC_etapipi_mumu = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_etaprime_etapipi_mumu"
  config.events        = 200_000
  config.decay_card    = decay_card_etapipi_mumu
  config.cross_section = :default
end

# ============================================================================
# Event selection (BOSS)
# ============================================================================
note_pid     = "pion/lepton separation uses p < 0.8 GeV for pions and p > 1.0 GeV for leptons, " \
               "with e/mu separation by E/p (electron: E/p > 0.8, muon: E/p < 0.4). BOSS " \
               "approximates this with identify_high_momentum_leptons (single 1.0 GeV momentum " \
               "threshold and an EMC eraw threshold of 0.8 GeV)."
note_species = "the l+ and l- forming the J/psi candidate must be of the same species " \
               "(e+e- or mu+mu-); the combined index_lp / index_lm lists cannot express this, " \
               "so it is enforced offline / via the per-species exclusive MC."

# ---- Mode I: eta' -> gamma pi+ pi-  (final state gamma pi+ pi- l+ l-) -------
alg_name_gpipi = "EtaPrimeGammaPiPi"
alg_gpipi = Algorithm.new(alg_name_gpipi)
alg_gpipi.set_header(["#{alg_name_gpipi}Alg/#{alg_name_gpipi}.h"])
         .set_constant({"ECMS" => [:double, 4.26]}) # per-run beam energy used by the framework; representative value
         .set_alias({"std::vector<double>" => "Vdouble"})
         .note(:pid_correction_method, note_pid)
         .note(:same_species_lepton_requirement, note_species)

gpipi_selection = Selection.new
    .select_track {
        cos_theta 0.93   # |cos(theta)| < 0.93
        Vz        10.0   # |Vz| < 10 cm
        Vr        1.0    # Vr < 1 cm
        nChrp     "==2"  # two positive tracks
        nChrn     "==2"  # two negative tracks
        nNet      "==0"  # net charge zero
    }
    .select_photon {
        tdc_emc_start     0      # EMC timing window
        tdc_emc_end       14
        angle_to_track    20.0   # >= 20 deg from any charged track
        energyThreshold_b 0.025  # > 25 MeV (barrel)
        energyThreshold_e 0.050  # > 50 MeV (endcap)
        nGam              ">=1"  # at least one photon for the gamma pi+ pi- channels
    }
    .pid(method: :probability) {
        prob_cut 0.001
        identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                       treat_as_electron_if_energy_above: 0.8
        identify :pion, against: [:kaon]
        npip "==1"
        npim "==1"
        nlp  "==1"
        nlm  "==1"
    }
    # 4C kinematic fit to gamma pi+ pi- l+ l- ; tight chi2 < 40 applied offline
    .kinematic_fit([:gamma, :pip, :pim, :lp, :lm]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
    }

alg_gpipi.with_decay_card(decay_card_gpipi_ee).apply(gpipi_selection)
alg_gpipi.execute_on(data_points + incMC_samples + exMC_gpipi_ee + exMC_gpipi_mumu)

# ---- Mode II: eta' -> eta pi+ pi-, eta -> gamma gamma ----------------------
alg_name_etapipi = "EtaPrimeEtaPiPi"
alg_etapipi = Algorithm.new(alg_name_etapipi)
alg_etapipi.set_header(["#{alg_name_etapipi}Alg/#{alg_name_etapipi}.h"])
          .set_constant({"ECMS" => [:double, 4.26]})
          .set_alias({"std::vector<double>" => "Vdouble"})
          .note(:pid_correction_method, note_pid)
          .note(:same_species_lepton_requirement, note_species)

etapipi_selection = Selection.new
    .select_track {
        cos_theta 0.93   # |cos(theta)| < 0.93
        Vz        10.0   # |Vz| < 10 cm
        Vr        1.0    # Vr < 1 cm
        nChrp     "==2"  # two positive tracks
        nChrn     "==2"  # two negative tracks
        nNet      "==0"  # net charge zero
    }
    .select_photon {
        tdc_emc_start     0      # EMC timing window
        tdc_emc_end       14
        angle_to_track    20.0   # >= 20 deg from any charged track
        energyThreshold_b 0.025  # > 25 MeV (barrel)
        energyThreshold_e 0.050  # > 50 MeV (endcap)
        nGam              ">=2"  # two photons from eta -> gamma gamma
    }
    .pid(method: :probability) {
        prob_cut 0.001
        identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                       treat_as_electron_if_energy_above: 0.8
        identify :pion, against: [:kaon]
        npip "==1"
        npim "==1"
        nlp  "==1"
        nlm  "==1"
    }
    # 5C kinematic fit to gamma gamma pi+ pi- l+ l- : 4C + M(gamma gamma) -> M(eta)
    .kinematic_fit([:gamma, :gamma, :pip, :pim, :lp, :lm]) {
        nominal
        constrain_four_momentum
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
        chi2_cut 200
    }

alg_etapipi.with_decay_card(decay_card_etapipi_ee).apply(etapipi_selection)
alg_etapipi.execute_on(data_points + incMC_samples + exMC_etapipi_ee + exMC_etapipi_mumu)