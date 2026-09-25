# Datasets: real data and inclusive MC for the specified BOSS versions in the
# continuum range 3.508–4.951 GeV (41 energy points, 26.1 fb⁻¹ total).
data_points = DatasetManager.real_data
  .where(cms_energy: { value: 3508..4951 })
  .all
  .select { |d| %w[703 705 706 709 712].include?(d.boss) }

incMC_points = DatasetManager.inclusive_mc
  .where(cms_energy: { value: 3508..4951 })
  .all
  .select { |d| %w[703 705 706 709 712].include?(d.boss) }

# Decay card for Mode I: e+e- -> phi eta', phi -> K+ K-, eta' -> gamma pi+ pi-
decay_card_modeI = <<~DECAYCARD
  Decay psi(4260)
  1.0000 phi eta' PHSP;
  Enddecay

  Decay phi
  1.000 K+ K- VSS;
  Enddecay

  Decay eta'
  1.000 gamma pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

# Decay card for Mode II: e+e- -> phi eta', phi -> K+ K-,
# eta' -> eta pi+ pi-, eta -> gamma gamma
decay_card_modeII = <<~DECAYCARD
  Decay psi(4260)
  1.0000 phi eta' PHSP;
  Enddecay

  Decay phi
  1.000 K+ K- VSS;
  Enddecay

  Decay eta'
  1.000 eta pi+ pi- PHSP;
  Enddecay

  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Exclusive MC for the continuum phi eta' signal with eta' -> gamma pi+ pi-
# (Mode I). 100k events are generated for each energy point in the scan.
exMC_modeI = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_phi_etaprime_gammapipim"
  config.events        = 100_000
  config.decay_card    = decay_card_modeI
  config.cross_section = :default
end

# Common base selection shared by both eta' decay modes:
#   - charged tracks: exactly two positive and two negative, net charge zero,
#     |cosθ| < 0.93, |Vz| < 10 cm, Vr < 1 cm
#   - photons: TDC 0–14 (0–700 ns), angle to any track > 10°,
#     E > 25 MeV (barrel) / 50 MeV (endcap), at least one photon
#   - PID: probability method, separate K from π, require one K+ and one K−
#   - remove identified kaons, then assign the remaining tracks as π+ and π−
base_selection = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==2"
    nChrn     "==2"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=1"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    nkp      "==1"
    nkm      "==1"
  }
  .remove([:kp <= :chrgp, :km <= :chrgn])
  .assign({ chrgp: :pip, chrgn: :pim })

# Mode I selection: 4C kinematic fit on K+ K− π+ π− γ, with a loose φ mass
# window applied during the fit. χ² cut 100, nominal.
modeI_selection = base_selection.dup
  .kinematic_fit([:kp, :km, :pip, :pim, :gamma]) {
    nominal
    constrain_four_momentum
    invariant_mass_of(:kp, :km).within(0.979, 1.059) # |M(K+K−) − m_φ| < 0.04
    chi2_cut 100
  }

# Mode II selection: first reconstruct η from two photons (1C Kalman fit,
# χ² < 200, at least one η), then perform a 5C kinematic fit on
# K+ K− π+ π− η (4C plus implicit η mass constraint), χ² cut 200, nominal.
modeII_selection = base_selection.dup
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 200
    neta     ">=1"
  }
  .kinematic_fit([:kp, :km, :pip, :pim, :eta]) {
    nominal
    constrain_four_momentum
    invariant_mass_of(:kp, :km).within(0.979, 1.059) # |M(K+K−) − m_φ| < 0.04
    chi2_cut 200
  }

# Algorithm for Mode I
alg_modeI = Algorithm.new("PhiEtaPrimeModeI")
alg_modeI.set_header(["PhiEtaPrimeModeIAlg/PhiEtaPrimeModeI.h"])
         .set_constant({ "ECMS" => [:double, 4.260] })
         .note(:ecms_per_energy_point,
               "ECMS varies over the 41 energy points; the correct value per " \
               "dataset is set in the BOSS jobOptions at execution time")
alg_modeI.with_decay_card(decay_card_modeI).apply(modeI_selection)
alg_modeI.execute_on(data_points + incMC_points + exMC_modeI)

# Algorithm for Mode II
alg_modeII = Algorithm.new("PhiEtaPrimeModeII")
alg_modeII.set_header(["PhiEtaPrimeModeIIAlg/PhiEtaPrimeModeII.h"])
          .set_constant({ "ECMS" => [:double, 4.260] })
          .note(:ecms_per_energy_point,
                "ECMS varies over the 41 energy points; the correct value per " \
                "dataset is set in the BOSS jobOptions at execution time")
          .note(:psi3770_search,
                "search for ψ(3770) → φη′ is performed in the ROOT analysis " \
                "using the same event selection")
          .note(:etaprime_mass_window,
                "η′ mass window is applied in the subsequent ROOT analysis, " \
                "not in BOSS")
alg_modeII.with_decay_card(decay_card_modeII).apply(modeII_selection)
alg_modeII.execute_on(data_points + incMC_points)