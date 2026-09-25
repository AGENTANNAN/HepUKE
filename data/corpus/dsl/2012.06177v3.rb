# BESIII DSL: Search for X(2370) and eta_c in J/psi -> gamma eta eta eta'
# Paper: 2012.06177v3
# J/psi at sqrt(s) = 3.097 GeV, two eta' decay modes

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# ============================================================
# Decay cards
# ============================================================

# Mode I: J/psi -> gamma eta eta eta', eta' -> gamma pi+ pi-, eta -> gamma gamma
decay_card_modeI = <<~DECAYCARD
  Decay J/psi
    1 gamma eta eta eta_prime PHSP;
  Enddecay
  Decay eta_prime
    1 gamma pi+ pi- PHSP;
  Enddecay
  Decay eta
    1 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# Mode II: J/psi -> gamma eta eta eta', eta' -> pi+ pi- eta, eta -> gamma gamma
decay_card_modeII = <<~DECAYCARD
  Decay J/psi
    1 gamma eta eta eta_prime PHSP;
  Enddecay
  Decay eta_prime
    1 pi+ pi- eta PHSP;
  Enddecay
  Decay eta
    1 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# Exclusive MC samples
sig_mc_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_eta3etap_modeI"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

sig_mc_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_eta3etap_modeII"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

# ============================================================
# Algorithm Mode I: eta' -> gamma pi+ pi-
# 6C kinematic fit (4C + 2x1C eta mass)
# >=6 photons, 2 charged tracks, no PID
# ============================================================

alg_modeI = Algorithm.new("Eta3EtapModeI", version: '00-00-01')
alg_modeI.set_header(["Eta3EtapModeIAlg/Eta3EtapModeI.h"])
          .set_constant({ "ECMS" => [:double, 3.097] })

modeI_selection = Selection.new
  .select_track do
    nChrp ">=1"       # pi+
    nChrn ">=1"       # pi-
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=6"     # at least 6 for 3 eta -> gamma gamma; plus 2 (radiative + eta' decay gamma)
  end
  # Reconstruct three etas from gamma gamma pairs
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=3"         # need 3 etas
  end
  # 6C kinematic fit: 4C + 2 eta mass constraints
  # Participants: radiative gamma (1), 3 etas, eta'-decay gamma (1), pi+, pi-
  .kinematic_fit([:gamma, :eta, :eta, :eta, :gamma, :pip, :pim]) do
    constrain_four_momentum
    invariant_mass_of(:eta).constrain_to_nominal_mass_of(:eta)
    invariant_mass_of(:eta).constrain_to_nominal_mass_of(:eta)
    chi2_cut 30
    nominal
  end

alg_modeI.note(:pi0_veto, "pi0 veto: reject if |M(gamma gamma) - M_pi0| < 0.015 GeV/c^2 for any gamma pair; applied in ROOT")
alg_modeI.note(:eta_veto, "eta veto: reject extra eta candidates beyond the 3 used; applied in ROOT")
alg_modeI.note(:pipi_mass_cut, "M(pi+pi-) > 0.5 GeV/c^2 applied in ROOT to suppress backgrounds from rho-like structures")
alg_modeI.note(:eta_mass_window, "eta mass sidebands used for background estimation via ROOT")
alg_modeI.note(:etap_mass_region, "eta' mass region for signal extraction via ROOT fitting")
alg_modeI.note(:ml_reweighting, "ML reweighting for simulation-data agreement; not expressible in BOSS DSL")
alg_modeI.note(:photon_identification, "Radiative photon identified as the one not from eta or eta' decay; combinatorics handled by kinematic fit")

alg_modeI.with_decay_card(decay_card_modeI).apply(modeI_selection)
alg_modeI.execute_on([jpsi_data, jpsi_incMC, sig_mc_modeI])

# ============================================================
# Algorithm Mode II: eta' -> pi+ pi- eta, eta -> gamma gamma
# 7C kinematic fit (4C + 3x1C eta mass)
# >=7 photons, 2 charged tracks, no PID
# ============================================================

alg_modeII = Algorithm.new("Eta3EtapModeII", version: '00-00-01')
alg_modeII.set_header(["Eta3EtapModeIIAlg/Eta3EtapModeII.h"])
           .set_constant({ "ECMS" => [:double, 3.097] })

modeII_selection = Selection.new
  .select_track do
    nChrp ">=1"       # pi+
    nChrn ">=1"       # pi-
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=7"     # 6 for 3 eta -> gamma gamma + 1 radiative
  end
  # Reconstruct three etas from gamma gamma pairs
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=3"         # need 3 etas
  end
  # 7C kinematic fit: 4C + 3 eta mass constraints
  # Participants: radiative gamma (1), 3 etas, pi+, pi-
  .kinematic_fit([:gamma, :eta, :eta, :eta, :pip, :pim]) do
    constrain_four_momentum
    invariant_mass_of(:eta).constrain_to_nominal_mass_of(:eta)
    invariant_mass_of(:eta).constrain_to_nominal_mass_of(:eta)
    invariant_mass_of(:eta).constrain_to_nominal_mass_of(:eta)
    chi2_cut 50
    nominal
  end

alg_modeII.note(:pi0_veto, "pi0 veto: reject if |M(gamma gamma) - M_pi0| < 0.015 GeV/c^2 for any gamma pair; applied in ROOT")
alg_modeII.note(:eta_veto, "eta veto: reject extra eta candidates beyond the 3 used; applied in ROOT")
alg_modeII.note(:etap_mass_region, "eta' mass region for signal extraction via ROOT fitting")
alg_modeII.note(:ml_reweighting, "ML reweighting for simulation-data agreement; not expressible in BOSS DSL")
alg_modeII.note(:photon_identification, "Radiative photon identified as the one not from eta or eta' decay; combinatorics handled by kinematic fit")

alg_modeII.with_decay_card(decay_card_modeII).apply(modeII_selection)
alg_modeII.execute_on([jpsi_data, jpsi_incMC, sig_mc_modeII])