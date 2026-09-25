# BESIII DSL: e+e- -> eta' pi+ pi- Born cross sections at 19 energies
# Paper: 2012.07360v3
# Energy range: 2.00 - 3.08 GeV, ConExc generator, two eta' decay modes

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# ============================================================
# 19 R-scan datasets at BOSS 713
# Energies: 2.000, 2.050, 2.100, 2.125, 2.150, 2.175, 2.200,
#           2.232, 2.309, 2.386, 2.396, 2.500, 2.644, 2.646,
#           2.700, 2.800, 2.900, 2.950, 3.080 GeV
# (Note: exact match to paper's energy list; some dataset names
#  use slightly different energy labels than exact c.m. energy)
# ============================================================

scan_datasets = [
  DatasetManager.real_data.find("713_Rscan_2000"),
  DatasetManager.real_data.find("713_Rscan_2050"),
  DatasetManager.real_data.find("713_Rscan_2100"),
  DatasetManager.real_data.find("713_Rscan_2125"),
  DatasetManager.real_data.find("713_Rscan_2150"),
  DatasetManager.real_data.find("713_Rscan_2175"),
  DatasetManager.real_data.find("713_Rscan_2200"),
  DatasetManager.real_data.find("713_Rscan_2232"),
  DatasetManager.real_data.find("713_Rscan_2309"),
  DatasetManager.real_data.find("713_Rscan_2386"),
  DatasetManager.real_data.find("713_Rscan_2396"),
  DatasetManager.real_data.find("713_Rscan_2500"),
  DatasetManager.real_data.find("713_Rscan_2644"),
  DatasetManager.real_data.find("713_Rscan_2646"),
  DatasetManager.real_data.find("713_Rscan_2700"),
  DatasetManager.real_data.find("713_Rscan_2800"),
  DatasetManager.real_data.find("713_Rscan_2900"),
  DatasetManager.real_data.find("713_Rscan_2950"),
  DatasetManager.real_data.find("713_Rscan_3080"),
]

scan_incMCs = scan_datasets.map do |ds|
  name = ds.sample_name
  DatasetManager.inclusive_mc.find(name)
rescue => e
  nil
end.compact

# ============================================================
# ConExc decay card
# Mode 40 = eta' pi+ pi- at sqrt(s) 1.58-3.42 GeV
# Multi-energy auto-injection: omit "Particle vpho" line
# ============================================================

conexc_card = <<~DECAYCARD
  Decay vpho
    1 ConExc 40;
  Enddecay
  Decay vhdr
    1 eta_prime pi+ pi- PHSP;
  Enddecay
  Decay eta_prime
    1 eta pi+ pi- PHSP;
  Enddecay
  Decay eta
    1 gamma gamma PHSP;
  Enddecay
  End
DECAYCARD

# ============================================================
# Exclusive MC created across the 19 energy scan points
# ============================================================

sig_mc_scan = DatasetManager.create_exclusive_mc_for(scan_datasets) do |config|
  config.sample_name   = "sig_etap_pipi_scan"
  config.events        = 200_000
  config.decay_card    = conexc_card
  config.cross_section = :default
end

# ============================================================
# Algorithm Mode I: eta' -> eta pi+ pi-, eta -> gamma gamma
# 4C kinematic fit: 2 charged (pi+ pi-) + 4 photons (eta' pi+ pi-)
# chi2_cut 100
# Competing vetoes: 2(pi+pi-) and 2(gamma pi+ pi-)
# ============================================================

alg_modeI = Algorithm.new("EtapPipiModeI", version: '00-00-01')
alg_modeI.set_header(["EtapPipiModeIAlg/EtapPipiModeI.h"])
          .set_constant({ "ECMS" => [:double, 2.125] })   # nominal, overridden per energy

modeI_selection = Selection.new
  .select_track do
    nChrp ">=2"       # 2 pi+ (eta' -> eta pi+ pi- has one pi+; plus the direct pi+)
    nChrn ">=2"       # 2 pi-
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
    nGam              ">=4"     # 2 from eta -> gamma gamma + 2 from eta' decay photons
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
  end
  # Reconstruct eta from gamma gamma
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=1"
  end
  # Nominal 4C kinematic fit: pi+, pi-, pi+, pi-, eta (4 charged + 1 composite)
  .kinematic_fit([:pip, :pim, :pip, :pim, :eta]) do
    constrain_four_momentum
    chi2_cut 100
    nominal
  end
  # Competing veto hypothesis 1: 2(pi+ pi-) -> 4 pions (no eta)
  .kinematic_fit([:pip, :pim, :pip, :pim]) do
    constrain_four_momentum
    # non-nominal, no chi2 cut; chi2 compared in ROOT
  end
  # Competing veto hypothesis 2: 2(gamma pi+ pi-) -> 2 photons + 4 pions (no eta)
  .kinematic_fit([:gamma, :gamma, :pip, :pim, :pip, :pim]) do
    constrain_four_momentum
    # non-nominal, no chi2 cut; chi2 compared in ROOT
  end

alg_modeI.note(:etap_mass_window, "eta' mass (0.94, 0.97) GeV/c^2 applied in ROOT")
alg_modeI.note(:eta_mass_window, "eta mass (0.518, 0.578) GeV/c^2 applied in ROOT; sideband method for background")
alg_modeI.note(:pipi_mass_cut, "M(pi+pi-) from eta' decay > 0.5 GeV/c^2 to suppress backgrounds")
alg_modeI.note(:two_pipi_veto, "Reject if chi2_2pipi < chi2_nominal; 2(pi+pi-) hypothesis veto")
alg_modeI.note(:two_gamma_pipi_veto, "Reject if chi2_2gamma_pipi < chi2_nominal; 2(gamma pi+pi-) hypothesis veto")
alg_modeI.note(:isr_correction, "ISR correction factor (1+delta) applied in ROOT to obtain Born cross section")
alg_modeI.note(:vp_correction, "Vacuum polarization correction applied in ROOT")

alg_modeI.with_decay_card(conexc_card).apply(modeI_selection)
alg_modeI.execute_on(scan_datasets + scan_incMCs + sig_mc_scan)

# ============================================================
# Algorithm Mode II: eta' -> gamma pi+ pi-
# 4C kinematic fit: 2 charged (pi+ pi-) + 3 photons (radiative gamma + 2 from eta')
# chi2_cut 50
# Competing vetoes: 2(pi+pi-) and 2(gamma pi+ pi-)
# ============================================================

alg_modeII = Algorithm.new("EtapPipiModeII", version: '00-00-01')
alg_modeII.set_header(["EtapPipiModeIIAlg/EtapPipiModeII.h"])
           .set_constant({ "ECMS" => [:double, 2.125] })   # nominal, overridden per energy

modeII_selection = Selection.new
  .select_track do
    nChrp ">=2"       # pi+, pi+ (from eta' decay + direct)
    nChrn ">=2"       # pi-, pi-
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
    nGam              ">=3"     # eta' -> gamma pi+ pi- gives 1 photon, plus 2 direct
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
  end
  # Nominal 4C kinematic fit: 4 charged pions + 3 photons
  .kinematic_fit([:pip, :pim, :pip, :pim, :gamma, :gamma, :gamma]) do
    constrain_four_momentum
    chi2_cut 50
    nominal
  end
  # Competing veto hypothesis 1: 2(pi+ pi-) -> 4 pions (minimal photons)
  .kinematic_fit([:pip, :pim, :pip, :pim]) do
    constrain_four_momentum
    # non-nominal, no chi2 cut
  end
  # Competing veto hypothesis 2: 2(gamma pi+ pi-) -> 2 photons + 4 pions
  .kinematic_fit([:gamma, :gamma, :pip, :pim, :pip, :pim]) do
    constrain_four_momentum
    # non-nominal, no chi2 cut
  end

alg_modeII.note(:etap_mass_window, "eta' mass (0.94, 0.98) GeV/c^2 applied in ROOT")
alg_modeII.note(:gamma_energy_cut, "Photon energy > 0.1 GeV (to suppress ISR photons); applied in ROOT")
alg_modeII.note(:pipi_mass_cut, "M(pi+pi-) from eta' decay for rho dominance check; applied in ROOT")
alg_modeII.note(:two_pipi_veto, "Reject if chi2_2pipi < chi2_nominal; 2(pi+pi-) hypothesis veto")
alg_modeII.note(:two_gamma_pipi_veto, "Reject if chi2_2gamma_pipi < chi2_nominal; 2(gamma pi+pi-) hypothesis veto")
alg_modeII.note(:isr_correction, "ISR correction factor (1+delta) applied in ROOT to obtain Born cross section")
alg_modeII.note(:vp_correction, "Vacuum polarization correction applied in ROOT")

alg_modeII.with_decay_card(conexc_card).apply(modeII_selection)
alg_modeII.execute_on(scan_datasets + scan_incMCs + sig_mc_scan)