# ============================================================
# Dataset preparation — e+e- -> eta pi+ pi-  (R-scan, 2.000-3.080 GeV)
# ============================================================

# 19 R-scan real-data points, BOSS release 713 (sample name = <BOSS>_<Ecms MeV>)
scan_samples = %w[
  713_2000 713_2050 713_2100 713_2125 713_2150 713_2175 713_2200 713_2232
  713_2309 713_2386 713_2396 713_2644 713_2646 713_2900 713_2950 713_2981
  713_3000 713_3020 713_3080
]
scan_data = scan_samples.map { |s| DatasetManager.real_data.find(s) }

# Inclusive MC samples at 2.125, 2.396 and 2.900 GeV
incMC_points = %w[713_2125 713_2396 713_2900].map { |s| DatasetManager.inclusive_mc.find(s) }

# ConExc decay card — continuum production of eta pi+ pi- with ISR modelling.
# ConExc mode 37 = eta pi+ pi-; the final state is generated through the vhdr
# decay to pi+ pi- eta (PHSP), eta -> gamma gamma. The literal token "ConExc"
# makes the DSL switch to the no-KKMC template and inject Particle vpho <ECMS>
# per energy point, so Particle vpho must NOT be written here.
decay_card_signal = <<~DECAYCARD
    Decay vpho
    1.000000 ConExc 37;
    Enddecay

    Decay vhdr
    1.000000 pi+ pi- eta PHSP;
    Enddecay

    Decay eta
    1.000000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive ConExc MC for every scan point (same card/cross-section,
# one signal sample per related real-data point)
exMCs = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name   = "exmc_etapipi_conexc"
  config.events        = 100_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

# ============================================================
# Event selection (BOSS)
# ============================================================
alg_name = "EtaPiPi"
eta_pipi_alg = Algorithm.new(alg_name)
eta_pipi_alg
  .set_header(["#{alg_name}Alg/#{alg_name}.h"])
  .set_constant({ "ECMS" => [:double, 2.900] })   # representative scan energy
  .set_alias({ "std::vector<double>" => "Vdouble" })

# --- BOSS-side procedures not expressible in the DSL are preserved as notes ---
eta_pipi_alg
  .note(:pid_correction_method,
        "charged pions with E/p >= 0.8 are removed as Bhabha background")
  .note(:background_veto,
        "eta mass window [0.523, 0.573] GeV/c^2 with sidebands [0.478, 0.503] and " \
        "[0.593, 0.618] GeV/c^2 applied at ROOT level for sideband background estimation")
  .note(:helicity_cut,
        "photon helicity-frame cut |cos(theta_gamma)| < 0.95 applied on the " \
        "kinematic-fit-corrected four-momenta from the nominal 4C fit (ROOT level)")
  .note(:conexc_mode_range,
        "ConExc mode 37 is tabulated for sqrt(s) = 1.025-2.975 GeV; the 3.000, 3.020 " \
        "and 3.080 GeV points lie above the range and need generator extrapolation")
  .note(:partial_wave_analysis,
        "a partial-wave analysis extracts the intermediate contributions (ROOT level)")

# Full selection chain
event_selection = Selection.new
  .select_track {                 # charged track selection
      cos_theta 0.93              # |cos(theta)| < 0.93
      Vz        10.0              # |Vz| < 10 cm
      Vr        1.0               # Vr < 1 cm
      nChrp     "==1"             # exactly one positive track
      nChrn     "==1"             # exactly one negative track
      nNet      "==0"             # net charge zero
  }
  .select_photon {                # photon selection
      tdc_emc_start     0         # TDC start
      tdc_emc_end       14        # TDC end
      energyThreshold_b 0.025     # 25 MeV in EMC barrel
      energyThreshold_e 0.050     # 50 MeV in EMC endcap
      angle_to_track    10.0      # more than 10 deg from any charged track
      nGam              ">=2"     # at least two photons
  }
  .pid(method: :probability) {    # particle identification
      prob_cut 0.001
      identify :pion, against: [:kaon, :proton]   # pi+ and pi- vs K, p
      npip "==1"                  # exactly one pi+
      npim "==1"                  # exactly one pi-
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {       # 1C fit: reconstruct eta -> gamma gamma
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
      chi2_cut 25
      neta ">=1"                  # at least one eta candidate
  }
  # 4C kinematic fit to pi+ pi- eta, correcting the four-momenta to the CMS energy.
  # The photon helicity cut and the eta mass window are applied after this final
  # fit and therefore belong to the ROOT stage (see notes above).
  .kinematic_fit([:pip, :pim, :eta]) {
      nominal                     # this is the nominal fit
      vertex_fit([0, 1])          # common vertex for pip(0) and pim(1)
      constrain_four_momentum     # 4-momentum conservation vs CMS energy
      chi2_cut 100
  }

# Generate the algorithm and run over data, inclusive MC and signal MC
eta_pipi_alg.with_decay_card(decay_card_signal).apply(event_selection)
root_files = eta_pipi_alg.execute_on(scan_data + incMC_points + exMCs)