# Core DSL classes and dependencies are loaded automatically at execution

### Dataset description ###
# Real data and inclusive MC over the whole scan window, 2.232-4.590 GeV (131 points)
scan_data  = DatasetManager.real_data.where(cms_energy: { value: 2232.0..4590.0 })
scan_incMC = DatasetManager.inclusive_mc.where(cms_energy: { value: 2232.0..4590.0 })

# ------------------------------------------------------------------ decay cards
# Bhabha signal e+e- -> (gamma) e+e- : top mother psi(4260) (KKMC convention),
# PHOTOS supplies the radiative (ISR/FSR) photons.
decay_card_bhabha = <<~DECAYCARD
    Decay psi(4260)
    1.0000 e+ e- gamma PHOTOS PHSP;
    Enddecay
    End
DECAYCARD

# Diphoton signal e+e- -> gamma gamma
decay_card_diphoton = <<~DECAYCARD
    Decay psi(4260)
    1.0000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

# tau+tau- QED background
decay_card_tau = <<~DECAYCARD
    Decay psi(4260)
    1.0000 tau+ tau- PHSP;
    Enddecay
    End
DECAYCARD

# Generic continuum-hadron background (PHSP because the generator is unspecified)
decay_card_hadron = <<~DECAYCARD
    Decay psi(4260)
    1.0000 pi+ pi- pi0 PHSP;
    Enddecay
    End
DECAYCARD

# ------------------------------------------ per-energy-point exclusive MC (131)
exMC_bhabha = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name   = "exmc_qed_bhabha"
  config.events        = 200_000
  config.decay_card    = decay_card_bhabha
  config.cross_section = :default
end

exMC_diphoton = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name   = "exmc_qed_diphoton"
  config.events        = 200_000
  config.decay_card    = decay_card_diphoton
  config.cross_section = :default
end

exMC_tau = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name   = "exmc_tautau_bkg"
  config.events        = 100_000
  config.decay_card    = decay_card_tau
  config.cross_section = :default
end

exMC_hadron = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name   = "exmc_continuum_hadron_bkg"
  config.events        = 100_000
  config.decay_card    = decay_card_hadron
  config.cross_section = :default
end

### Event selection (BOSS) ###
# ------------------------------------------------------------ Bhabha chain ---
alg_name_bhabha = "QEDBhabha"
alg_bhabha = Algorithm.new(alg_name_bhabha)
alg_bhabha.set_header(["#{alg_name_bhabha}Alg/#{alg_name_bhabha}.h"])
          .set_constant({"ECMS" => [:double, 3.773]})   # placeholder; true beam energy is per point
          .set_alias({"std::vector<double>" => "Vdouble"})
          .note(:ecms_per_point, "the measurement spans 131 c.m. energy points (2.232-4.590 GeV); the ECMS constant is only a code-generation placeholder, the actual per-point beam energy is taken from each dataset at execution time")

bhabha_selection = Selection.new
bhabha_selection.select_track {          # two charged tracks: exactly one + and one -
                    cos_theta 0.8        # |cos(theta)| < 0.8
                    Vz 10.0              # |Vz| < 10 cm
                    Vr 1.0               # Vr < 1 cm
                    nChrp "==1"          # exactly one positively charged track
                    nChrn "==1"          # exactly one negatively charged track
                    nNet  "==0"          # net charge zero
                  }
                 # per-track EMC energy over momentum > 0.65 (no PID is used)
                 .remove(:chrgp) { condition "ep_ratio_of(:chrgp) <= 0.65" }
                 .remove(:chrgn) { condition "ep_ratio_of(:chrgn) <= 0.65" }
                 .assign({:chrgp => :ep, :chrgn => :em})   # positive -> e+, negative -> e-
                 .kinematic_fit([:ep, :em]) {              # nominal 4C fit to the e+e- four-momentum
                    nominal
                    constrain_four_momentum
                    chi2_cut 200
                 }

# ---------------------------------------------------------- Diphoton chain ---
alg_name_diphoton = "QEDDiphoton"
alg_diphoton = Algorithm.new(alg_name_diphoton)
alg_diphoton.set_header(["#{alg_name_diphoton}Alg/#{alg_name_diphoton}.h"])
           .set_constant({"ECMS" => [:double, 3.773]})   # placeholder; true beam energy is per point
           .set_alias({"std::vector<double>" => "Vdouble"})
           .note(:ecms_per_point, "the measurement spans 131 c.m. energy points (2.232-4.590 GeV); the ECMS constant is only a code-generation placeholder, the actual per-point beam energy is taken from each dataset at execution time")

diphoton_selection = Selection.new
diphoton_selection.select_track {        # veto charged tracks
                    cos_theta 0.8
                    Vz 10.0
                    Vr 1.0
                    nChrp "==0"          # zero positively charged tracks
                    nChrn "==0"          # zero negatively charged tracks
                    nNet  "==0"          # net charge zero
                  }
                  .select_photon {       # at least two good photons
                    tdc_emc_start 0      # EMC timing window 0-700 ns
                    tdc_emc_end   14
                    energyThreshold_b 0.025   # > 25 MeV in the barrel
                    energyThreshold_e 0.050   # > 50 MeV in the endcap
                    angle_to_track 10.0       # at least 10 degrees from any charged track
                    nGam ">=2"
                  }
                  .kinematic_fit([:gamma, :gamma]) {   # nominal 4C fit to the two-photon four-momentum
                    nominal
                    constrain_four_momentum
                    chi2_cut 200
                  }

# ------------------------------------------------------- generate the code ---
alg_bhabha.with_decay_card(decay_card_bhabha).apply(bhabha_selection)
alg_diphoton.with_decay_card(decay_card_diphoton).apply(diphoton_selection)

# --------------------------------------------------------------- execution ---
# Selection chains are energy independent and are applied to real data, inclusive
# MC, the shared QED backgrounds, and each chain's own signal exclusive MC.
common_datasets = scan_data + scan_incMC + exMC_tau + exMC_hadron

root_files_bhabha   = alg_bhabha.execute_on(common_datasets + exMC_bhabha)
root_files_diphoton = alg_diphoton.execute_on(common_datasets + exMC_diphoton)