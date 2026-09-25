# =============================================================================
# BESIII integrated-luminosity measurement, 3.810 - 4.600 GeV
#   main channel : e+e- -> (gamma) e+e-   (large-angle Bhabha scattering)
#   cross-check  : e+e- -> gamma gamma
# Real data (18 distinct XYZ samples) + inclusive MC (4.260 GeV) + exclusive MC
# =============================================================================

### Dataset description ###
# 18 distinct real-data samples between 3.810 and 4.600 GeV; the two acquisitions
# at 4.230, 4.260 and 4.420 GeV are merged into a single sample name each.
data_samples = [
  "703_3810", "703_3872", "703_3900", "703_4009", "703_4090",
  "703_4180", "703_4190", "703_4200", "703_4210", "703_4220",
  "703_4230", "703_4260", "703_4280", "703_4360",
  "703_4420", "703_4530", "703_4575", "703_4600"
].map { |name| DatasetManager.real_data.find(name) }

# Inclusive MC at 4.260 GeV (= 500 pb-1 equivalent), used for background studies
# and for the optimisation of the selection.
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")

### Decay cards (EvtGen format) ###
# Main channel: Bhabha scattering, e+e- with an additional FSR photon.
# No resonance is produced directly, so psi(4260) is used as the top mother
# (BESIII KKMC simulation convention).
decay_card_bhabha = <<~DECAYCARD
    Decay psi(4260)
    1.0000  e+  e-  gamma     PHOTOS PHSP;
    Enddecay
    End
DECAYCARD

# Cross-check channel: e+e- -> gamma gamma.
decay_card_gg = <<~DECAYCARD
    Decay psi(4260)
    1.0000  gamma  gamma      PHSP;
    Enddecay
    End
DECAYCARD

### Exclusive MC: 1M Bhabha and 1M gamma gamma at EVERY energy point ###
# create_exclusive_mc_for runs the same signal MC once per energy point,
# automatically re-deriving the per-point configuration from the related dataset.
exMC_bhabha = DatasetManager.create_exclusive_mc_for(data_samples) do |config|
  config.sample_name   = "exmc_bhabha"     # auto-suffixed per energy point
  config.events        = 1_000_000
  config.decay_card    = decay_card_bhabha
  config.cross_section = :default
end

exMC_gg = DatasetManager.create_exclusive_mc_for(data_samples) do |config|
  config.sample_name   = "exmc_gg"         # auto-suffixed per energy point
  config.events        = 1_000_000
  config.decay_card    = decay_card_gg
  config.cross_section = :default
end

# Persist the generated MC configurations (call save on each element of the Array).
exMC_bhabha.each { |m| m.save_to_config(format: :yaml, file_path: "temp_for_test") }
exMC_gg.each     { |m| m.save_to_config(format: :yaml, file_path: "temp_for_test") }

# =============================================================================
# Main channel: e+e- -> (gamma) e+e-   (Bhabha)
# =============================================================================
alg_bhabha = Algorithm.new("BhabhaLumi")
alg_bhabha.set_header(["BhabhaLumiAlg/BhabhaLumi.h"])
          .set_constant({"ECMS" => [:double, 4.260]})   # 4.260 GeV reference values used in BOSS
          .set_alias({"std::vector<double>" => "Vdouble"})
          .note(:energy_rescaled_cuts,
                "EMC energy and momentum cuts: E_dep > sqrt(s)/4.26 x 1.55 GeV and " \
                "p > sqrt(s)/4.26 x 2.0 GeV/c. BOSS applies the 4.260 GeV reference " \
                "values (1.55 GeV and 2.0 GeV/c); the per-energy sqrt(s)/4.26 rescaling " \
                "is applied at the ROOT analysis level.")

bhabha_selection = Selection.new
  .select_track {
    cos_theta  0.8     # |cos(theta)| < 0.8
    Vz         10.0    # |Vz| < 10 cm
    Vr         1.0     # Vr < 1 cm (transverse plane)
    nChrp      "==1"   # exactly one positively charged track
    nChrn      "==1"   # exactly one negatively charged track
    nNet       "==0"   # net charge zero
  }
  # no PID: assign e+/e- purely by charge (positively charged -> e+, negatively -> e-)
  .assign({:chrgp => :ep, :chrgn => :em})
  # momentum requirement p > 2.0 GeV/c (4.260 GeV reference value)
  .remove(:ep) { condition "three_momentum_of(:ep) < 2.0" }
  .remove(:em) { condition "three_momentum_of(:em) < 2.0" }
  # nominal 4C kinematic fit of the e+e- pair to the c.m. four-momentum
  .kinematic_fit([:ep, :em]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_bhabha.with_decay_card(decay_card_bhabha).apply(bhabha_selection)

# =============================================================================
# Cross-check channel: e+e- -> gamma gamma
# =============================================================================
alg_gg = Algorithm.new("GammaGammaLumi")
alg_gg.set_header(["GammaGammaLumiAlg/GammaGammaLumi.h"])
      .set_constant({"ECMS" => [:double, 4.260]})
      .set_alias({"std::vector<double>" => "Vdouble"})
      .note(:charged_track_veto,
            "All charged tracks vetoed: after the |Vz| < 10 cm and Vr < 1 cm quality " \
            "cuts the event is required to contain zero charged tracks with net charge zero.")
      .note(:energy_rescaled_cuts,
            "The two most energetic photons are required to have E > sqrt(s)/4.26 x 1.8 GeV. " \
            "BOSS applies the 4.260 GeV reference value (1.8 GeV); the per-energy " \
            "sqrt(s)/4.26 rescaling is applied at the ROOT analysis level.")

gg_selection = Selection.new
  .select_track {
    Vz    10.0     # |Vz| < 10 cm
    Vr    1.0      # Vr < 1 cm
    nChrp "==0"    # veto all positive tracks
    nChrn "==0"    # veto all negative tracks
    nNet  "==0"    # zero net charge
  }
  .select_photon {
    tdc_emc_start     0      # TDC window 0 - 14 (700 ns units)
    tdc_emc_end       14
    angle_to_track    10.0   # angle to nearest charged track > 10 degrees
    energyThreshold_b 0.025  # barrel EMC energy threshold: 25 MeV
    energyThreshold_e 0.050  # endcap EMC energy threshold: 50 MeV
    nGam              ">=2"  # at least two photons
  }
  # nominal 4C kinematic fit of the two photons to the c.m. four-momentum
  .kinematic_fit([:gamma, :gamma]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_gg.with_decay_card(decay_card_gg).apply(gg_selection)

# =============================================================================
# Execute both algorithms on real data, inclusive MC and the exclusive MC
# =============================================================================
root_files_bhabha = alg_bhabha.execute_on(data_samples + [incMC_4260] + exMC_bhabha)
root_files_gg     = alg_gg.execute_on(data_samples + [incMC_4260] + exMC_gg)