# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset preparation ###
# BOSS 7.1.3 R-scan samples: 20 energy points from 2.05 to 3.08 GeV
scan_energies = %w[
  2050 2100 2150 2175 2200 2232 2309 2386 2396 2500
  2644 2646 2700 2800 2900 2950 2981 3000 3020 3080
]
data_points  = scan_energies.map { |e| DatasetManager.real_data.find("713_Rscan_#{e}") }     # real data at each scan point
incmc_points = scan_energies.map { |e| DatasetManager.inclusive_mc.find("713_Rscan_#{e}") } # inclusive MC at each scan point

# ConExc decay card: e+e- -> phi eta' with ISR modelled by ConExc
# (continuum / Born-cross-section scan; no KKMC psi(4260) top mother, no Particle vpho line for multi-energy)
decay_card_phi_etap = <<~DECAYCARD
    Decay vpho
    1.0 ConExc 8;
    Enddecay

    Decay phi
    1.0 K+ K- VSS;
    Enddecay

    Decay eta_prime
    1.0 pi+ pi- gamma PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive MC per energy point, sharing the same ConExc card across the scan
exMCs = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_phi_etap"
  config.events        = 100_000
  config.decay_card    = decay_card_phi_etap
  config.cross_section = :default
end
exMCs.each { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }

### Event selection (BOSS) ###
alg_name = "PhiEtaPrime"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 2.561]})   # nominal scan energy; overridden per point via ConExc vpho injection
   .set_alias({"std::vector<double>" => "Vdouble"})

# Fully reconstructed 4-track topology: e+e- -> phi eta' -> K+K- pi+pi- gamma
event_selection = Selection.new
  .select_track {
    cos_theta 0.93    # |cos(theta)| < 0.93
    Vz        10.0    # |Vz| < 10 cm
    Vr        1.0     # Vr < 1 cm
    nChrp     "==2"   # exactly 2 positive tracks
    nChrn     "==2"   # exactly 2 negative tracks
    nNet      "==0"   # net charge zero (4 charged tracks total)
  }
  .select_photon {
    energyThreshold_b 0.025   # 25 MeV barrel energy threshold
    energyThreshold_e 0.050   # 50 MeV endcap energy threshold
    tdc_emc_start     0       # EMC TDC start
    tdc_emc_end       14      # EMC TDC end
    angle_to_track    10.0    # at least 10 deg from any charged track
    nGam              ">=1"   # at least one photon (from eta' -> pi+pi-gamma)
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]   # separate pi from K/p (both charges)
    identify :kaon, against: [:pion, :proton]   # separate K from pi/p (both charges)
    npip "==1"
    npim "==1"
    nkp  "==1"
    nkm  "==1"
  }
  .kinematic_fit([:kp, :km, :pip, :pim, :gamma]) {
    nominal                  # nominal fit: corrected four-momenta are the ones stored
    constrain_four_momentum  # 4C energy-momentum constraint
    chi2_cut 200             # loose chi2 cut; tight cut applied in ROOT
  }

# BOSS-side procedures consciously left out of the encoded selection, preserved for downstream use
alg.note(:partial_reconstruction, "the 3-track topology with one missing kaon, reconstructed by a 1C kinematic fit, is a separate analysis and is not encoded in this fully reconstructed 4-track selection")
   .note(:background_veto, "initial-state-radiation photons with energy below 70 MeV are suppressed to reduce the ISR / continuum background; this photon veto is applied outside this BOSS spec")

alg.with_decay_card(decay_card_phi_etap).apply(event_selection)
root_files = alg.execute_on(data_points + incmc_points + exMCs)