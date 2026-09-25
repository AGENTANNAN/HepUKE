# ==================================================================
# Born cross section of e+e- -> Xi- anti-Xi+ via single-baryon tagging
# Eight BESIII R-scan points (BOSS 713, round17 2024)
# ==================================================================

### Dataset preparation ###

# Real data at the eight scan points
data_points = [
  DatasetManager.real_data.find("713_Rscan_2644"),  # 2.644 GeV, 34.0 pb^-1
  DatasetManager.real_data.find("713_Rscan_2646"),  # 2.646 GeV, 33.7 pb^-1
  DatasetManager.real_data.find("713_Rscan_2900"),  # 2.900 GeV, 105  pb^-1
  DatasetManager.real_data.find("713_Rscan_2950"),  # 2.950 GeV, 15.9 pb^-1
  DatasetManager.real_data.find("713_Rscan_2981"),  # 2.981 GeV, 16.1 pb^-1
  DatasetManager.real_data.find("713_Rscan_3000"),  # 3.000 GeV, 15.9 pb^-1
  DatasetManager.real_data.find("713_Rscan_3020"),  # 3.020 GeV, 17.3 pb^-1
  DatasetManager.real_data.find("713_Rscan_3080")   # 3.080 GeV, 126  pb^-1
]

# Matching inclusive MC (continuum background) at each point
incMC_points = [
  DatasetManager.inclusive_mc.find("713_Rscan_2644"),
  DatasetManager.inclusive_mc.find("713_Rscan_2646"),
  DatasetManager.inclusive_mc.find("713_Rscan_2900"),
  DatasetManager.inclusive_mc.find("713_Rscan_2950"),
  DatasetManager.inclusive_mc.find("713_Rscan_2981"),
  DatasetManager.inclusive_mc.find("713_Rscan_3000"),
  DatasetManager.inclusive_mc.find("713_Rscan_3020"),
  DatasetManager.inclusive_mc.find("713_Rscan_3080")
]

# ConExc decay card for the continuum Born process (ISR modelled to O(alpha^2)).
# Detected by the literal token `ConExc`: no KKMC top mother and no `Particle vpho`
# line (the DSL injects `Particle vpho <ECMS> 0.0` per energy point).
decay_card_signal = <<~DECAYCARD
    Decay vpho
    1.0 ConExc vhdr Xi- anti-Xi+;
    Enddecay

    Decay Xi-
    1.0 Lambda0 pi- PHSP;
    Enddecay

    Decay anti-Xi+
    1.0 anti-Lambda0 pi+ PHSP;
    Enddecay

    Decay Lambda0
    1.0 p+ pi- PHSP;
    Enddecay

    Decay anti-Lambda0
    1.0 anti-p- pi+ PHSP;
    Enddecay

    End
DECAYCARD

# One 100k-event exclusive signal MC per scan point (same card / xs / events,
# differing only by related_dataset)
exMCs_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_ximinus_xibarplus"
  config.events        = 100_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###

alg_name = "XiSingleTag"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.000]})   # representative point; beam energy set per dataset
   .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
  .select_track {                                # charged track selection
    cos_theta 0.93                               # |cos(theta)| < 0.93
    Vz 10.0                                      # |Vz| < 10 cm
    Vr 1.0                                       # Vr < 1 cm
    nChrp ">=2"                                  # at least 2 positive tracks
    nChrn ">=3"                                  # at least 3 negative tracks
  }
  .select_photon {                               # photon selection
    tdc_emc_start 0                              # EMC TDC window 0-14
    tdc_emc_end 14
    angle_to_track 10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
  }
  .pid(method: :probability) {                   # probability-method PID
    prob_cut 0.001                               # prob > 0.001
    identify :proton, against: [:kaon, :pion]    # p+ separated from K / pi
    identify :pion,   against: [:kaon, :proton]  # pi separated from K / p
    nprp ">=1"                                   # >= 1 proton (positive charge)
    npim ">=2"                                   # >= 2 pi-
  }
  .secondary_vertex_fit([:prp, :pim]) {           # Lambda -> p pi-
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .secondary_vertex_fit([:Lambda, :pim]) {        # Xi- -> Lambda pi-
    build_virtual_particle(:Xi_m).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:Xi_m, :pim]) {                 # nominal 4C fit on the tagged Xi- side
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg.note(:isr_correction,
         "ISR correction iterated to the 1% level, evaluated from the ConExc generator log at each scan point")
   .note(:recoil_mass_definition,
         "anti-Xi+ is not reconstructed; the single-tag recoil mass M_recoil = sqrt((E_cm - E_(Lambda pi))^2 - |p_cm - p_(Lambda pi)|^2) is formed against the tagged Xi- and used for the unbinned maximum-likelihood yield fit and the 90% CL upper limits at ROOT level")
   .note(:background_veto,
         "double-counting of the single-tag candidates (~19%) corrected at ROOT level")
   .with_decay_card(decay_card_signal)
   .apply(event_selection)

# Run on all eight data points, their inclusive MC and the signal MC
root_files = alg.execute_on(data_points + incMC_points + exMCs_signal)