# Paper: 2211.11253v2
# Title: Normalized Differential Cross Sections of Inclusive pi0 and KS0 Production
# Energy: 2.2324, 2.4000, 2.8000, 3.0500, 3.4000, 3.6710 GeV (6 energy points)
# This is an INCLUSIVE hadron production measurement, not exclusive channel study
# Uses luarlw generator for inclusive hadronic MC; no EvtGen decay card for signal

# Note: This is a fundamentally different analysis type from exclusive channel studies.
# The DSL is designed for exclusive analyses with decay cards and kinematic fits.
# We provide the best-effort translation focusing on dataset preparation and
# the charged track / photon / pi0 / KS0 reconstruction selection.

### Dataset preparation ###
# 6 energy points, using closest available BESIII data samples
data_2232 = DatasetManager.load_real_data.find("713_Rscan_2232")
data_2400 = DatasetManager.load_real_data.find("713_Rscan_2396")  # closest to 2400 MeV
data_2800 = DatasetManager.load_real_data.find("713_Rscan_2800")
data_3080 = DatasetManager.load_real_data.find("713_Rscan_3080")  # closest to 3050 MeV

# For 3400 and 3671 MeV points, use available samples
data_3650 = DatasetManager.load_real_data.find("709_3650")  # closest to 3400 MeV
data_3686 = DatasetManager.load_real_data.find("709_3686")  # closest to 3671 MeV

incMC_2232 = DatasetManager.load_inclusive_mc.find("713_Rscan_2232")
incMC_2396 = DatasetManager.load_inclusive_mc.find("713_Rscan_2396")
incMC_2800 = DatasetManager.load_inclusive_mc.find("713_Rscan_2800")
incMC_3080 = DatasetManager.load_inclusive_mc.find("713_Rscan_3080")
incMC_3650 = DatasetManager.load_inclusive_mc.find("709_3650")
incMC_3686 = DatasetManager.load_inclusive_mc.find("709_3686")

all_data  = [data_2232, data_2400, data_2800, data_3080, data_3650, data_3686]
all_incMC = [incMC_2232, incMC_2396, incMC_2800, incMC_3080, incMC_3650, incMC_3686]

# No exclusive signal MC generated via decay card - this analysis uses
# luarlw inclusive hadronic generator. Signal processes e+e- -> pi0/KS0 + X
# are embedded in the inclusive MC.

### Event selection ###
alg = Algorithm.new("InclPi0KS0")
alg.set_header(["InclPi0KS0Alg/InclPi0KS0.h"])
   .set_constant({"ECMS" => [:double, 2.800]})
   .set_alias({"std::vector<double>" => "Vdouble"})

# Note: Because ECMS varies across 6 energy points, the constant is set
# to an approximate central value. The actual ECMS used per energy point
# is handled by BOSS at runtime.

# Charged track selection for inclusive hadronic events
# Paper: |cos(theta)| < 0.93
# At least 2 good charged tracks (events with >3 tracks kept directly)
event_selection = Selection.new
event_selection.select_track {
                  cos_theta 0.93
                  nTot ">=2"
                }
                # Photon selection for pi0 reconstruction
                # Energy > 25 MeV barrel (|cos(theta)|<0.80)
                # Energy > 50 MeV endcap (0.86<|cos(theta)|<0.92)
                # Angle to nearest charged track > 10 deg
                # EMC time [0, 700] ns
                .select_photon {
                  tdc_emc_start 0
                  tdc_emc_end   700
                  angle_to_track 10.0
                  energyThreshold_b 0.025
                  energyThreshold_e 0.050
                }
                # Build pi0 from photon pairs via Kalman fit
                .kalman_kinematic_fit([:gamma, :gamma]) {
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                  chi2_cut 200
                  npi0 ">=1"
                }
                # Build KS0 from secondary vertex fit of pi+ pi-
                .secondary_vertex_fit([:pip, :pim]) {
                  build_virtual_particle(:K_S0).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
                }

# No kinematic fit - this is an inclusive measurement.
# Cross sections are normalized by total hadronic cross section.
# Signal yields extracted by fitting M(gammagamma) and M(pi+pi-) spectra.

alg.note(:luarlw_inclusive_mc,
  "Uses luarlw inclusive hadronic MC generator (not exclusive EvtGen decay cards). Signal e+e- -> pi0/KS0 + X embedded in inclusive MC.")

alg.note(:qed_background,
  "QED backgrounds (e+e- -> e+e-, mu+mu-, gammagamma) generated with babayaga3.5. tau+tau- at 3.671 GeV via KKMC.")

alg.note(:no_kinematic_fit,
  "No kinematic fit used. Signal/background separated by fitting M(gammagamma) and M(pi+pi-) spectra in momentum bins.")

alg.note(:pi0_angular_cut,
  "pi0 candidates: |cos(theta_gamma)| < 0.8 for p_pi0 < 0.3 GeV/c; |cos(theta_gamma)| < 0.95 for p_pi0 > 0.3 GeV/c. Applied in ROOT.")

alg.note(:ks0_detector_cuts,
  "KS0 tracks: |cos(theta)|<0.93, Vz<30cm, Vxy<10cm. Decay length > 2*sigma. No PID applied on KS0 daughter tracks.")

alg.note(:energy_sample_notes,
  "6 energy points: 2.2324(713/2232), 2.4000(~2396), 2.8000(713/2800), 3.0500(~3080), 3.4000(~3650), 3.6710(~3686). Actual samples may differ.")

event_selection  # Return the selection chain

alg.apply(event_selection)

all_datasets = all_data + all_incMC
root_files = alg.execute_on(all_datasets)