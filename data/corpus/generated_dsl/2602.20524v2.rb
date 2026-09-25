# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset preparation ###
jpsi_data  = DatasetManager.real_data.find("708_3097")       # J/psi (3.097 GeV) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")    # Corresponding inclusive MC sample

# Decay card for the signal process e+e- -> J/psi -> Xi- anti-Xi+, Xi- -> Lambda pi-, Lambda -> p pi-
# (EvtGen syntax / EvtGen particle names). Generated in phase space (PHSP) as requested.
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 Xi- anti-Xi+ PHSP;
    Enddecay

    Decay Xi-
    1.0000 Lambda0 pi- PHSP;
    Enddecay

    Decay anti-Xi+
    1.0000 anti-Lambda0 pi+ PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi- PHSP;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ PHSP;
    Enddecay

    End
DECAYCARD

# Phase-space exclusive MC (2M events) for the signal chain, used to extract the normalisation factor
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_xim_xipbar_phsp"
  config.related_dataset = jpsi_data          # Associated real dataset (better simulation)
  config.events          = 2_000_000          # 2M phase-space events
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "JpsiXiXibar"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})      # CMS energy of the J/psi run
            .set_alias({"std::vector<double>" => "Vdouble"})

# BOSS-side procedures that have no DSL construct are preserved as notes
my_algorithm
  .note(:xi_delta_cut, "per-Xi candidate cut delta < 0.016 GeV/c^2 applied after the secondary-vertex "
        "fits, with delta = sqrt((m(Lambda pi-) - m_Xi-)^2 + R^2 * (m(p pi-) - m_Lambda)^2) and "
        "R = sigma_Xi- / sigma_Lambda; the combined mass-difference/variable-resolution quantity is "
        "not expressible in the current DSL and is applied in the generated selection code")
  .note(:decay_length_cut, "Xi candidates with negative decay length rejected after the "
        "secondary-vertex fit")

event_selection = Selection.new
event_selection.select_track {                 # Charged track selection
                  cos_theta 0.93               # |cos(theta)| < 0.93
                  Vz        100.0              # |Vz| < 100 cm along the beam direction
                  Vr        10.0               # |Vr| < 10 mm in the transverse plane
                  nChrp     ">=3"              # At least 3 positively charged tracks
                  nChrn     ">=3"              # At least 3 negatively charged tracks
                }
               .pid(method: :probability) {    # PID with the probability method
                  prob_cut 0.001               # PID probability > 0.001
                  identify :proton, against: [:kaon, :pion]  # p+ / p-bar vs K, pi
                  identify :pion,   against: [:kaon, :proton] # pi+ / pi- vs K, p
                  nprp ">=1"                   # At least one proton
                  nprm ">=1"                   # At least one anti-proton
                  npip ">=2"                   # At least two pi+
                  npim ">=2"                   # At least two pi-
                }
               # Momentum-dependent PID requirement: protons above 0.3 GeV/c, pions below 0.3 GeV/c
               .remove(:prp) { condition "three_momentum_of(:prp) < 0.3" }
               .remove(:prm) { condition "three_momentum_of(:prm) < 0.3" }
               .remove(:pip) { condition "three_momentum_of(:pip) > 0.3" }
               .remove(:pim) { condition "three_momentum_of(:pim) > 0.3" }
               # Reconstruct Lambda -> p pi- through a secondary-vertex fit
               .secondary_vertex_fit([:prp, :pim]) {
                  build_virtual_particle(:Lambda).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
               }
               # Reconstruct anti-Lambda -> p-bar pi+ through a secondary-vertex fit
               .secondary_vertex_fit([:prm, :pip]) {
                  build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
               }
               # Reconstruct Xi- -> Lambda pi- through a secondary-vertex fit
               .secondary_vertex_fit([:Lambda, :pim]) {
                  build_virtual_particle(:Xi_minus).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
               }
               # Reconstruct anti-Xi+ -> anti-Lambda pi+ through a secondary-vertex fit
               .secondary_vertex_fit([:Lambda_bar, :pip]) {
                  build_virtual_particle(:Xi_plus_bar).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
               }
               # 4C kinematic fit to the Xi- anti-Xi+ system
               .kinematic_fit([:Xi_minus, :Xi_plus_bar]) {
                  nominal                       # Nominal fit: corrected four-momenta from this fit are kept
                  constrain_four_momentum       # 4C energy-momentum constraint
                  chi2_cut 200                  # Loose chi2 < 200 (tight cut applied in ROOT)
               }
               # NOTE: the two-dimensional sideband in m(Lambda pi-) vs m(anti-Lambda pi+)
               # (1.274-1.306 and 1.338-1.370 GeV/c2) is used for background subtraction and is
               # applied at the ROOT level, i.e. after the nominal kinematic fit.

my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])