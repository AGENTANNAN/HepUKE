# Analysis: Precise measurement of matter-antimatter asymmetry with entangled
# hyperon-antihyperon pairs at BESIII.
# Process: e+ e- -> J/psi -> Xi- Xibar+, Xi- -> Lambda pi-, Lambda -> p pi-.

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Signal exclusive decay card
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 Xi-  anti-Xi+                    PHSP;
    Enddecay

    Decay Xi-
    1.0000 Lambda0 pi-                      PHSP;
    Enddecay

    Decay anti-Xi+
    1.0000 anti-Lambda0 pi+                 PHSP;
    Enddecay

    Decay Lambda0
    1.0000 p+ pi-                           PHSP;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+                      PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive phase-space uniform MC sample used to determine the normalisation
# factor in the maximum-likelihood fit.
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_XimXibarp_signalMC"
  config.related_dataset = jpsi_data
  config.events          = 2000000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) ###
alg_name = "JpsiXiXibar"
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({ "ECMS" => [:double, 3.097] })
            .set_alias({ "std::vector<double>" => "Vdouble" })

# Selection chain: charged tracks, proton/pion PID, Lambda / Lambda_bar
# reconstruction via secondary vertex fits, then Xi- / Xibar+ via secondary
# vertex fits, and finally a 4C kinematic fit to the whole system.
event_selection = Selection.new
event_selection.select_track {
                  cos_theta 0.93        # |cos(theta_LAB)| < 0.93 (MDC acceptance)
                  Vz        100.0
                  Vr        10.0
                  nChrp    ">=3"        # at least 3 positively charged tracks
                  nChrn    ">=3"        # at least 3 negatively charged tracks
                }
                .pid(method: :probability) {
                  prob_cut 0.001
                  identify :proton, against: [:kaon, :pion]
                  identify :pion,   against: [:kaon, :proton]
                  nprp ">=1"
                  nprm ">=1"
                  npip ">=2"
                  npim ">=2"
                }
                .secondary_vertex_fit([:prp, :pim]) {
                  build_virtual_particle(:Lambda).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
                }
                .secondary_vertex_fit([:prm, :pip]) {
                  build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
                }
                .secondary_vertex_fit([:Lambda, :pim]) {
                  build_virtual_particle(:Xim).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
                }
                .secondary_vertex_fit([:Lambda_bar, :pip]) {
                  build_virtual_particle(:Xip).by_minimizing_mass_difference
                  remove_used_particle_from_candidate_list
                }
                # 4C kinematic fit imposing energy-momentum conservation on the
                # e+ e- -> J/psi -> Xi- Xibar+ system.
                .kinematic_fit([:Xim, :Xip]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 200
                }

# Non-expressible BOSS-side procedures captured as notes.
my_Algorithm
  .note(:proton_momentum_selection,
        "Proton candidates are required to have momentum > 0.3 GeV/c and pion " \
        "candidates to have momentum < 0.3 GeV/c, corresponding to the " \
        "non-overlapping momentum ranges of protons and pions in the signal " \
        "process.")
  .note(:xi_delta_cut,
        "Additional kinematic constraint delta < 0.016 GeV/c^2 applied per Xi " \
        "candidate, where delta = sqrt((m_{Lambda pi-} - m_{Xi-})^2 + " \
        "R^2 * (m_{p pi-} - m_{Lambda})^2), with R = sigma_{Xi-}/sigma_{Lambda} " \
        "the ratio of the m_{Lambda pi-} to m_{p pi-} mass resolutions. " \
        "The threshold is validated by S/sqrt(S+B) optimisation.")
  .note(:negative_decay_length_veto,
        "Events for which the Xi/Lambda reconstruction algorithm returns a " \
        "negative decay-length reading are rejected.")
  .note(:sideband_definition,
        "Two-dimensional sideband in m_{Lambda pi-} vs m_{Lambda_bar pi+} used " \
        "for background subtraction: lower band 1.274 < m < 1.306 GeV/c^2, " \
        "upper band 1.338 < m < 1.370 GeV/c^2.")

my_Algorithm.with_decay_card(decay_card_signal).apply(event_selection)
root_files = my_Algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])
