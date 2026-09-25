# Dataset preparation
jpsi_data   = DatasetManager.real_data.find("708_3097")
jpsi_incMC  = DatasetManager.inclusive_mc.find("708_3097")

# Decay card for signal: J/psi -> p pbar; the secondary scattering on beam pipe
# material is applied to the generated proton during Geant4 detector simulation.
decay_card_signal = <<~DECAYCARD
  Decay J/psi
  1.0000 p+ anti-p-      PHSP;
  Enddecay

  End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "jpsi_to_ppbar_signal"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

# Event selection (BOSS)
alg_name = "JpsiPPbarRelPhase"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({ "ECMS" => [:double, 3.097] })

event_selection = Selection.new
event_selection.select_track {
                 cos_theta 0.93        # |cos(theta)| < 0.93 for MDC tracks
                 nChrp    ">=2"        # two protons (p + p_pipe) expected
                 nChrn    ">=1"        # one anti-proton expected
                 nTot     "==3"        # exactly three charged tracks
               }
               .pid(method: :probability) {
                 prob_cut 0.001
                 # Proton hypothesis has the greatest likelihood: L(p) > L(K) and L(p) > L(pi)
                 identify :proton, against: [:kaon, :pion]
                 nprp ">=2"            # incident proton candidate + target proton knocked out
                 nprm ">=1"            # anti-proton from J/psi decay
               }

alg.note(:secondary_scattering_vertex,
         "Secondary elastic pp scattering vertex fit performed on the two final-state protons. " \
         "The reconstructed vertex is required to lie at the beam-pipe mineral-oil layer " \
         "(3.0 < Rxy < 3.5 cm) or at the MDC inner wall (6.0 < Rxy < 6.8 cm), corresponding " \
         "to +/-4 sigma windows around the two scattering surfaces.")
   .note(:recoil_kinematics,
         "Incident-proton four-momentum determined by recoil kinematics: " \
         "p_p = p_e+ + p_e- - p_pbar; then the target-proton momentum is computed " \
         "from four-momentum conservation between the two final-state protons and " \
         "the incident proton.")
   .note(:target_momentum_cut,
         "Require |p_target| < 50 MeV/c to suppress quasi-free nuclear collisions " \
         "(where the incident proton scatters off a bound nucleon with sizable Fermi motion).")
   .note(:pbar_recoil_mass_window,
         "Recoil mass against the anti-proton required in 0.85 < M_pbar^Recoil < 1.02 GeV/c^2 " \
         "(approximately +/-4 sigma proton-mass window).")
   .note(:spin_precession,
         "Proton spin precession in the 1.0 T BESIII solenoidal field between the production " \
         "point and the scattering vertex (Rxy = 3.2 or 6.3 cm) is simulated with Geant4 and " \
         "accounted for as a rotation of the helicity-frame coordinate system.")
   .note(:analyzing_power,
         "Analyzing power A_N(theta) for elastic pp scattering taken from the SAID database " \
         "and applied event-by-event in the downstream unbinned ML fit of sin(Delta Phi).")

alg.with_decay_card(decay_card_signal).apply(event_selection)
alg.execute_on([jpsi_data, jpsi_incMC, exMC_signal])
