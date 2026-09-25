# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")        # J/psi (3.097 GeV) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")     # Corresponding inclusive MC sample

# Decay card for the signal process: J/psi -> p+ anti-p- phase space.
# The secondary elastic pp scattering on beam-pipe material is a Geant4 detector
# effect (not a decay), so it is NOT written in the decay card.
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 p+ anti-p- PHSP;
    Enddecay
    End
DECAYCARD

# 500k exclusive signal-MC events; secondary scattering of the generated proton
# on beam-pipe material is applied during the Geant4 simulation.
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_ppbar_secondary_scatter"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "JpsiPPbarScatter"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})   # double ECMS = 3.097 GeV

event_selection = Selection.new
event_selection.select_track {
                  cos_theta   0.93      # |cos(theta)| < 0.93
                  nChrp       "==2"     # exactly two positive tracks
                  nChrn       "==1"     # exactly one negative track -> three charged tracks in total
                }                       # no photon, Vz/Vr or net-charge criterion encoded
               .pid(method: :probability) {
                  prob_cut   0.001      # PID probability > 0.001
                  identify :proton, against: [:kaon, :pion]   # identify p+ and p-bar vs K and pi
                  nprp       ">=2"      # at least two proton candidates (scattered + target proton)
                  nprm       ">=1"      # at least one antiproton candidate
                }
               .secondary_vertex_fit([:prp, :prp]) {   # elastic pp-scattering vertex from the two final-state protons
                  by_minimizing_verfit_chi2
                  remove_used_particle_from_candidate_list
                }
               .partial_rec([2]) {   # reconstruct anti-p- (recID 2); infer the incident proton (recID 1) by recoil
                  require_recoil_mass 0.85, 1.02   # antiproton recoil-mass window (GeV/c^2)
                }
               # No global nC kinematic fit is encoded for this analysis.

my_algorithm
  .note(:secondary_scattering_simulation,
        "secondary elastic pp scattering of the generated proton on beam-pipe
         material (mineral-oil layer and MDC inner wall) is applied in the
         Geant4 simulation of the J/psi -> p+ anti-p- exclusive signal MC")
  .note(:scattering_vertex_rxy,
        "the elastic pp-scattering secondary vertex is required at
         Rxy 3.0-3.5 cm (mineral-oil layer) or 6.0-6.8 cm (MDC inner wall),
         i.e. +-4 sigma around the two scattering surfaces; this vertex-radius
         window has no dedicated DSL method and is applied at BOSS level")
  .note(:target_proton_momentum,
        "incident proton reconstructed by recoil p_p = p_e+ + p_e- - p_pbar and
         the target proton from four-momentum conservation; |p_target| < 50 MeV/c
         cut to suppress quasi-free nuclear collisions (no dedicated DSL method)")
  .note(:spin_precession,
        "proton spin precession in the 1.0 T solenoid between production and the
         scattering vertex (Rxy = 3.2 or 6.3 cm) is simulated with Geant4 and
         treated as a helicity-frame rotation")
  .with_decay_card(decay_card_signal).apply(event_selection)

# Execute the algorithm on real data, inclusive MC and the exclusive signal MC
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])