# frozen_string_literal: true

### Dataset preparation ###
# 37 R-scan energy points from 4.009 to 4.951 GeV (BOSS 703-707)
scan_point_names = %w[
  703_4009 703_4180 703_4190 703_4200 703_4210 703_4220 703_4230
  703_4237 703_4245 703_4246 703_4260 703_4270 703_4280 703_4310
  703_4360 703_4390 703_4420 703_4470 703_4530 703_4575 703_4600
  705_4130 705_4160 705_4290 705_4315 705_4340 705_4380 705_4400 705_4440
  706_4610 706_4620 706_4640 706_4660 706_4680 706_4700
  707_4740 707_4750 707_4780 707_4840 707_4914 707_4946
]
scan_points = scan_point_names.map { |name| DatasetManager.real_data.find(name) }

# Inclusive MC at 4.178 GeV (~40x the data luminosity)
incMC = DatasetManager.inclusive_mc.find("703_4180")

# ConExc decay card: continuum production with ISR modelling across the whole R-scan.
# The DSL auto-detects "ConExc", switches to the no-KKMC template and injects
# "Particle vpho <ECMS> 0.0" per energy point, so Particle vpho is omitted here.
# PHSP is a placeholder for the full amplitude model of e+e- -> p K- anti-Lambda.
decay_card_signal = <<~DECAYCARD
    Decay vpho
    1.0 ConExc -1 p+ K- anti-Lambda0;
    Enddecay

    Decay anti-Lambda0
    1.0 anti-p- pi+ PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive signal MC at every R-scan energy point
exMCs = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "exmc_pKLambdaBar"   # becomes exmc_pKLambdaBar_<boss>_<energy>
  config.events        = 100_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name  = "pKLambdaBar"
algorithm = Algorithm.new(alg_name)
algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 4.26]})   # per-point beam energy handled by the scan template

event_selection = Selection.new
  .select_track {
    cos_theta  0.93          # |cos(theta)| < 0.93
    Vz         20.0          # |Vz| < 20 cm
    nChrp      "==2"         # exactly two positive tracks
    nChrn      "==2"         # exactly two negative tracks
    nNet       "==0"         # net charge zero
  }
  # No PID: momenta are high, the p/pbar/pi+- assignment is resolved by the
  # secondary-vertex fit and the 4C fit.
  .assign({:chrgp => :prp, :chrgn => :pim})
  # Reconstruct the (anti-)Lambda from a p pi- pair; same chain for the charge-conjugate
  # mode (pbar pi+). Keep the candidate closest in mass to the nominal Lambda mass.
  .secondary_vertex_fit([:prp, :pim]) {
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # Keep exactly the best candidate inside 1.10 < M(p pi-) < 1.13 GeV/c^2
  .for_each(:Lambda) do
    where { mass < 1.10 || mass > 1.13 }
    remove
  end
  # Kaon candidate = remaining negative track; veto |cos(theta_K)| > 0.83
  # (suppresses beam-induced backgrounds)
  .assign({:chrgn => :km})
  .remove(:km) { condition "abs(cos_theta_of(:km)) > 0.83" }
  # 4C kinematic fit on the anti-Lambda p K- final state (charge-conjugate partner shares this chain)
  .kinematic_fit([:Lambda, :prp, :km]) {
    nominal                    # nominal fit: corrected four-momenta are the ones saved
    constrain_four_momentum    # 4C energy-momentum constraint against the measured CMS four-vector
    chi2_cut 100               # chi^2 < 100
  }

algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on the R-scan data points, the inclusive MC and the exclusive signal MC
root_files = algorithm.execute_on(scan_points + [incMC] + exMCs)