# =============================================================================
# Born cross section of e+e- -> p pbar via ConExc
# Scan: 3.510 - 4.946 GeV (47 c.m. energies, ~26 fb-1 total)
# This file currently configures two representative points (4.260 GeV Y(4260)
# and 3.773 GeV psi(3770)); expand `scan_points` / `incMC_points` to the full
# 47-point scan for production.
# =============================================================================

### Dataset description ###
data_4260  = DatasetManager.real_data.find("703_4260")     # Y(4260), 4.260 GeV
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")  # matching inclusive MC
data_3773  = DatasetManager.real_data.find("712_3773")     # psi(3770), 3.773 GeV
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")  # matching inclusive MC

# Representative energy points; extended to the full 47-point scan for production.
scan_points  = [data_4260, data_3773]
incMC_points = [incMC_4260, incMC_3773]

### Decay card (ConExc) ###
# ConExc supplies the ISR (up to second order) and the vacuum-polarisation /
# measured sigma0(m) treatment. The DSL auto-detects the `ConExc` token, switches
# to the no-KKMC simulation template and injects `Particle vpho <ECMS> 0.0` per
# energy point -- so NO `Particle vpho` line is written here.
decay_card_ppbar = <<~DECAYCARD
    Decay vpho
    1.000 ConExc 24;
    Enddecay
    End
DECAYCARD

### Signal Monte Carlo: 100k-event exclusive ConExc MC per scan point ###
# Same signal MC (same decay card / cross section / event count) run over the
# distinct energy points -> use create_exclusive_mc_for.
exMCs_ppbar = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "exmc_ppbar_conexc"   # auto-suffixed per dataset
  config.events        = 100_000
  config.decay_card    = decay_card_ppbar
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "PpbarBornXS"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.260]})   # representative value; see note
            .set_alias({"std::vector<double>" => "Vdouble"})
            .note(:measured_cms_energy, "the nominal 4C kinematic fit is performed with the per-run measured c.m. energy (energy scan from 3.510 to 4.946 GeV); the ECMS constant is only a fallback / representative value")
            .note(:proton_ep_ratio, "the positive proton candidate is additionally required to satisfy E/p < 0.5 (per-track EMC energy over MDC momentum)")
            .note(:ppbar_opening_angle, "the opening angle between the p and pbar three-momenta is required to be larger than 3.1 rad")
            .note(:muc_hit_depth, "both the p and pbar tracks are required to have a muon-counter hit depth < 40 cm")
            .note(:efficiency_curve, "detection efficiency and the ISR / vacuum-polarisation correction factors are extracted iteratively from the ConExc signal MC at each energy point")

event_selection = Selection.new
event_selection.select_track {          # charged-track quality cuts (no photon cuts in this analysis)
                  cos_theta 0.93        # |cos(theta)| < 0.93
                  Vz        10.0        # |Vz| < 10 cm
                  Vr        1.0         # Vr < 1 cm
                  nChrp     "==1"       # exactly one positive track
                  nChrn     "==1"       # exactly one negative track
                  nNet      "==0"       # net charge zero
                }
               .pid(method: :probability) {   # probability PID, 0.001 threshold
                  prob_cut 0.001
                  identify :proton, against: [:kaon, :pion]  # separate p+/p- from K and pi
                  nprp "==1"        # exactly one proton
                  nprm "==1"        # exactly one anti-proton
                }
               .kinematic_fit([:prp, :prm]) {   # nominal 4C fit to p+ p-
                  nominal
                  constrain_four_momentum       # constrain to the (per-run measured) c.m. four-momentum
                  chi2_cut 200                  # loose chi2 cut; tight cut applied in ROOT
                }

# Attach decay card (defines the kinematic variables) and render the selection
my_algorithm.with_decay_card(decay_card_ppbar).apply(event_selection)

# Execute on real data + inclusive MC (per scan point) + ConExc signal MC
root_files = my_algorithm.execute_on(scan_points + incMC_points + exMCs_ppbar)