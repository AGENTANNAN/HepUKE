# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description — R-scan e+e- -> n nbar (2.0 - 3.08 GeV) ###
# Center-of-mass energies (MeV) of the 713 R-scan points covering 2.0 - 3.08 GeV.
scan_energies = [2000, 2050, 2100, 2125, 2150, 2175, 2200, 2232,
                 2309, 2386, 2396, 2500, 2644, 2646, 2700, 2800,
                 2900, 2950, 2981, 3000, 3020, 3080]

# Real data and the matching inclusive MC sample at every scan point.
scan_data  = scan_energies.map { |e| DatasetManager.real_data.find("713_#{e}") }
scan_incMC = scan_energies.map { |e| DatasetManager.inclusive_mc.find("713_#{e}") }

# ConExc decay card for the continuum process e+e- -> n nbar (ISR modelling, mode 79).
# `Particle vpho` is intentionally omitted: the DSL injects it per energy point for a multi-energy scan.
decay_card_nnbar = <<~DECAYCARD
    Decay vpho
    1   ConExc 79;
    Enddecay
    End
DECAYCARD

# Signal exclusive MC: 100k events at each energy point, straight-line cross-section assumption.
exMCs_signal = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name   = "exmc_nnbar"        # -> exmc_nnbar_713_2000, exmc_nnbar_713_2050, ...
  config.events        = 100_000             # 100k events per energy point
  config.decay_card    = decay_card_nnbar    # ConExc n nbar card
  config.cross_section = :straight_line      # straight-line cross-section assumption
end

### Event selection (BOSS) ###
alg_name = "NNbarAnalysis"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.08]}) # nominal CMS energy; per-point value is taken from each dataset
            .set_alias({"std::vector<double>" => "Vdouble"})

# Fully neutral final state: require a charge-free event and at least two good photons.
# Neither PID nor a kinematic fit is applied (n / nbar are not charged tracks and no intermediate
# resonance is mass-constrained), so the selection chain ends at the photon step.
event_selection = Selection.new
event_selection.select_track {   # neutral-only: no charged tracks at all
                  nChrp  "==0"   # zero positively charged tracks
                  nChrn  "==0"   # zero negatively charged tracks
                  nNet   "==0"   # total charge zero
                }
               .select_photon {  # good EMC photons
                  tdc_emc_start     0      # EMC TDC window start
                  tdc_emc_end       14     # EMC TDC window end
                  energyThreshold_b 0.025  # barrel shower energy > 25 MeV
                  energyThreshold_e 0.050  # endcap shower energy > 50 MeV
                  nGam ">=2"               # at least two photons
                }

# Detector-level neutron/antineutron identification cannot be expressed as a track-level DSL
# cut; captured here so it can be propagated into the systematic-uncertainty evaluation.
my_algorithm
  .note(:neutron_identification,
        "neutron and antineutron identified from detector response (TOF timing and EMC shower "
        "energy/shape); a BDT trained on these inputs separates the n / nbar candidates from the "
        "photon and beam-induced backgrounds. No explicit track-level cut is applied in BOSS.")
  .note(:cosmic_veto,
        "MUC hit-pattern cosmic-ray veto applied to reject cosmic events that mimic the neutral "
        "n nbar final state.")
  .note(:event_categories,
        "selected candidates are split into three statistically independent categories A / B / C "
        "according to the n / nbar interaction topology; the categories are combined only at the "
        "cross-section extraction stage for systematic evaluation.")
  .with_decay_card(decay_card_nnbar)
  .apply(event_selection)

# Execute on every scan point (real data + inclusive MC) and the per-point signal MC.
root_files = my_algorithm.execute_on(scan_data + scan_incMC + exMCs_signal)