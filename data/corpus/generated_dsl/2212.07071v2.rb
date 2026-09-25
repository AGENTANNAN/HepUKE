# frozen_string_literal: true

### Dataset preparation ###
# Neutron EM form factors via e+e- -> n nbar in the continuum.
# 12 scan points (2.0000 - 2.9500 GeV); the extraction groups them into 5 c.m. energy intervals.
scan_names = %w[
  713_Rscan_2000
  713_Rscan_2050
  713_Rscan_2125
  713_Rscan_2150
  713_Rscan_2175
  713_Rscan_2200
  713_Rscan_2232
  713_Rscan_2309
  713_Rscan_2386
  713_Rscan_2396
  713_Rscan_2644
  713_Rscan_2950
]

data_points = scan_names.map { |name| DatasetManager.real_data.find(name) }   # real data at each point
incMC_points = scan_names.map { |name| DatasetManager.inclusive_mc.find(name) } # inclusive MC at each point

# ConExc decay card (NLO: ISR + vacuum polarisation, Coulomb factor 1) for e+e- -> n nbar.
# The DSL auto-detects the `ConExc` token, switches to the no-KKMC template and injects
# `Particle vpho <ECMS>` for every energy point of the scan -> no explicit Particle vpho line.
decay_card_nnbar = <<~DECAYCARD
    Decay vpho
    1.0000   ConExc   1   n0   anti-n0 ;
    Enddecay
    End
DECAYCARD

# 500k-event signal exclusive MC per energy point (single shared card / cross section)
exMC_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "nnbar_signal_exclusive_mc"
  config.events        = 500_000
  config.decay_card    = decay_card_nnbar
  config.cross_section = :default
end
exMC_signal.each { |m| m.save_to_config(format: :yaml, file_path: 'nnbar_signal_mc') }

### Event selection (BOSS) ###
alg_name = "NeutronFF"
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({ "ECMS" => [:double, 2.0] }) # representative value; per-point ECMS is injected by the ConExc mechanism
            .note(:neutron_four_momentum, "The neutron and antineutron four-momenta are taken directly from the EMC shower energy and direction; no BOSS kinematic fit is applied (fully neutral, counting measurement).")
            .note(:neutron_tof_selection, "TOF requirement: hits within the azimuthal span of six scintillators along the neutron momentum; not expressible in the DSL.")
            .note(:event_category, "Events are classified into category A (TOF and EMC shower from both n and nbar), B (EMC from both, TOF from nbar only), C (EMC from both, no TOF).")
            .note(:emc_shower_window, "The nbar (n) candidate EMC shower must have deposited energy in [0.5, 2.0] GeV and |cos(theta)| < 0.7; the upper energy bound and the cos(theta) cut are not expressible in select_photon.")
            .note(:efficiency_calibration, "Neutron/antineutron detection efficiencies are calibrated data-driven with J/psi -> pbar pi+ n (p+ pi- nbar) control samples (10087M J/psi events, sample 708_3097).")
            .note(:trigger_correction, "Trigger efficiency correction derived from the EMC-based neutral trigger.")
            .note(:efficiency_curve, "Form-factor-model neutron efficiency determined iteratively, converging within 1%.")
            .note(:background_veto, "Beam-related background and e+e- -> gamma gamma background are accounted for in the composite maximum-likelihood (NLL/MIGRAD) yield fit per cos(theta_nbar) bin (7 equidistant bins within -0.7 < cos(theta_nbar) < 0.7).")

event_selection = Selection.new
event_selection
  .select_track {              # fully neutral final state
    nChrp  "==0"               # zero positively charged tracks in the MDC
    nChrn  "==0"               # zero negatively charged tracks in the MDC
    nTot   "==0"               # no charged tracks at all
  }
  .select_photon {             # EMC showers from n and nbar
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    energyThreshold_b 0.5      # EMC barrel shower energy > 0.5 GeV
    energyThreshold_e 0.5      # EMC endcap shower energy > 0.5 GeV
    nGam ">=2"                 # at least two EMC showers (n and nbar candidates)
  }

my_Algorithm.with_decay_card(decay_card_nnbar).apply(event_selection)
root_files = my_Algorithm.execute_on(data_points + incMC_points + exMC_signal)