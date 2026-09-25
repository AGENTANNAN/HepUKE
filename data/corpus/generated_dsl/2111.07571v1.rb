# Core DSL classes and dependencies are loaded automatically at execution.

### Dataset preparation ###
# Real data and inclusive MC
jpsi_data  = DatasetManager.real_data.find("708_3097")     # 2017 J/psi real data, sqrt(s) = 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # J/psi inclusive MC
cont_data  = DatasetManager.real_data.find("708_3080")     # 3.08 GeV continuum real data (background)
psip_data  = DatasetManager.real_data.find("709_3686")     # psi(3686) real data
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # psi(3686) inclusive MC

# Decay card for the efficiency channel: psi(3686) -> pi+ pi- J/psi, J/psi -> gamma e+ e-
decay_card_psip_jpsi = <<~DECAYCARD
  Decay psi(2S)
  1.0000 pi+ pi- J/psi PHSP;
  Enddecay

  Decay J/psi
  1.0000 gamma e+ e- PHSP;
  Enddecay

  End
DECAYCARD

# Decay card for the luminosity process e+e- -> gamma gamma.
# No intermediate resonance is specified, so the BESIII KKMC convention top mother psi(4260) is used.
decay_card_gg = <<~DECAYCARD
  Decay psi(4260)
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Exclusive MC for the efficiency study: 10M psi(3686) -> pi+ pi- J/psi, J/psi -> gamma e+ e-
exMC_psip_jpsi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_pipi_jpsi_gammaee"
  config.related_dataset = psip_data
  config.events          = 10_000_000
  config.decay_card      = decay_card_psip_jpsi
  config.cross_section   = :default
end

### Event selection (BOSS) - inclusive J/psi yield at 3.097 GeV ###
alg_name_inc = "InclusiveJpsiYield"
alg_inc = Algorithm.new(alg_name_inc)
alg_inc.set_header(["#{alg_name_inc}Alg/#{alg_name_inc}.h"])
       .set_constant({"ECMS" => [:double, 3.097]})

# Inclusive selection: every event passing the quality cuts is counted.
sel_inc = Selection.new
sel_inc.select_track {
          cos_theta 0.93      # |cos(theta)| < 0.93
          Vz        15.0      # |Vz| < 15 cm
          Vr        1.0       # Vr < 1 cm
          nChrp     ">=1"     # at least one positive track
          nChrn     ">=1"     # at least one negative track
          nTot      ">=2"     # at least two tracks in total
        }
        .select_photon {
          tdc_emc_start     0      # EMC time window 0-700 ns
          tdc_emc_end       14
          energyThreshold_b 0.025  # barrel |cos(theta)| < 0.83, E > 25 MeV
          energyThreshold_e 0.050  # endcap 0.86 < |cos(theta)| < 0.93, E > 50 MeV
          angle_to_track    10.0   # more than 10 deg from any charged track
        }

alg_inc
  .note(:track_momentum_cut, "all charged-track momenta required < 2.0 GeV/c; for events with exactly two tracks, additionally each momentum < 1.5 GeV/c and each EMC shower energy < 1.0 GeV to reject Bhabha and dimuon backgrounds. These per-track momentum/energy conditions are not expressible in select_track.")
  .note(:visible_energy_cut, "total visible energy required > 1.0 GeV; no dedicated selection primitive for this global quantity.")
  .note(:no_pid_no_kinematic_fit, "inclusive J/psi selection applies no PID and no kinematic fit; the J/psi is counted inclusively. Accordingly no nominal kinematic fit is declared.")
  .with_decay_card(decay_card_psip_jpsi)
  .apply(sel_inc)

alg_inc.execute_on([jpsi_data, jpsi_incMC, cont_data])

### Event selection (BOSS) - J/psi detection efficiency from psi(3686) -> pi+ pi- J/psi ###
alg_name_eff = "PsipJpsiEff"
alg_eff = Algorithm.new(alg_name_eff)
alg_eff.set_header(["#{alg_name_eff}Alg/#{alg_name_eff}.h"])
       .set_constant({"ECMS" => [:double, 3.686]})

# Tag psi(3686) -> pi+ pi- J/psi with at least two oppositely charged soft pions
sel_eff = Selection.new
sel_eff.select_track {
          cos_theta 0.93      # |cos(theta)| < 0.93
          Vz        15.0      # |Vz| < 15 cm
          Vr        1.0       # Vr < 1 cm
          nChrp     ">=1"     # at least one soft pi+
          nChrn     ">=1"     # at least one soft pi-
        }

alg_eff
  .note(:soft_pion_selection, "efficiency measured on psi(3686) -> pi+ pi- J/psi requiring at least two oppositely charged soft pions with momentum < 0.4 GeV/c; the soft-pion momentum condition is not expressible in select_track.")
  .note(:recoil_mass_fit, "the J/psi detection efficiency is extracted by fitting the pi+ pi- recoil-mass spectrum (ROOT-level step, applied after this BOSS selection).")
  .with_decay_card(decay_card_psip_jpsi)
  .apply(sel_eff)

alg_eff.execute_on([psip_data, psip_incMC, exMC_psip_jpsi])

### Event selection (BOSS) - luminosity from e+e- -> gamma gamma ###
alg_name_lumi = "LuminosityGG"
alg_lumi = Algorithm.new(alg_name_lumi)
alg_lumi.set_header(["#{alg_name_lumi}Alg/#{alg_name_lumi}.h"])
        .set_constant({"ECMS" => [:double, 3.097]})

sel_lumi = Selection.new
sel_lumi.select_photon {
           nGam ">=2"    # at least two EMC showers
         }

alg_lumi
  .note(:gg_luminosity_selection, "luminosity from e+e- -> gamma gamma: at least two EMC showers, |cos(theta)| < 0.8, the second-most-energetic shower energy between 1.2 and 1.6 GeV, and |Delta phi| < 2.5 deg. The shower-energy ordering, per-shower angle and azimuthal-correlation conditions are not expressible with the current selection primitives.")
  .with_decay_card(decay_card_gg)
  .apply(sel_lumi)

alg_lumi.execute_on([jpsi_data, cont_data, psip_data])