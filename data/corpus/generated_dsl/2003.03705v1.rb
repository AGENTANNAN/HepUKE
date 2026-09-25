# frozen_string_literal: true

### Dataset preparation ###
# e+e- -> eta J/psi Born cross section scan, 3.810 - 4.600 GeV (BOSS 7.0.3 samples)
scan_names = %w[
  703_3810 703_3872 703_3900 703_4009 703_4090 703_4180 703_4190
  703_4200 703_4210 703_4220 703_4230 703_4237 703_4245 703_4246
  703_4260 703_4270 703_4280 703_4310 703_4360 703_4390 703_4420
  703_4470 703_4530 703_4575 703_4600
]
scan_data  = scan_names.map { |name| DatasetManager.real_data.find(name) }     # real data at each energy point
scan_incMC = scan_names.map { |name| DatasetManager.inclusive_mc.find(name) }  # matching inclusive MC

# Decay card for the signal, written for the ConExc generator (continuum / ISR model).
# ConExc mode 80 corresponds to the eta J/psi final state.
# `Particle vpho` is intentionally omitted: for a multi-energy scan the DSL injects
# `Particle vpho <ECMS> 0.0` once per energy point automatically.
decay_card_signal = <<~DECAYCARD
    Decay vpho
    1.000 ConExc 80;
    Enddecay

    Decay J/psi
    1.000 e+ e- PHOTOS VLL;
    Enddecay

    Decay eta
    0.3935 gamma gamma PHSP;
    0.2268 pi+ pi- pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# One 100k-event exclusive MC per scan point; the same sample serves both modes.
exMCs_signal = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name   = "exmc_etajpsi_conexc80_ee"   # auto-suffixed per energy point
  config.events        = 100_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
# Common track + photon selection shared by Mode I and Mode II
common_selection = Selection.new
    .select_track {
        cos_theta 0.93      # |cos(theta)| < 0.93
        Vz        10.0      # |Vz| < 10 cm
        Vr        1.0       # Vr < 1 cm
        nNet      "==0"     # net charge zero
    }
    .select_photon {
        tdc_emc_start     0      # EMC TDC window [0, 14]
        tdc_emc_end       14
        energyThreshold_b 0.025  # E > 25 MeV in the barrel
        energyThreshold_e 0.050  # E > 50 MeV in the endcap
        angle_to_track    10.0   # at least 10 degrees away from any charged track
        nGam              ">=2"  # at least two good photons
    }

# ---------------- Mode I: eta -> gamma gamma, J/psi -> e+ e- ----------------
alg_name_modeI = "EtaJpsiEEGG"
alg_modeI = Algorithm.new(alg_name_modeI)
alg_modeI.set_header(["#{alg_name_modeI}Alg/#{alg_name_modeI}.h"])
         .set_constant({"ECMS" => [:double, 4.260]})   # nominal scan energy; overridden per energy point at execution
         .set_alias({"std::vector<double>" => "Vdouble"})

modeI_selection = common_selection.dup
    .select_track {
        nChrp "==1"   # exactly one positive track
        nChrn "==1"   # exactly one negative track
        nNet  "==0"
    }
    .pid(method: :probability) {
        # high-momentum (p > 1.0 GeV) tracks treated as leptons:
        # electron if EMC energy > 0.6 GeV, otherwise muon
        identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                       treat_as_electron_if_energy_above: 0.6
        nlp "==1"     # exactly one positive lepton
        nlm "==1"     # exactly one negative lepton
    }
    .kinematic_fit([:lp, :lm, :gamma, :gamma]) {
        nominal                                          # nominal fit (corrected four-momenta)
        constrain_four_momentum                          # 4C energy-momentum conservation
        invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)  # M(e+e-) constrained to J/psi
        chi2_cut 200                                     # loose chi2 cut; tight cut applied in ROOT
    }

alg_modeI.with_decay_card(decay_card_signal).apply(modeI_selection)

# ---------------- Mode II: eta -> pi+ pi- pi0, pi0 -> gamma gamma, J/psi -> e+ e- ----------------
alg_name_modeII = "EtaJpsiEEPiPiPi0"
alg_modeII = Algorithm.new(alg_name_modeII)
alg_modeII.set_header(["#{alg_name_modeII}Alg/#{alg_name_modeII}.h"])
          .set_constant({"ECMS" => [:double, 4.260]})
          .set_alias({"std::vector<double>" => "Vdouble"})

modeII_selection = common_selection.dup
    .select_track {
        nChrp "==2"   # exactly two positive tracks
        nChrn "==2"   # exactly two negative tracks
        nNet  "==0"
    }
    .pid(method: :probability) {
        prob_cut 0.001
        # high-momentum (p > 1.0 GeV) tracks treated as leptons:
        # electron if EMC energy > 0.6 GeV, otherwise muon
        identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                       treat_as_electron_if_energy_above: 0.6
        identify :pion, against: [:kaon, :proton]   # remaining tracks: pi+ / pi- vs K and p
        nlp  "==1"    # exactly one positive lepton
        nlm  "==1"    # exactly one negative lepton
        npip ">=1"    # at least one pi+
        npim ">=1"    # at least one pi-
    }
    .remove([:lp <= :chrgp, :lm <= :chrgn])   # remove the identified leptons from the charged lists
    .kalman_kinematic_fit([:gamma, :gamma]) { # reconstruct pi0 from a photon pair (1C mass constraint)
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25
        npi0 ">=1"   # at least one pi0 candidate
    }
    .kinematic_fit([:lp, :lm, :pip, :pim, :pi0]) {
        nominal                                          # nominal fit (corrected four-momenta)
        constrain_four_momentum                          # 4C energy-momentum conservation
        invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)  # M(e+e-) constrained to J/psi
        chi2_cut 200                                     # loose chi2 cut; tight cut applied in ROOT
    }

alg_modeII.with_decay_card(decay_card_signal).apply(modeII_selection)

### Execution on data, inclusive MC and signal exclusive MC ###
root_files_modeI  = alg_modeI.execute_on(scan_data + scan_incMC + exMCs_signal)
root_files_modeII = alg_modeII.execute_on(scan_data + scan_incMC + exMCs_signal)