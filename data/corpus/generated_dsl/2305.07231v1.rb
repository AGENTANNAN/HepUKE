### ================================================================
###  J/psi -> Xi0 anti-Xi0  (double-tag) -- BOSS event selection
###    tag side : anti-Xi0 -> anti-Lambda0 pi0,
###               anti-Lambda0 -> anti-p- pi+,  pi0 -> gamma gamma
###    signal   : Xi0 -> K- e+   (Delta(B-L) = 0)
###               Xi0 -> K+ e-   (|Delta(B-L)| = 2)
### ================================================================

### Dataset preparation ###
jpsi_data  = DatasetManager.real_data.find("708_3097")       # J/psi @ 3.097 GeV real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")    # matching inclusive MC

# Decay card -- signal mode I : Xi0 -> K- e+
decay_card_modeI = <<~DECAYCARD
    Decay J/psi
    1.0000 Xi0 anti-Xi0 PHSP;
    Enddecay

    Decay anti-Xi0
    1.0000 anti-Lambda0 pi0 PHSP;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay Xi0
    1.0000 K- e+ PHSP;
    Enddecay

    End
DECAYCARD

# Decay card -- signal mode II : Xi0 -> K+ e-
decay_card_modeII = <<~DECAYCARD
    Decay J/psi
    1.0000 Xi0 anti-Xi0 PHSP;
    Enddecay

    Decay anti-Xi0
    1.0000 anti-Lambda0 pi0 PHSP;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+ HypWK;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay Xi0
    1.0000 K+ e- PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC samples (500k events per signal mode)
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_Xi0Xibar0_Ke"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_modeI
  config.cross_section   = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3097_Xi0Xibar0_Kecc"
  config.related_dataset = jpsi_data
  config.events          = 500000
  config.decay_card      = decay_card_modeII
  config.cross_section   = :default
end

### Event selection (tag side, shared by both signal modes) ###
event_selection_common = Selection.new
event_selection_common
  .select_track {                  # charged-track quality cuts
    cos_theta 0.93                 # |cos(theta)| < 0.93
    Vz        10.0                 # |Vz| < 10 cm
    Vr        1.0                  # Vr < 1 cm
    nChrp     "==2"                # pbar + pi+ + K + e  ->  2 positive tracks
    nChrn     "==2"                #                           2 negative tracks
    nNet      "==0"
  }
  .select_photon {                 # photons from pi0 -> gamma gamma
    energyThreshold_b 0.025        # E_gamma > 25 MeV in the barrel (|cos(theta)| < 0.80)
    energyThreshold_e 0.050
    nGam ">=2"
  }
  .pid(method: :probability) {     # probability method with 0.001 cut
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.6
    identify :proton, against: [:pion, :kaon]     # pbar of anti-Lambda0
    identify :pion,   against: [:kaon, :proton]   # pi+ of anti-Lambda0
    identify :kaon,   against: [:pion, :proton]   # K of the signal Xi0
    nprm ">=1"
    npip ">=1"
  }
  .secondary_vertex_fit([:prm, :pip]) {           # anti-Lambda0 -> anti-p- pi+ (secondary vertex)
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
  }
  .kinematic_fit([:prm, :pip, :gamma, :gamma]) {  # 4C kinematic fit of the tag side
    nominal
    constrain_four_momentum                                                   # nominal 4-momentum constraint
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)      # gamma gamma -> pi0 mass constraint
    invariant_mass_of(:prm, :pip).within(1.1107, 1.1207)                      # |M(pbar pi+) - M_Lambda| < 5 MeV
    invariant_mass_of(:gamma, :gamma).within(0.115, 0.150)                    # 115 < M(gamma gamma) < 150 MeV
    invariant_mass_of(:prm, :pip, :gamma, :gamma).within(1.2949, 1.3349)      # |M(Lambdabar pi0) - M_Xi0| < 20 MeV
    chi2_cut 200
  }

### Algorithm -- signal mode I : Xi0 -> K- e+ ###
alg_name_I = "Xi0Xibar0Ke"
alg_I = Algorithm.new(alg_name_I)
alg_I.set_header(["#{alg_name_I}Alg/#{alg_name_I}.h"])
     .set_constant({"ECMS" => [:double, 3.097]})
     .set_alias({"std::vector<double>" => "Vdouble"})
     .note(:pid_correction_method,
           "signal-side electron selected by the confidence-level criterion CL_e > 0.001 and
            CL_e/(CL_e+CL_pi+CL_K) > 0.8, and the signal kaon is taken as the hypothesis with the
            highest confidence level; these non-standard PID criteria are approximated here by
            identify_high_momentum_leptons and identify :kaon, which are not equivalent")
selection_I = event_selection_common.dup
alg_I.with_decay_card(decay_card_modeI).apply(selection_I)
root_files_I = alg_I.execute_on([jpsi_data, jpsi_incMC, exMC_modeI])

### Algorithm -- signal mode II : Xi0 -> K+ e- ###
alg_name_II = "Xi0Xibar0Kecc"
alg_II = Algorithm.new(alg_name_II)
alg_II.set_header(["#{alg_name_II}Alg/#{alg_name_II}.h"])
      .set_constant({"ECMS" => [:double, 3.097]})
      .set_alias({"std::vector<double>" => "Vdouble"})
      .note(:pid_correction_method,
            "signal-side electron selected by the confidence-level criterion CL_e > 0.001 and
             CL_e/(CL_e+CL_pi+CL_K) > 0.8, and the signal kaon is taken as the hypothesis with the
             highest confidence level; approximated here by identify_high_momentum_leptons and
             identify :kaon")
selection_II = event_selection_common.dup
alg_II.with_decay_card(decay_card_modeII).apply(selection_II)
root_files_II = alg_II.execute_on([jpsi_data, jpsi_incMC, exMC_modeII])