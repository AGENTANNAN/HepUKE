### Dataset description ###
# J/psi (3.097 GeV): real data and inclusive MC
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # corresponding inclusive MC

# ---- Decay card: J/psi -> gamma K_S0 K_S0 pi0 ----
decay_card_modeI = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma K_S0 K_S0 pi0 PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# ---- Decay card: J/psi -> gamma pi0 pi0 eta ----
decay_card_modeII = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma pi0 pi0 eta PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC: 200k events for each mode
exMC_modeI = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_gKSKSpi0"
  config.related_dataset = jpsi_data
  config.events         = 200000
  config.decay_card     = decay_card_modeI
  config.cross_section  = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_jpsi_gpi0pi0eta"
  config.related_dataset = jpsi_data
  config.events         = 200000
  config.decay_card     = decay_card_modeII
  config.cross_section  = :default
end

### Event selection (BOSS) ###

# ============================================================
# Mode I:  J/psi -> gamma K_S0 K_S0 pi0
#          (K_S0 -> pi+ pi-, pi0 -> gamma gamma)
# ============================================================
alg_name_I = "JpsiGammaKSKSpi0"
alg_I = Algorithm.new(alg_name_I)
alg_I.set_header(["#{alg_name_I}Alg/#{alg_name_I}.h"])
     .set_constant({"ECMS" => [:double, 3.097]})   # J/psi CMS energy (GeV)

sel_modeI = Selection.new
sel_modeI.select_track {
           cos_theta 0.93          # |cos(theta)| < 0.93
           nChrp "==2"             # two K_S0 -> pi+ pi- : 2 pi+
           nChrn "==2"             # two K_S0 -> pi+ pi- : 2 pi-
           nNet  "==0"             # net charge zero
         }
        .select_photon {
           tdc_emc_start 0         # EMC timing window (0 - 14 -> within 700 ns)
           tdc_emc_end   14
           angle_to_track 10.0     # opening angle > 10 deg to nearest charged track
           energyThreshold_b 0.025 # E > 25 MeV (barrel)
           energyThreshold_e 0.050 # E > 50 MeV (endcap)
           nGam ">=3"              # at least 3 photons (radiative gamma + pi0 -> gamma gamma)
         }
        .assign({:chrgp => :pip, :chrgn => :pim})   # no PID: all charged tracks taken as pions
        .secondary_vertex_fit([:pip, :pim]) {        # first K_S0 -> pi+ pi-
           build_virtual_particle(:K_S0).by_minimizing_mass_difference
           remove_used_particle_from_candidate_list
         }
        .secondary_vertex_fit([:pip, :pim]) {        # second K_S0 -> pi+ pi-
           build_virtual_particle(:K_S0).by_minimizing_mass_difference
           remove_used_particle_from_candidate_list
         }
        .kalman_kinematic_fit([:gamma, :gamma]) {    # pi0 -> gamma gamma (1C mass fit)
           invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
           chi2_cut 200
           npi0 ">=1"
         }
        .kinematic_fit([:gamma, :K_S0, :K_S0, :pi0]) {  # 4C to gamma K_S0 K_S0 pi0 (with the
           nominal                                       # mass-constrained pi0 -> effectively 5C)
           constrain_four_momentum
           chi2_cut 200    # loose cut in BOSS; tight 40 applied offline in ROOT
         }

alg_I.with_decay_card(decay_card_modeI).apply(sel_modeI)

# ============================================================
# Mode II: J/psi -> gamma pi0 pi0 eta
#          (pi0 -> gamma gamma, eta -> gamma gamma)
# ============================================================
alg_name_II = "JpsiGammaPi0Pi0Eta"
alg_II = Algorithm.new(alg_name_II)
alg_II.set_header(["#{alg_name_II}Alg/#{alg_name_II}.h"])
      .set_constant({"ECMS" => [:double, 3.097]})   # J/psi CMS energy (GeV)

sel_modeII = Selection.new
sel_modeII.select_track {
            cos_theta 0.93        # |cos(theta)| < 0.93
            nChrp "==0"           # no charged tracks allowed
            nChrn "==0"
            nNet  "==0"           # net charge zero
          }
         .select_photon {
            tdc_emc_start 0       # EMC timing window (0 - 14 -> within 700 ns)
            tdc_emc_end   14
            angle_to_track 10.0   # opening angle > 10 deg to nearest charged track
            energyThreshold_b 0.025
            energyThreshold_e 0.050
            nGam ">=7"            # at least 7 photons (radiative gamma + 2 pi0 + eta)
          }
         .kalman_kinematic_fit([:gamma, :gamma]) {    # pi0 -> gamma gamma (1C mass fit)
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
            chi2_cut 10
            npi0 ">=2"
          }
         .kalman_kinematic_fit([:gamma, :gamma]) {    # eta -> gamma gamma (1C mass fit)
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
            chi2_cut 200
            neta ">=1"
          }
         .kinematic_fit([:gamma, :pi0, :pi0, :eta]) {  # 4C to gamma pi0 pi0 eta
            nominal
            constrain_four_momentum
            chi2_cut 200  # loose cut in BOSS; tight 40 applied offline in ROOT
          }

alg_II.with_decay_card(decay_card_modeII).apply(sel_modeII)

### Execute on datasets ###
root_files_I  = alg_I.execute_on([jpsi_data, jpsi_incMC, exMC_modeI])
root_files_II = alg_II.execute_on([jpsi_data, jpsi_incMC, exMC_modeII])