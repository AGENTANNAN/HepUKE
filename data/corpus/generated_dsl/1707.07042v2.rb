# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# psi(3686) signal sample (448.1 x 10^6 events) + inclusive MC
psip_data  = DatasetManager.real_data.find("709_3686")       # psi(2S) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")    # psi(2S) inclusive MC at 3.686 GeV

# 48 pb^-1 continuum sample at 3.65 GeV used as background control
cont_data  = DatasetManager.real_data.find("709_3650")       # continuum real data at 3.65 GeV
cont_incMC = DatasetManager.inclusive_mc.find("709_3650")    # continuum inclusive MC at 3.65 GeV

# ------------------------------------------------------------------
# Decay cards (EvtGen format): psi(2S) -> gamma chi_c0, chi_c0 -> eta(') eta(')
# ------------------------------------------------------------------

# Mode 1: both eta' -> gamma pi+ pi-
decay_card_mode1 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c0 PHSP;
    Enddecay

    Decay chi_c0
    1.000 eta' eta' PHSP;
    Enddecay

    Decay eta'
    1.000 gamma pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Mode 2: both eta' -> eta pi+ pi- (eta -> gamma gamma)
decay_card_mode2 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c0 PHSP;
    Enddecay

    Decay chi_c0
    1.000 eta' eta' PHSP;
    Enddecay

    Decay eta'
    1.000 eta pi+ pi- PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode 3: mixed eta' eta' (one eta' -> gamma pi+ pi-, the other -> eta pi+ pi-)
decay_card_mode3 = <<~DECAYCARD
    Alias etap2 eta'

    Decay psi(2S)
    1.000 gamma chi_c0 PHSP;
    Enddecay

    Decay chi_c0
    1.000 eta' etap2 PHSP;
    Enddecay

    Decay eta'
    1.000 gamma pi+ pi- PHSP;
    Enddecay

    Decay etap2
    1.000 eta pi+ pi- PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode 4: eta eta' with eta' -> gamma pi+ pi-
decay_card_mode4 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c0 PHSP;
    Enddecay

    Decay chi_c0
    1.000 eta eta' PHSP;
    Enddecay

    Decay eta'
    1.000 gamma pi+ pi- PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode 5: eta eta' with eta' -> eta pi+ pi- (both eta -> gamma gamma)
decay_card_mode5 = <<~DECAYCARD
    Alias eta2 eta

    Decay psi(2S)
    1.000 gamma chi_c0 PHSP;
    Enddecay

    Decay chi_c0
    1.000 eta eta' PHSP;
    Enddecay

    Decay eta'
    1.000 eta2 pi+ pi- PHSP;
    Enddecay

    Decay eta2
    1.000 gamma gamma PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

### Exclusive MC samples (200k events each, five modes) ###
exMC_mode1 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_psip_chic_etapetap_gpipi_gpipi"
    config.related_dataset = psip_data
    config.events          = 200_000
    config.decay_card      = decay_card_mode1
    config.cross_section   = :default
end

exMC_mode2 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_psip_chic_etapetap_etapipi_etapipi"
    config.related_dataset = psip_data
    config.events          = 200_000
    config.decay_card      = decay_card_mode2
    config.cross_section   = :default
end

exMC_mode3 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_psip_chic_etapetap_mixed"
    config.related_dataset = psip_data
    config.events          = 200_000
    config.decay_card      = decay_card_mode3
    config.cross_section   = :default
end

exMC_mode4 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_psip_chic_etaetap_gpipi"
    config.related_dataset = psip_data
    config.events          = 200_000
    config.decay_card      = decay_card_mode4
    config.cross_section   = :default
end

exMC_mode5 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_psip_chic_etaetap_etapipi"
    config.related_dataset = psip_data
    config.events          = 200_000
    config.decay_card      = decay_card_mode5
    config.cross_section   = :default
end

### Event selection (BOSS) — five independent chains ###
ECMS = [:double, 3.686]   # psi(3686) centre-of-mass energy (GeV)

# Part 1: Mode 1 — both eta' -> gamma pi+ pi-
# Final state: 4 charged tracks + >=3 photons (transition gamma + 2 gamma from eta')
alg1_name = "ChicEtaP2GPiPi"
alg1 = Algorithm.new(alg1_name)
alg1.set_header(["#{alg1_name}Alg/#{alg1_name}.h"])
    .set_constant({"ECMS" => ECMS})

sel1 = Selection.new
sel1.select_track {
        cos_theta 0.93      # |cos(theta)| < 0.93
        Vz        10.0      # |Vz| < 10 cm
        Vr        1.0       # Vr < 1 cm
        nChrp     "==2"     # 2 positive tracks
        nChrn     "==2"     # 2 negative tracks
        nNet      "==0"     # net charge zero
    }
    .select_photon {
        tdc_emc_start     0        # EMC time window 0 ...
        tdc_emc_end       14       # ... to 700 ns
        energyThreshold_b 0.025    # > 25 MeV in barrel
        energyThreshold_e 0.050    # > 50 MeV in endcaps
        nGam              ">=3"    # at least 3 photons
    }
    .assign({:chrgp => :pip, :chrgn => :pim})   # no PID: every track assumed pion
    .kinematic_fit([:gamma, :gamma, :gamma, :pip, :pip, :pim, :pim]) {
        nominal                  # nominal 4C fit (published tight chi2 cut applied later in ROOT)
        constrain_four_momentum  # constrain to beam four-momentum
        chi2_cut 200             # loose BOSS-level chi2 cut
    }
alg1.with_decay_card(decay_card_mode1).apply(sel1)

# Part 2: Mode 2 — both eta' -> eta pi+ pi- (eta -> gamma gamma)
# Final state: 4 charged tracks + >=5 photons (transition gamma + 4 gamma from 2 eta)
alg2_name = "ChicEtaP2EtaPiPi"
alg2 = Algorithm.new(alg2_name)
alg2.set_header(["#{alg2_name}Alg/#{alg2_name}.h"])
    .set_constant({"ECMS" => ECMS})

sel2 = Selection.new
sel2.select_track {
        cos_theta 0.93
        Vz        10.0
        Vr        1.0
        nChrp     "==2"
        nChrn     "==2"
        nNet      "==0"
    }
    .select_photon {
        tdc_emc_start     0
        tdc_emc_end       14
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam              ">=5"    # at least 5 photons
    }
    .assign({:chrgp => :pip, :chrgn => :pim})
    .kalman_kinematic_fit([:gamma, :gamma]) {   # build eta from photon pair
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
        chi2_cut 25
        neta ">=2"                              # require at least two eta candidates
    }
    .kinematic_fit([:gamma, :eta, :eta, :pip, :pip, :pim, :pim]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
    }
alg2.with_decay_card(decay_card_mode2).apply(sel2)

# Part 3: Mode 3 — mixed eta' eta' (one eta' -> gamma pi+ pi-, one -> eta pi+ pi-)
# Final state: 4 charged tracks + >=4 photons (transition gamma + eta' gamma + 2 gamma from eta)
alg3_name = "ChicEtaP2Mixed"
alg3 = Algorithm.new(alg3_name)
alg3.set_header(["#{alg3_name}Alg/#{alg3_name}.h"])
    .set_constant({"ECMS" => ECMS})

sel3 = Selection.new
sel3.select_track {
        cos_theta 0.93
        Vz        10.0
        Vr        1.0
        nChrp     "==2"
        nChrn     "==2"
        nNet      "==0"
    }
    .select_photon {
        tdc_emc_start     0
        tdc_emc_end       14
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam              ">=4"    # at least 4 photons
    }
    .assign({:chrgp => :pip, :chrgn => :pim})
    .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
        chi2_cut 25
        neta ">=1"                              # require at least one eta candidate
    }
    .kinematic_fit([:gamma, :gamma, :eta, :pip, :pip, :pim, :pim]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
    }
alg3.with_decay_card(decay_card_mode3).apply(sel3)

# Part 4: Mode 4 — eta eta' with eta' -> gamma pi+ pi-
# Final state: 2 charged tracks + >=4 photons (transition gamma + eta' gamma + 2 gamma from eta)
alg4_name = "ChicEtaEtaPGPiPi"
alg4 = Algorithm.new(alg4_name)
alg4.set_header(["#{alg4_name}Alg/#{alg4_name}.h"])
    .set_constant({"ECMS" => ECMS})

sel4 = Selection.new
sel4.select_track {
        cos_theta 0.93
        Vz        10.0
        Vr        1.0
        nChrp     "==1"     # 1 positive track
        nChrn     "==1"     # 1 negative track
        nNet      "==0"     # net charge zero
    }
    .select_photon {
        tdc_emc_start     0
        tdc_emc_end       14
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam              ">=4"    # at least 4 photons
    }
    .assign({:chrgp => :pip, :chrgn => :pim})
    .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
        chi2_cut 25
        neta ">=1"                              # require at least one eta candidate
    }
    .kinematic_fit([:gamma, :gamma, :eta, :pip, :pim]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
    }
alg4.with_decay_card(decay_card_mode4).apply(sel4)

# Part 5: Mode 5 — eta eta' with eta' -> eta pi+ pi- (both eta -> gamma gamma)
# Final state: 2 charged tracks + >=5 photons (transition gamma + 4 gamma from 2 eta)
alg5_name = "ChicEtaEtaPEtaPiPi"
alg5 = Algorithm.new(alg5_name)
alg5.set_header(["#{alg5_name}Alg/#{alg5_name}.h"])
    .set_constant({"ECMS" => ECMS})

sel5 = Selection.new
sel5.select_track {
        cos_theta 0.93
        Vz        10.0
        Vr        1.0
        nChrp     "==1"
        nChrn     "==1"
        nNet      "==0"
    }
    .select_photon {
        tdc_emc_start     0
        tdc_emc_end       14
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam              ">=5"    # at least 5 photons
    }
    .assign({:chrgp => :pip, :chrgn => :pim})
    .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
        chi2_cut 25
        neta ">=2"                              # require at least two eta candidates
    }
    .kinematic_fit([:gamma, :eta, :eta, :pip, :pim]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
    }
alg5.with_decay_card(decay_card_mode5).apply(sel5)

### Execute on datasets (signal data/MC, continuum control, and the exclusive MC sample) ###
# The published decay-mode-dependent chi2 cut (25-90) is applied AFTER the fit, i.e. in ROOT.
root_files_1 = alg1.execute_on([psip_data, psip_incMC, cont_data, cont_incMC, exMC_mode1])
root_files_2 = alg2.execute_on([psip_data, psip_incMC, cont_data, cont_incMC, exMC_mode2])
root_files_3 = alg3.execute_on([psip_data, psip_incMC, cont_data, cont_incMC, exMC_mode3])
root_files_4 = alg4.execute_on([psip_data, psip_incMC, cont_data, cont_incMC, exMC_mode4])
root_files_5 = alg5.execute_on([psip_data, psip_incMC, cont_data, cont_incMC, exMC_mode5])