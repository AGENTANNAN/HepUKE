### Dataset preparation ###
# Real data and inclusive MC at the two energy points
jpsi_data  = DatasetManager.real_data.find("708_3097")      # J/psi real data (225.2 M events) at 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # J/psi inclusive MC (Lund model) at 3.097 GeV
psip_data  = DatasetManager.real_data.find("709_3686")      # psi(2S) real data (1.06e8 events) at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # psi(2S) inclusive MC at 3.686 GeV

# Decay cards (EvtGen syntax) for the four modes
decay_card_jpsi_gamma = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma p+ anti-p- PHSP;
    Enddecay
    End
DECAYCARD

decay_card_jpsi_pi0 = <<~DECAYCARD
    Decay J/psi
    1.0000 pi0 p+ anti-p- PHSP;
    Enddecay
    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

decay_card_psip_gamma = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma p+ anti-p- PHSP;
    Enddecay
    End
DECAYCARD

decay_card_psip_pi0 = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi0 p+ anti-p- PHSP;
    Enddecay
    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

# Exclusive MC: 200k events per mode (the pi0 p pbar modes are the dominant backgrounds)
exMC_jpsi_gamma = DatasetManager.create_exclusive_mc do |config|
    config.sample_name = "exmc_jpsi_gamma_ppbar"
    config.related_dataset = jpsi_data
    config.events = 200000
    config.decay_card = decay_card_jpsi_gamma
    config.cross_section = :default
end

exMC_jpsi_pi0 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name = "exmc_jpsi_pi0_ppbar"    # background: J/psi -> pi0 p pbar
    config.related_dataset = jpsi_data
    config.events = 200000
    config.decay_card = decay_card_jpsi_pi0
    config.cross_section = :default
end

exMC_psip_gamma = DatasetManager.create_exclusive_mc do |config|
    config.sample_name = "exmc_psip_gamma_ppbar"
    config.related_dataset = psip_data
    config.events = 200000
    config.decay_card = decay_card_psip_gamma
    config.cross_section = :default
end

exMC_psip_pi0 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name = "exmc_psip_pi0_ppbar"    # background: psi(2S) -> pi0 p pbar
    config.related_dataset = psip_data
    config.events = 200000
    config.decay_card = decay_card_psip_pi0
    config.cross_section = :default
end

### Event selection (BOSS) ###
# Common selection chain shared by the J/psi and psi(2S) gamma p pbar channels
event_selection = Selection.new
    .select_track {                 # exactly one proton and one antiproton, net charge zero
        cos_theta 0.93              # |cos(theta)| < 0.93
        Vz 10.0                     # |Vz| < 10 cm
        Vr 1.0                      # Vr < 1 cm
        nChrp "==1"                 # exactly one positive track
        nChrn "==1"                 # exactly one negative track
        nNet  "==0"                 # net charge zero
    }
    .select_photon {                # at least one good photon
        tdc_emc_start 0             # EMC TDC in [0, 14]
        tdc_emc_end 14
        angle_to_track 10.0         # standard photon-track separation
        energyThreshold_b 0.025     # > 25 MeV in the barrel
        energyThreshold_e 0.050     # > 50 MeV in the endcap
        nGam ">=1"                  # at least one photon
    }
    .pid(method: :probability) {    # identify proton and antiproton against kaon and pion
        prob_cut 0.001
        identify :proton, against: [:kaon, :pion]   # p+ and anti-p- at once
        nprp "==1"
        nprm "==1"
    }
    .select_isolated_photon {       # photon isolated from the p and pbar tracks by > 30 deg
        angle_to_prp_track 30.0
        angle_to_prm_track 30.0
        nGam ">=1"
    }
    .kinematic_fit([:gamma, :prp, :prm]) {   # 4C kinematic fit to gamma p pbar
        nominal                     # nominal fit; fitted 4-momenta are used
        constrain_four_momentum     # 4-momentum conservation against the CMS energy
        chi2_cut 20                 # chi2 < 20
    }

# J/psi channel
alg_name_jpsi = "JpsiGammaPPbar"
alg_jpsi = Algorithm.new(alg_name_jpsi)
alg_jpsi.set_header(["#{alg_name_jpsi}Alg/#{alg_name_jpsi}.h"])
        .set_constant({"ECMS" => [:double, 3.097]})   # J/psi CMS energy
        .set_alias({"std::vector<double>" => "Vdouble"})
        .note(:missing_energy_cut, "|U_miss| < 0.05 GeV, U_miss computed from the two charged tracks (p and pbar); applied as a BOSS-level event selection")
        .note(:photon_transverse_momentum_cut, "P_t(gamma)^2 < 0.0005 (GeV/c)^2 for the radiative photon; BOSS-level selection")
        .note(:low_momentum_veto, "reject any event containing a charged track with momentum below 0.3 GeV/c")

alg_jpsi.with_decay_card(decay_card_jpsi_gamma).apply(event_selection)

# psi(2S) channel (same selection chain as the J/psi channel, different CMS energy)
alg_name_psip = "PsipGammaPPbar"
alg_psip = Algorithm.new(alg_name_psip)
alg_psip.set_header(["#{alg_name_psip}Alg/#{alg_name_psip}.h"])
        .set_constant({"ECMS" => [:double, 3.686]})   # psi(2S) CMS energy
        .set_alias({"std::vector<double>" => "Vdouble"})
        .note(:missing_energy_cut, "|U_miss| < 0.05 GeV, U_miss computed from the two charged tracks (p and pbar); applied as a BOSS-level event selection")
        .note(:photon_transverse_momentum_cut, "P_t(gamma)^2 < 0.0005 (GeV/c)^2 for the radiative photon; BOSS-level selection")
        .note(:low_momentum_veto, "reject any event containing a charged track with momentum below 0.3 GeV/c")

alg_psip.with_decay_card(decay_card_psip_gamma).apply(event_selection.dup)

# Execute on real data, inclusive MC, and the exclusive MC samples for each energy point
root_files_jpsi = alg_jpsi.execute_on([jpsi_data, jpsi_incMC, exMC_jpsi_gamma, exMC_jpsi_pi0])
root_files_psip = alg_psip.execute_on([psip_data, psip_incMC, exMC_psip_gamma, exMC_psip_pi0])