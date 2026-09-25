### Dataset description ###
# psi(2S) real data and inclusive MC at 3.686 GeV
psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")
# 3.65 GeV continuum real data and inclusive MC
cont_data  = DatasetManager.real_data.find("709_3650")
cont_incMC = DatasetManager.inclusive_mc.find("709_3650")

# Decay cards for the three signal modes: psi' -> gamma chi_c1 -> gamma gamma V
decay_card_phi = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c1 P2GC1;
    Enddecay
    Decay chi_c1
    1.000 gamma phi HELAMP 1 0 1.7 0 1.7 0 1 0;
    Enddecay
    Decay phi
    1.000 K+ K- VSS;
    Enddecay
    End
DECAYCARD

decay_card_rho = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c1 P2GC1;
    Enddecay
    Decay chi_c1
    1.000 gamma rho0 HELAMP 1 0 1.7 0 1.7 0 1 0;
    Enddecay
    Decay rho0
    1.000 pi+ pi- VSS;
    Enddecay
    End
DECAYCARD

decay_card_omega = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c1 P2GC1;
    Enddecay
    Decay chi_c1
    1.000 gamma omega HELAMP 1 0 1.7 0 1.7 0 1 0;
    Enddecay
    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay
    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

# Exclusive MC samples (200k events each), pegged to the psi(2S) real data
exMC_phi = DatasetManager.create_exclusive_mc do |config|
    config.sample_name    = "exmc_psip_gg_phi"
    config.related_dataset = psip_data
    config.events          = 200000
    config.decay_card      = decay_card_phi
    config.cross_section   = :default
end

exMC_rho = DatasetManager.create_exclusive_mc do |config|
    config.sample_name    = "exmc_psip_gg_rho0"
    config.related_dataset = psip_data
    config.events          = 200000
    config.decay_card      = decay_card_rho
    config.cross_section   = :default
end

exMC_omega = DatasetManager.create_exclusive_mc do |config|
    config.sample_name    = "exmc_psip_gg_omega"
    config.related_dataset = psip_data
    config.events          = 200000
    config.decay_card      = decay_card_omega
    config.cross_section   = :default
end

### Event selection (BOSS) — one Algorithm + one Selection chain per mode ###

# ---------- phi mode: psi' -> gamma chi_c1 -> gamma gamma phi, phi -> K+ K- ----------
alg_name_phi = "ChiCJGammaPhi"
alg_phi = Algorithm.new(alg_name_phi)
alg_phi.set_header(["#{alg_name_phi}Alg/#{alg_name_phi}.h"])
       .set_constant({"ECMS" => [:double, 3.686]})
       .note(:pid_correction_method, "PID is applied to the lower-momentum charged track only; " \
             "the partner track is taken as the other kaon without an independent PID decision")

phi_selection = Selection.new
    .select_track {                 # Charged track selection
        cos_theta 0.93              # |cos(theta)| < 0.93
        Vz        10.0              # |Vz| < 10 cm
        Vr        10.0              # Vr < 1 cm (10 mm)
        nChrp     "==1"             # exactly one positive track
        nChrn     "==1"             # exactly one negative track
        nNet      "==0"             # net charge zero
    }
    .select_photon {                # Photon selection
        tdc_emc_start     0
        tdc_emc_end       14
        angle_to_track    10.0      # at least 10 deg from any charged track
        energyThreshold_b 0.025     # 25 MeV in the barrel
        energyThreshold_e 0.050     # 50 MeV in the endcap
        nGam              ">=2"     # the two photons of the psi' / chi_c1 cascade
    }
    .pid(method: :probability) {    # PID (probability method)
        prob_cut 0.001
        identify :kaon, against: [:pion, :proton]   # K+ and K-
        nkp ">=1"
        nkm ">=1"
    }
    .kinematic_fit([:gamma, :gamma, :kp, :km]) {    # 4C fit to gamma gamma K+ K-
        nominal
        constrain_four_momentum
        chi2_cut 200            # loose cut; tighter chi2 <= 100 and mass sidebands applied in ROOT
    }

alg_phi.with_decay_card(decay_card_phi).apply(phi_selection)

# ---------- rho0 mode: psi' -> gamma chi_c1 -> gamma gamma rho0, rho0 -> pi+ pi- ----------
alg_name_rho = "ChiCJGammaRho"
alg_rho = Algorithm.new(alg_name_rho)
alg_rho.set_header(["#{alg_name_rho}Alg/#{alg_name_rho}.h"])
       .set_constant({"ECMS" => [:double, 3.686]})
       .note(:pid_correction_method, "PID is applied to the lower-momentum charged track only; " \
             "the partner track is taken as the other pion without an independent PID decision")

rho_selection = Selection.new
    .select_track {
        cos_theta 0.93
        Vz        10.0
        Vr        10.0
        nChrp     "==1"
        nChrn     "==1"
        nNet      "==0"
    }
    .select_photon {
        tdc_emc_start     0
        tdc_emc_end       14
        angle_to_track    10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam              ">=2"
    }
    .pid(method: :probability) {
        prob_cut 0.001
        identify :pion, against: [:kaon, :proton]   # pi+ and pi-
        npip ">=1"
        npim ">=1"
    }
    .kinematic_fit([:gamma, :gamma, :pip, :pim]) {  # 4C fit to gamma gamma pi+ pi-
        nominal
        constrain_four_momentum
        chi2_cut 200            # loose cut; tighter chi2 <= 100 and mass sidebands applied in ROOT
    }

alg_rho.with_decay_card(decay_card_rho).apply(rho_selection)

# ---------- omega mode: psi' -> gamma chi_c1 -> gamma gamma omega, omega -> pi+ pi- pi0, pi0 -> gamma gamma ----------
alg_name_omega = "ChiCJGammaOmega"
alg_omega = Algorithm.new(alg_name_omega)
alg_omega.set_header(["#{alg_name_omega}Alg/#{alg_name_omega}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})
         .note(:pid_correction_method, "PID is applied to the lower-momentum charged track only; " \
               "the partner track is taken as the other pion without an independent PID decision")

omega_selection = Selection.new
    .select_track {
        cos_theta 0.93
        Vz        10.0
        Vr        10.0
        nChrp     "==1"
        nChrn     "==1"
        nNet      "==0"
    }
    .select_photon {
        tdc_emc_start     0
        tdc_emc_end       14
        angle_to_track    10.0
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        nGam              ">=4"     # 2 cascade photons + 2 photons from pi0 -> gamma gamma
    }
    .pid(method: :probability) {
        prob_cut 0.001
        identify :pion, against: [:kaon, :proton]
        npip ">=1"
        npim ">=1"
    }
    .kalman_kinematic_fit([:gamma, :gamma]) {       # 1C mass-constrained fit reconstructing the pi0
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25
        npi0 ">=1"
    }
    .kinematic_fit([:gamma, :gamma, :pi0, :pip, :pim]) {  # 5C fit: 4C + pi0 mass constraint
        nominal
        constrain_four_momentum
        chi2_cut 200
    }

alg_omega.with_decay_card(decay_card_omega).apply(omega_selection)

### Execute on data, inclusive MC, continuum and exclusive MC ###
root_files_phi   = alg_phi.execute_on([psip_data, psip_incMC, cont_data, cont_incMC, exMC_phi])
root_files_rho   = alg_rho.execute_on([psip_data, psip_incMC, cont_data, cont_incMC, exMC_rho])
root_files_omega = alg_omega.execute_on([psip_data, psip_incMC, cont_data, cont_incMC, exMC_omega])