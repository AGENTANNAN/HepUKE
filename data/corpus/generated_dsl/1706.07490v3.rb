# ============================================================================
#  BOSS-side DSL
#  e+e- -> phi phi omega   and   e+e- -> phi phi phi
#  (partial reconstruction: the phi phi pair is reconstructed and the third
#   particle -- the omega or the third phi -- is tagged only through its
#   recoil mass against the phi phi system; no kinematic fit is applied)
# ============================================================================

### ------------------------------ Datasets ------------------------------ ###
# Six BESIII CMS energies: 4008, 4226, 4258, 4358, 4416 and 4600 MeV
data_4008 = DatasetManager.real_data.find("703_4009")   # sqrt(s) ~ 4008 MeV
data_4226 = DatasetManager.real_data.find("703_4230")   # sqrt(s) ~ 4226 MeV
data_4258 = DatasetManager.real_data.find("703_4260")   # sqrt(s) ~ 4258 MeV
data_4358 = DatasetManager.real_data.find("703_4360")   # sqrt(s) ~ 4358 MeV
data_4416 = DatasetManager.real_data.find("703_4420")   # sqrt(s) ~ 4416 MeV
data_4600 = DatasetManager.real_data.find("703_4600")   # sqrt(s) ~ 4600 MeV

incMC_4008 = DatasetManager.inclusive_mc.find("703_4009")
incMC_4226 = DatasetManager.inclusive_mc.find("703_4230")
incMC_4258 = DatasetManager.inclusive_mc.find("703_4260")
incMC_4358 = DatasetManager.inclusive_mc.find("703_4360")
incMC_4416 = DatasetManager.inclusive_mc.find("703_4420")
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")

data_points  = [data_4008, data_4226, data_4258, data_4358, data_4416, data_4600]
incMC_points = [incMC_4008, incMC_4226, incMC_4258, incMC_4358, incMC_4416, incMC_4600]

### ----------------------------- Decay cards ---------------------------- ###
# Signal: e+e- -> phi phi omega, generated flat in phase space.
# KKMC ISR and the psi(4260) top-mother convention are the DSL defaults.
decay_card_phiphiomega = <<~DECAYCARD
    Decay psi(4260)
    1.0 phi phi omega PHSP;
    Enddecay

    Decay phi
    1.0 K+ K- VSS;
    Enddecay

    Decay omega
    1.0 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.0 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Signal: e+e- -> phi phi phi, generated flat in phase space.
decay_card_phiphiphi = <<~DECAYCARD
    Decay psi(4260)
    1.0 phi phi phi PHSP;
    Enddecay

    Decay phi
    1.0 K+ K- VSS;
    Enddecay

    End
DECAYCARD

# Peaking background: e+e- -> K+ K- phi phi
decay_card_bkg_KpKm_phiphi = <<~DECAYCARD
    Decay psi(4260)
    1.0 K+ K- phi phi PHSP;
    Enddecay

    Decay phi
    1.0 K+ K- VSS;
    Enddecay

    End
DECAYCARD

# Peaking background: e+e- -> K+ K- K+ K- phi
decay_card_bkg_KpKmKpKm_phi = <<~DECAYCARD
    Decay psi(4260)
    1.0 K+ K- K+ K- phi PHSP;
    Enddecay

    Decay phi
    1.0 K+ K- VSS;
    Enddecay

    End
DECAYCARD

# Peaking background: e+e- -> K+ K- phi omega
decay_card_bkg_KpKm_phiomega = <<~DECAYCARD
    Decay psi(4260)
    1.0 K+ K- phi omega PHSP;
    Enddecay

    Decay phi
    1.0 K+ K- VSS;
    Enddecay

    Decay omega
    1.0 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.0 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Peaking background: e+e- -> K+ K- K+ K- omega
decay_card_bkg_KpKmKpKm_omega = <<~DECAYCARD
    Decay psi(4260)
    1.0 K+ K- K+ K- omega PHSP;
    Enddecay

    Decay omega
    1.0 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.0 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

### ---------------------------- Exclusive MC ---------------------------- ###
# Signals: 500k events per mode per energy point (same signal MC at each of
# the six energy points -> create_exclusive_mc_for over the six datasets).
exMC_phiphiomega = DatasetManager.create_exclusive_mc_for(data_points) do |config|
    config.sample_name   = "exmc_phiphiomega"
    config.events        = 500_000
    config.decay_card    = decay_card_phiphiomega
    config.cross_section = :default
end

exMC_phiphiphi = DatasetManager.create_exclusive_mc_for(data_points) do |config|
    config.sample_name   = "exmc_phiphiphi"
    config.events        = 500_000
    config.decay_card    = decay_card_phiphiphi
    config.cross_section = :default
end

# Peaking backgrounds: 100k events per mode per energy point.
exMC_bkg_KpKm_phiphi = DatasetManager.create_exclusive_mc_for(data_points) do |config|
    config.sample_name   = "exmc_bkg_KpKm_phiphi"
    config.events        = 100_000
    config.decay_card    = decay_card_bkg_KpKm_phiphi
    config.cross_section = :default
end

exMC_bkg_KpKmKpKm_phi = DatasetManager.create_exclusive_mc_for(data_points) do |config|
    config.sample_name   = "exmc_bkg_KpKmKpKm_phi"
    config.events        = 100_000
    config.decay_card    = decay_card_bkg_KpKmKpKm_phi
    config.cross_section = :default
end

exMC_bkg_KpKm_phiomega = DatasetManager.create_exclusive_mc_for(data_points) do |config|
    config.sample_name   = "exmc_bkg_KpKm_phiomega"
    config.events        = 100_000
    config.decay_card    = decay_card_bkg_KpKm_phiomega
    config.cross_section = :default
end

exMC_bkg_KpKmKpKm_omega = DatasetManager.create_exclusive_mc_for(data_points) do |config|
    config.sample_name   = "exmc_bkg_KpKmKpKm_omega"
    config.events        = 100_000
    config.decay_card    = decay_card_bkg_KpKmKpKm_omega
    config.cross_section = :default
end

# Flatten: each create_exclusive_mc_for call returns an Array<ExclusiveMC>
all_exMC = exMC_phiphiomega + exMC_phiphiphi +
           exMC_bkg_KpKm_phiphi + exMC_bkg_KpKmKpKm_phi +
           exMC_bkg_KpKm_phiomega + exMC_bkg_KpKmKpKm_omega

### ------------------------ Event selection (BOSS) ---------------------- ###
alg_name_phiphiomega = "PhiPhiOmegaRecoil"
alg_phiphiomega = Algorithm.new(alg_name_phiphiomega)
alg_phiphiomega.set_header(["#{alg_name_phiphiomega}Alg/#{alg_name_phiphiomega}.h"])
               .set_constant({"ECMS" => [:double, 4.600]})
               .note(:per_energy_ecms, "the analysis runs at six CMS energies " \
                     "(4008, 4226, 4258, 4358, 4416, 4600 MeV); the ECMS constant " \
                     "must be set to the matching energy point for each job so that " \
                     "the recoil mass against the phi phi system is correct")

alg_name_phiphiphi = "PhiPhiPhiRecoil"
alg_phiphiphi = Algorithm.new(alg_name_phiphiphi)
alg_phiphiphi.set_header(["#{alg_name_phiphiphi}Alg/#{alg_name_phiphiphi}.h"])
             .set_constant({"ECMS" => [:double, 4.600]})
             .note(:per_energy_ecms, "the analysis runs at six CMS energies " \
                   "(4008, 4226, 4258, 4358, 4416, 4600 MeV); the ECMS constant " \
                   "must be set to the matching energy point for each job so that " \
                   "the recoil mass against the phi phi system is correct")

# Common selection for both channels (identical phi phi reconstruction).
event_selection = Selection.new
    .select_track {
        cos_theta  0.93      # |cos(theta)| < 0.93
        Vz         10.0      # |Vz| < 10 cm
        Vr         1.0       # Vr < 1 cm
        nChrp      "==2"     # exactly two positive tracks
        nChrn      "==2"     # exactly two negative tracks
        nNet       "==0"     # net charge zero
    }
    .pid(method: :probability) {
        prob_cut   0.001                            # probability method, 0.001 cut
        identify :kaon, against: [:pion, :proton]   # separate K from pi/p (K+ and K- at once)
        nkp        "==2"                            # exactly two K+
        nkm        "==2"                            # exactly two K-
    }
    # Partial reconstruction (no kinematic fit): the two phi -> K+K- are
    # reconstructed and the third particle (omega, or the third phi, recID 3)
    # is left undetected and tagged only through its recoil mass against the
    # phi phi system.
    .partial_miss([3]) {
        # retain only the phi phi pairing that minimises
        # (M(K+K-) - M_phi)^2 summed over both phi candidates (PDG phi mass).
        best_combination_by_mass :phi, 1.019
        # No recoil-mass window is imposed: RM(phi phi) is the fit observable.
    }

alg_phiphiomega.with_decay_card(decay_card_phiphiomega).apply(event_selection)
alg_phiphiphi.with_decay_card(decay_card_phiphiphi).apply(event_selection)

### ------------------------------- Execute ------------------------------ ###
all_datasets = data_points + incMC_points + all_exMC
root_files_phiphiomega = alg_phiphiomega.execute_on(all_datasets)
root_files_phiphiphi    = alg_phiphiphi.execute_on(all_datasets)