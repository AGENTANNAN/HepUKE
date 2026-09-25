### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")      # psi(2S) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # corresponding inclusive MC

# ---------------------------------------------------------------------------
# Decay cards (EvtGen format) — one per signal mode.
# The whole chain psi(2S) -> gamma eta_c', eta_c' -> VV, V -> PP is generated
# with phase space, as described.
# ---------------------------------------------------------------------------
# Mode I: eta_c' -> rho0 rho0, rho0 -> pi+ pi-   (final state 2(pi+ pi-))
decay_card_rho = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma eta_c'        PHSP;
    Enddecay

    Decay eta_c'
    1.000 rho0 rho0           PHSP;
    Enddecay

    Decay rho0
    1.000 pi+ pi-             PHSP;
    Enddecay

    End
DECAYCARD

# Mode II: eta_c' -> K*0 anti-K*0, K*0 -> K+ pi-, anti-K*0 -> K- pi+
#          (final state pi+ pi- K+ K-)
decay_card_kstar = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma eta_c'        PHSP;
    Enddecay

    Decay eta_c'
    1.000 K*0 anti-K*0        PHSP;
    Enddecay

    Decay K*0
    1.000 K+ pi-              PHSP;
    Enddecay

    Decay anti-K*0
    1.000 K- pi+              PHSP;
    Enddecay

    End
DECAYCARD

# Mode III: eta_c' -> phi phi, phi -> K+ K-   (final state 2(K+ K-))
decay_card_phi = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma eta_c'        PHSP;
    Enddecay

    Decay eta_c'
    1.000 phi phi             PHSP;
    Enddecay

    Decay phi
    1.000 K+ K-               PHSP;
    Enddecay

    End
DECAYCARD

# ---------------------------------------------------------------------------
# Exclusive MC: 200k events for each of the three modes
# ---------------------------------------------------------------------------
exMC_rho = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_psip_gam_etacp_rho0rho0"
    config.related_dataset = psip_data
    config.events          = 200_000
    config.decay_card      = decay_card_rho
    config.cross_section   = :default
end

exMC_kstar = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_psip_gam_etacp_kstar0kstar0"
    config.related_dataset = psip_data
    config.events          = 200_000
    config.decay_card      = decay_card_kstar
    config.cross_section   = :default
end

exMC_phi = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_psip_gam_etacp_phiphi"
    config.related_dataset = psip_data
    config.events          = 200_000
    config.decay_card      = decay_card_phi
    config.cross_section   = :default
end

# ---------------------------------------------------------------------------
# Event selection (BOSS)
# The three channels share one common selection chain; only the final-state
# particle combination entering the kinematic fits differs.
# ---------------------------------------------------------------------------
common_selection = Selection.new
    .select_track {                 # charged track quality cuts
        cos_theta 0.93              # |cos(theta)| < 0.93
        Vz        10.0              # |Vz| < 10 cm
        Vr        1.0               # Vr < 1 cm
        nChrp     "==2"             # exactly two positive tracks
        nChrn     "==2"             # exactly two negative tracks
        nNet      "==0"             # net charge zero
    }
    .select_photon {                # photon selection
        tdc_emc_start      0        # TDC start time
        tdc_emc_end       14        # TDC end time
        angle_to_track    10.0      # at least 10 degrees from any charged track
        energyThreshold_b  0.025    # 25 MeV threshold in the barrel
        energyThreshold_e  0.025    # 25 MeV threshold in the endcap
        nGam               ">=1"    # at least one photon (radiative photon)
    }
    .pid(method: :chi2_sum) {       # combined chi2-sum PID over pi/K/p hypotheses
        chi_min_cut 4               # chi2_PID < 4
        identify :pion, :kaon, :proton
    }

# ---------------- Mode I: psi(2S) -> gamma eta_c', eta_c' -> rho0 rho0 ----------------
alg_name_rho = "PsipGamEtacpRhoRho"
alg_rho = Algorithm.new(alg_name_rho)
alg_rho.set_header(["#{alg_name_rho}Alg/#{alg_name_rho}.h"])
       .set_constant({"ECMS" => [:double, 3.686]})
       .set_alias({"std::vector<double>" => "Vdouble"})

sel_rho = common_selection.dup
    .kinematic_fit([:gamma, :pip, :pim, :pip, :pim]) {   # 4C fit: picks the radiative photon
        constrain_four_momentum
        chi2_cut 40
    }
    .kinematic_fit([:gamma, :pip, :pim, :pip, :pim]) {   # nominal 3C fit -> fitted M_VV
        nominal
        constrain_three_momentum
    }

alg_rho.with_decay_card(decay_card_rho).apply(sel_rho)
root_files_rho = alg_rho.execute_on([psip_data, psip_incMC, exMC_rho])

# ---------------- Mode II: psi(2S) -> gamma eta_c', eta_c' -> K*0 anti-K*0 ----------------
alg_name_kstar = "PsipGamEtacpKstar0Kstar0"
alg_kstar = Algorithm.new(alg_name_kstar)
alg_kstar.set_header(["#{alg_name_kstar}Alg/#{alg_name_kstar}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})
         .set_alias({"std::vector<double>" => "Vdouble"})

sel_kstar = common_selection.dup
    .kinematic_fit([:gamma, :kp, :pim, :km, :pip]) {     # 4C fit: picks the radiative photon
        constrain_four_momentum
        chi2_cut 40
    }
    .kinematic_fit([:gamma, :kp, :pim, :km, :pip]) {     # nominal 3C fit -> fitted M_VV
        nominal
        constrain_three_momentum
    }

alg_kstar.with_decay_card(decay_card_kstar).apply(sel_kstar)
root_files_kstar = alg_kstar.execute_on([psip_data, psip_incMC, exMC_kstar])

# ---------------- Mode III: psi(2S) -> gamma eta_c', eta_c' -> phi phi ----------------
alg_name_phi = "PsipGamEtacpPhiPhi"
alg_phi = Algorithm.new(alg_name_phi)
alg_phi.set_header(["#{alg_name_phi}Alg/#{alg_name_phi}.h"])
       .set_constant({"ECMS" => [:double, 3.686]})
       .set_alias({"std::vector<double>" => "Vdouble"})

sel_phi = common_selection.dup
    .kinematic_fit([:gamma, :kp, :km, :kp, :km]) {       # 4C fit: picks the radiative photon
        constrain_four_momentum
        chi2_cut 40
    }
    .kinematic_fit([:gamma, :kp, :km, :kp, :km]) {       # nominal 3C fit -> fitted M_VV
        nominal
        constrain_three_momentum
    }

alg_phi.with_decay_card(decay_card_phi).apply(sel_phi)
root_files_phi = alg_phi.execute_on([psip_data, psip_incMC, exMC_phi])