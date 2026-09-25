# BESIII: tau lepton mass measurement from an energy scan near the tau+tau- threshold
# Four scan points with CM energies 3542.4, 3553.8, 3561.1 and 3600.2 MeV (~24 pb^-1 total).
# The first point lies below the tau pair production threshold and is used as a pure
# background control sample.
# 13 two-prong tau pair final states are analysed: ee, e mu, e pi, e K, mu mu, mu pi, mu K,
# pi K, pi pi, K K, e rho, mu rho and pi rho (neutrinos implied).
# Scope: dataset preparation + event selection (BOSS side) only.

### Dataset preparation ###
# tau-threshold scan points (BOSS name convention [version]_[CMS energy in MeV]).
tau_scan_data  = DatasetManager.real_data.where(cms_energy: {value: 3540..3610})
tau_scan_incMC = DatasetManager.inclusive_mc.where(cms_energy: {value: 3540..3610})

# --- Signal decay card: ten two-prong channels without photons ---
# Representative charge combination e mu; the other nine channels (e pi, e K, mu mu, mu pi,
# mu K, pi pi, pi K, K K, ee) differ only in the identity of the two charged daughters and
# are generated with the analogous card.
decay_card_twoprong = <<~DECAYCARD
    Decay psi(4260)
    1.000 tau+ tau- PHSP;
    Enddecay

    Decay tau+
    1.000 e+ nu_e anti-nu_tau PHSP;
    Enddecay

    Decay tau-
    1.000 mu- anti-nu_mu nu_tau PHSP;
    Enddecay

    End
DECAYCARD

# --- Signal decay card: the three X rho channels (X = e, mu, pi) ---
# The rho candidate is reconstructed from pi+ pi0, pi0 -> gamma gamma, hence two photons.
# Representative charge combination e rho.
decay_card_rho = <<~DECAYCARD
    Decay psi(4260)
    1.000 tau+ tau- PHSP;
    Enddecay

    Decay tau+
    1.000 e+ nu_e anti-nu_tau PHSP;
    Enddecay

    Decay tau-
    1.000 pi- anti-nu_tau rho+ PHSP;
    Enddecay

    Decay rho+
    1.000 pi+ pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC at every scan point (one sample per scan point, same card).
exMC_twoprong = DatasetManager.create_exclusive_mc_for(tau_scan_data) do |config|
  config.sample_name   = "tau_twoprong_exclusive_mc"    # auto-suffixed per energy point
  config.events        = 200_000
  config.decay_card    = decay_card_twoprong
  config.cross_section = :default
end

exMC_rho = DatasetManager.create_exclusive_mc_for(tau_scan_data) do |config|
  config.sample_name   = "tau_rho_exclusive_mc"
  config.events        = 200_000
  config.decay_card    = decay_card_rho
  config.cross_section = :default
end

# Inclusive tau pair MC (all tau decay channels) at each scan point, used for the
# selection efficiency eps_i = sum_j Br_j * eps_ij.
exMC_tau_inclusive = DatasetManager.create_exclusive_mc_for(tau_scan_data) do |config|
  config.sample_name   = "tau_inclusive_mc"
  config.events        = 500_000
  config.decay_card    = decay_card_twoprong
  config.cross_section = :default
end

### Event selection (BOSS) — ten photon-less two-prong channels ###
alg_name_twoprong = "TauMassTwoProng"
alg_twoprong = Algorithm.new(alg_name_twoprong)
alg_twoprong.set_header(["#{alg_name_twoprong}Alg/#{alg_name_twoprong}.h"])
# Multi-energy scan: ECMS is NOT set here, it is injected per scan point at run time.

sel_twoprong = Selection.new
sel_twoprong.select_track do            # good charged track selection
             cos_theta 0.93             # |cos(theta)| < 0.93
             Vz        10.0             # |Vz| < 10 cm
             Vr        1.0              # Vr = sqrt(Vx^2 + Vy^2) < 1 cm
             nChrp     "==1"            # exactly two good charged tracks -> one positive
             nChrn     "==1"            #                                   one negative
             nTot      "==2"            # total number of charged tracks is also required to be two
           end
           .select_photon do            # no photons are allowed in these ten channels
             tdc_emc_start     0
             tdc_emc_end       14       # 0 < t < 750 ns EMC timing window
             angle_to_track    20.0     # angle between the cluster and the nearest charged particle > 20 deg
             energyThreshold_b 0.025    # E > 25 MeV in the barrel EMC (|cos(theta)| < 0.8)
             energyThreshold_e 0.050    # E > 50 MeV in the endcap EMC (0.86 < |cos(theta)| < 0.92)
             nGam              "==0"    # no good photons in the event
           end
           .pid(method: :probability) do
             prob_cut 0.001
             # e/mu are separated by the high-momentum lepton path (E/p and MUC depth),
             # hadrons by the dE/dx + TOF probability.
             identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.2,
                                             treat_as_electron_if_energy_above: 0.6
             identify :pion, against: [:kaon]
             identify :kaon, against: [:pion]
           end
           # Only the two charged tracks are detected; the tau decay neutrinos are undetected.
           # The tau decay recIDs are, in card order:
           #   0 psi(4260), 1 tau+, 2 e+, 3 nu_e, 4 anti-nu_tau, 5 tau-, 6 mu-, 7 anti-nu_mu, 8 nu_tau
           .partial_miss([3, 4, 7, 8]) do
             require_recoil_mass 1.60, 1.90    # recoil mass consistent with the tau pair system
           end

alg_twoprong.note(:channel_assignment, "The ten photon-less two-prong channels (ee, e mu, e pi, e K, mu mu, mu pi, mu K, pi pi, pi K, K K) share the same final-state topology (exactly two good charged tracks, no photons) and the same selection; the channel of an event is determined offline by the PID hypothesis of each of the two tracks. The listed ALGORITHM is applied per channel with the corresponding tau decay card.")
             .note(:pid_table, "Per-track PID uses p, E, E/cp, the time-of-flight value, the MUC depth D and the number of MUC hits Nh together (Table III of the paper). The momentum windows p_min < p < p_max are channel and scan-point dependent (determined from signal MC): at the first scan point p_min (p_max) = 0.20 (0.92) GeV/c for e, 0.20 (0.90) GeV/c for mu, 0.84 (0.93) GeV/c for pi and 0.76 (0.88) GeV/c for K; for pi from rho the momentum requirement is removed.")
             .note(:ptem, "PTEM = |c p1 + c p2|_T / (W - |c p1| - |c p2|), the ratio of the net observed transverse momentum to the maximum possible missing energy, is used to reject two-photon QED background e+e- -> e+(e-e+)e-. The retained windows are channel dependent: ee PTEM > 0.3, e mu PTEM > 0.1, e pi PTEM > 0.1 (Table IV).")
             .note(:acoplanarity, "The acoplanarity angle theta_acop between the two final-state charged tracks (angle between their transverse momentum vectors) is required to be ee > 10 deg, e mu < 160 deg, e pi < 170 deg, e K < 170 deg, mu mu < 140 deg, mu h < 140 deg, h h < 160 deg (Table IV). The nominal windows were obtained from the below-threshold first scan point data, which contains only background.")
             .note(:ptem_acop_implementation, "PTEM and theta_acop are event-level variables built from both charged tracks; they cannot be expressed with the per-candidate DSL primitives and are applied as offline selections on the NTuple.")

alg_twoprong.with_decay_card(decay_card_twoprong).apply(sel_twoprong)

### Event selection (BOSS) — the three X rho channels (e rho, mu rho, pi rho) ###
alg_name_rho = "TauMassRho"
alg_rho = Algorithm.new(alg_name_rho)
alg_rho.set_header(["#{alg_name_rho}Alg/#{alg_name_rho}.h"])

sel_rho = Selection.new
sel_rho.select_track do
          cos_theta 0.93
          Vz        10.0
          Vr        1.0
          nChrp     "==1"
          nChrn     "==1"
          nTot      "==2"
        end
        .select_photon do              # exactly two photons from pi0 -> gamma gamma
          tdc_emc_start     0
          tdc_emc_end       14
          angle_to_track    20.0
          energyThreshold_b 0.025
          energyThreshold_e 0.050
          nGam              "==2"
        end
        .pid(method: :probability) do
          prob_cut 0.001
          identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.2,
                                          treat_as_electron_if_energy_above: 0.6
          identify :pion, against: [:kaon]
        end
        .kalman_kinematic_fit([:gamma, :gamma]) do   # pi0 -> gamma gamma
          invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
          chi2_cut 25
          npi0 ">=1"
        end
        # The rho candidate is built from the pi candidate and the pi0; neutrinos are undetected.
        # Tau decay recIDs in card order:
        #   0 psi(4260), 1 tau+, 2 e+, 3 nu_e, 4 anti-nu_tau,
        #   5 tau-, 6 pi-, 7 anti-nu_tau, 8 rho+, 9 pi+, 10 pi0, 11 gamma, 12 gamma
        .partial_miss([3, 4, 7]) do
          require_recoil_mass 1.55, 1.90
        end

alg_rho.note(:channel_assignment, "The three X rho channels (e rho, mu rho, pi rho) share the same final-state topology (exactly two good charged tracks and exactly two photons) and the same selection; the channel is determined offline from the PID hypothesis of the non-pion track (e, mu or pi) accompanying the rho candidate.")
       .note(:pi0_rho_windows, "The two photons must satisfy the pi0 mass window 112.8 MeV/c^2 < M(gamma gamma) < 146.4 MeV/c^2. They are then combined with a charged pi candidate to form the rho, and the rho candidate must satisfy 376.5 MeV/c^2 < M(pi pi0) < 1195.5 MeV/c^2 and p_min^rho < |p(rho)| < p_max^rho, where the momentum limits are determined from the rho momentum distribution in the signal MC at each scan point.")
       .note(:pid_table, "Per-track PID uses p, E, E/cp, the time-of-flight value, the MUC depth D and the number of MUC hits Nh (Table III). Momentum windows p_min < p < p_max are scan-point dependent; the pi from the rho has no momentum requirement.")
       .note(:ptem_acoplanarity, "Channel-dependent PTEM and acoplanarity windows are applied offline: e rho theta_acop < 170 deg, mu rho theta_acop < 150 deg (Table IV).")
       .note(:ptem_acop_implementation, "PTEM and theta_acop are event-level variables built from both charged tracks; they cannot be expressed with the per-candidate DSL primitives and are applied as offline selections on the NTuple.")

alg_rho.with_decay_card(decay_card_rho).apply(sel_rho)

### Shared BOSS-side notes ###
[alg_twoprong, alg_rho].each do |alg|
  alg.note(:scan_point_energy, "The scan point CM energies (3542.4, 3553.8, 3561.1 and 3600.2 MeV) are obtained from the beam energy measurement system (BEMS) through the Compton back-scattering edge, with a systematic uncertainty of 2e-5 on the beam energy; the CM energy of a point is E_CM = 2 sqrt(Ebar_e- * Ebar_e+) cos(theta_ee/2) with a crossing angle theta_ee = 0.022 rad.")
  alg.note(:beam_energy_spread, "The total beam energy spread of a scan point is delta_w = sqrt(delta_e-^2 + delta_e+^2) from the BEMS, determined to be (1.469 +- 0.064) MeV in the tau scan region; it enters the convolution of the tau pair cross section used in the mass fit.")
  alg.note(:luminosity, "The luminosity at each scan point is determined from two-gamma events e+e- -> gamma gamma (gamma) using the Babayaga 3.5 generator; Bhabha events provide a cross check and consistent luminosities within 2%.")
  alg.note(:tau_mass_fit, "The tau mass is extracted from a maximum likelihood fit to the CM energy dependence of the tau pair production cross section over the four scan points; the first point (below threshold) constrains the background. This fit is performed on the ROOT-level NTuple and lies outside the BOSS event selection.")
end

### Execution ###
alg_twoprong.execute_on(tau_scan_data + tau_scan_incMC + exMC_twoprong + exMC_tau_inclusive)
alg_rho.execute_on(tau_scan_data + tau_scan_incMC + exMC_rho)
