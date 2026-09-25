### Dataset description ###
psi3770_data  = DatasetManager.real_data.find("712_3773")     # ψ(3770) real data at √s = 3.773 GeV (2.93 fb⁻¹)
psi3770_incMC = DatasetManager.inclusive_mc.find("712_3773")  # corresponding inclusive MC at 3.773 GeV

# Decay card for the full signal chain:
#   ψ(3770) → D+ D−, D+ → K0 e+ νe, K0 → K_S0 → π0 π0
decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.0000   D+  D-        PHSP;
    Enddecay

    Decay D+
    1.0000   K0  e+  nu_e  PHSP;
    Enddecay

    Decay K0
    1.0000   K_S0         PHSP;
    Enddecay

    Decay K_S0
    1.0000   pi0  pi0     PHSP;
    Enddecay

    Decay pi0
    1.0000   gamma  gamma PHSP;
    Enddecay

    End
DECAYCARD

# 500k exclusive MC events for the full ψ(3770) → D+D−, D+ → K0 e+νe, K0 → K_S0 → π0π0 chain
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_Dp_K0enu_Ks_pi0pi0"
  config.related_dataset = psi3770_data
  config.events          = 500_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Tag-based event selection (BOSS) ###
# This is a tag analysis: the D− tag comes from the pre-stored DTag candidates (six hadronic
# single-tag modes); the signal side is the semileptonic D+ → K_S0(→π0π0) e+ νe, reconstructed
# from the tracks/showers the tag did not use, with the neutrino left missing.
alg_name      = "DpToK0EnuTagDm"
my_algorithm  = TagAnalysis.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.773]})   # 4C fit against the measured lab energy

# --- Tag side: hadronic D− single tag in six modes ---
my_algorithm.tag_side(:Dm) do |t|
  t.modes :DmtoKPiPi,        # D− → K+ π− π−
          :DmtoKsPi,         # D− → K_S0 π−
          :DmtoKPiPiPi0,     # D− → K+ π− π− π0
          :DmtoKsPiPi0,      # D− → K_S0 π− π0
          :DmtoKsPiPiPi,     # D− → K_S0 π+ π− π−
          :DmtoKKPi          # D− → K+ K− π−
  t.charm(-1)                # pin the tagged side to D−
end

# --- Signal side: D+ → K0 e+ νe, K0 → K_S0 → π0 π0 ---
my_algorithm.signal_side do |s|
  s.charged ep: 1            # exactly one good charged track not used by the tag, identified as e+
  s.photons 4                # at least four good photons (two π0 → γγ)
  s.missing :nu_e            # massless neutrino of the semileptonic decay
  s.require_charge 1         # signal track carries charge opposite to the tag D−
  s.min_photon_angle 10.0    # photon opening angle to the closest charged track
  s.min_photon_energy 0.025  # EMC shower energy floor (barrel value)
end

# --- 4C kinematic fit: tag D−, four photons, e+ and the massless ν to the lab energy,
#     with both π0 masses constrained to their nominal values ---
my_algorithm.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # first π0
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # second π0
  f.chi2_cut 200
end

# --- BOSS-side procedures that have no dedicated DSL construct ---
my_algorithm
  .note(:tag_reconstruction_cuts,
        "DTagAlg tag-side reconstruction uses: charged tracks with |cosθ|<0.93 and, unless "
        "they come from K_S0, V_xy<1 cm and |V_z|<10 cm; K/π separation from combined dE/dx "
        "and TOF confidence levels (CL_K>CL_π for a kaon, CL_π>CL_K for a pion); K_S0 daughters "
        "with |V_z|<20 cm, assigned π+π− without PID, forming a common vertex with "
        "|M(π+π−)−m_K_S0|<12 MeV/c² and L/σ_L>2; photons with shower time within 700 ns of the "
        "event start, E>25 MeV (barrel)/50 MeV (endcap) and opening angle >10° to the closest "
        "charged track; π0 from γγ with M(γγ) in (0.115,0.150) GeV/c² plus a mass-constrained fit.")
  .note(:tag_deltae_window,
        "Per tag mode only the combination with minimum |ΔE| is kept, requiring |ΔE| within "
        "±25 MeV for K+π−π−, K_S0π−, K_S0π+π−π− and K+K−π− and within (−55,+40) MeV for "
        "K+π−π−π0 and K_S0π−π0; the windows are mode-dependent, so ΔE is stored and windowed "
        "in the ROOT stage rather than cut in BOSS.")
  .note(:tag_mbc_signal_region,
        "ST signal region taken as 1.863 < M_BC < 1.877 GeV/c²; M_BC is stored unconditionally "
        "by the tag DSL and the signal-region window is applied in ROOT.")
  .note(:ks_to_pi0pi0_selection,
        "Signal K0 → K_S0 → π0π0 candidate must have M(π0π0) in (0.45,0.51) GeV/c², choosing "
        "the combination with the smallest sum of the two π0 mass-constrained χ² values.")
  .note(:fsr_recovery,
        "Photons within 5° of the electron direction are added to the electron four-momentum "
        "to partially recover FSR and bremsstrahlung before the 4C fit.")
  .note(:electron_pid_selection,
        "Signal-side electron identified via combined dE/dx, TOF and EMC information with "
        "CL_e>0.001 and CL_e/(CL_e+CL_π+CL_K)>0.8; this criterion differs from the fixed "
        "SimplePIDSvc lepton thresholds exposed by the tag DSL charged-key (`ep`) recipe.")

# Render the tag specification (no Selection argument) and run
my_algorithm.with_decay_card(decay_card_signal).apply
root_files = my_algorithm.execute_on([psi3770_data, psi3770_incMC, exMC_signal])