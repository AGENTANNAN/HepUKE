# ============================================================================
#  ψ(2S) → D_s^- π^+   and   ψ(2S) → D_s^- ρ^+
#     D_s^- → φ e^- ν̄_e ,  φ → K^+K^- ,  ρ^+ → π^+π^0 ,  π^0 → γγ
#  Inputs: ψ(2S) 3.686 GeV real data + inclusive ψ(2S) MC
#          + 100k-event exclusive MC for each of the two decay modes
# ============================================================================

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")      # ψ(2S) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # corresponding inclusive ψ(2S) MC

# --- Decay card, Mode I: ψ(2S) → D_s^- π^+, D_s^- → φ e^- ν̄_e, φ → K^+K^- ---
decay_card_ds_pi = <<~DECAYCARD
  Decay psi(2S)
  1.0000 D_s- pi+                  PHSP;
  Enddecay

  Decay D_s-
  1.0000 phi e- anti-nu_e          PHSP;
  Enddecay

  Decay phi
  1.0000 K+ K-                     VSS;
  Enddecay

  End
DECAYCARD

# --- Decay card, Mode II: ψ(2S) → D_s^- ρ^+, ρ^+ → π^+π^0, π^0 → γγ ---
decay_card_ds_rho = <<~DECAYCARD
  Decay psi(2S)
  1.0000 D_s- rho+                 PHSP;
  Enddecay

  Decay D_s-
  1.0000 phi e- anti-nu_e          PHSP;
  Enddecay

  Decay phi
  1.0000 K+ K-                     VSS;
  Enddecay

  Decay rho+
  1.0000 pi+ pi0                   VSS;
  Enddecay

  Decay pi0
  1.0000 gamma gamma               PHSP;
  Enddecay

  End
DECAYCARD

# --- 100k-event exclusive MC, one per decay mode ---
exMC_ds_pi = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_ds_pi_phi_e_nu"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_ds_pi
  config.cross_section   = :default
end

exMC_ds_rho = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_ds_rho_phi_e_nu"
  config.related_dataset = psip_data
  config.events          = 100_000
  config.decay_card      = decay_card_ds_rho
  config.cross_section   = :default
end

### Baseline event selection shared by both modes (tracks + photons + PID) ###
common_selection = Selection.new
  .select_track {                       # charged-track quality and multiplicity
      cos_theta  0.93                   # |cosθ| < 0.93
      Vr         1.0                    # Vr < 1 cm
      Vz         10.0                   # |Vz| < 10 cm
      nTot       "==4"                  # exactly four charged tracks
      nChrp      ">=1"                  # at least one positive track
      nChrn      ">=1"                  # at least one negative track
      nNet       "==0"                  # net charge zero
  }
  .select_photon {                      # good photon (EMC shower) selection
      tdc_emc_start     0               # EMC TDC window 0–14 (700 ns)
      tdc_emc_end       14
      angle_to_track    10.0            # > 10° to nearest charged track
      energyThreshold_b 0.025           # barrel  25 MeV
      energyThreshold_e 0.050           # endcap  50 MeV
  }
  .pid(method: :probability) {          # probability-method PID
      prob_cut 0.0                      # zero probability cut
      identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                     treat_as_electron_if_energy_above: 0.6
      identify :kaon, against: [:pion]              # K/π separation
      identify :pion, against: [:kaon, :electron]   # π/(K,e) separation
      nkp  ">=1"                        # at least one K^+
      nkm  ">=1"                        # at least one K^-
      npip ">=1"                        # at least one π^+
      nlm  ">=1"                        # at least one negative lepton
  }

# ============================================================================
#  Mode I :  ψ(2S) → D_s^- π^+ ,  D_s^- → φ e^- ν̄_e ,  φ → K^+K^-
# ============================================================================
alg_name_ds_pi = "DsPiPhiE"
alg_ds_pi = Algorithm.new(alg_name_ds_pi)
alg_ds_pi.set_header(["#{alg_name_ds_pi}Alg/#{alg_name_ds_pi}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})

sel_ds_pi = common_selection.dup
  .kinematic_fit([:kp, :km, :pip, :em]) {      # antineutrino treated as missing
      nominal
      miss_track_of :nu                        # ν̄_e is not reconstructed
      invariant_mass_of(:kp, :km).within(1.005, 1.035)   # φ → K^+K^- candidate window
      constrain_four_momentum
      chi2_cut 200                             # loose BOSS-level cut, tightened in ROOT
  }

alg_ds_pi
  .note(:electron_pid_selection,
        "electron candidates are further required to satisfy 0.86 < E/p < 1.03, " \
        "L(e) > 0.001 and L(e)/(L(e)+L(pi)+L(K)) > 0.8; the per-track probability " \
        "PID chain in BOSS does not expose E/p or likelihood-ratio cuts, so this " \
        "additional electron selection is applied outside the DSL chain")
  .note(:ds_recoil_mass_window,
        "D_s^- is reconstructed by recoil against the pi^+: " \
        "M_recoil(pi^+) = sqrt((P_CMS - p_pi+)^2) is required to lie in " \
        "(1.91, 2.02) GeV/c^2")
  .note(:missing_momentum_selection,
        "with the antineutrino treated as missing, require |p_miss| > 0.02 GeV/c " \
        "and |U_miss| < 0.064 GeV, where U_miss = E_miss - |p_miss| evaluated in " \
        "the CMS frame from the measured visible four-momenta")
  .with_decay_card(decay_card_ds_pi)
  .apply(sel_ds_pi)

root_files_ds_pi = alg_ds_pi.execute_on([psip_data, psip_incMC, exMC_ds_pi])

# ============================================================================
#  Mode II :  ψ(2S) → D_s^- ρ^+ ,  ρ^+ → π^+π^0 ,  π^0 → γγ
#             D_s^- → φ e^- ν̄_e ,  φ → K^+K^-
# ============================================================================
alg_name_ds_rho = "DsRhoPhiE"
alg_ds_rho = Algorithm.new(alg_name_ds_rho)
alg_ds_rho.set_header(["#{alg_name_ds_rho}Alg/#{alg_name_ds_rho}.h"])
          .set_constant({"ECMS" => [:double, 3.686]})

sel_ds_rho = common_selection.dup
  .select_photon { nGam ">=2" }                # at least two photons required in the ρ mode
  .kalman_kinematic_fit([:gamma, :gamma]) {    # π^0 → γγ with nominal-mass constraint
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 200
      npi0 ">=1"                               # at least one π^0 candidate
  }
  .kinematic_fit([:kp, :km, :pip, :pi0, :em]) {   # antineutrino treated as missing
      nominal
      miss_track_of :nu
      invariant_mass_of(:kp, :km).within(1.005, 1.035)   # φ → K^+K^- candidate window
      invariant_mass_of(:pip, :pi0).within(0.61, 0.93)   # ρ^+ → π^+π^0 candidate window
      constrain_four_momentum
      chi2_cut 200                             # loose BOSS-level cut, tightened in ROOT
  }

alg_ds_rho
  .note(:electron_pid_selection,
        "electron candidates are further required to satisfy 0.86 < E/p < 1.03, " \
        "L(e) > 0.001 and L(e)/(L(e)+L(pi)+L(K)) > 0.8; the per-track probability " \
        "PID chain in BOSS does not expose E/p or likelihood-ratio cuts, so this " \
        "additional electron selection is applied outside the DSL chain")
  .note(:ds_recoil_mass_window,
        "after treating the antineutrino as missing, the D_s^- recoil mass is " \
        "required to lie in (1.89, 2.08) GeV/c^2")
  .note(:missing_momentum_selection,
        "with the antineutrino treated as missing, require |p_miss| > 0.02 GeV/c " \
        "and |U_miss| < 0.10 GeV, where U_miss = E_miss - |p_miss| evaluated in " \
        "the CMS frame from the measured visible four-momenta")
  .note(:background_veto,
        "the energy of any extra photon in the rest frame of the event is required " \
        "to be below 0.12 GeV to suppress background with additional neutral clusters")
  .with_decay_card(decay_card_ds_rho)
  .apply(sel_ds_rho)

root_files_ds_rho = alg_ds_rho.execute_on([psip_data, psip_incMC, exMC_ds_rho])