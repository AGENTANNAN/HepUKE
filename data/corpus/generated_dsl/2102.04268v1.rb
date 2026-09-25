# ============================================================================
# Dataset preparation
# ============================================================================
# 17 BESIII scan points spanning 3.7730 - 4.5995 GeV
# (BESIII sample-name convention: [BOSS_version]_[ECMS in MeV])
scan_names = %w[
  712_3773 703_4009 703_4180 703_4190 703_4200 703_4210 703_4220
  703_4230 703_4237 703_4246 703_4260 703_4270 703_4280 703_4360
  703_4420 703_4530 703_4600
]

scan_data  = scan_names.map { |n| DatasetManager.real_data.find(n) }     # real data at each scan point
scan_incMC = scan_names.map { |n| DatasetManager.inclusive_mc.find(n) }  # matching inclusive MC at each energy

# --- Decay card: e+e- -> p pbar eta  (ConExc continuum mode 51, eta -> gamma gamma) ---
decay_card_eta = <<~DECAYCARD
    Decay vpho
    1.000 p+ anti-p- eta ConExc 51;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# --- Decay card: e+e- -> p pbar omega  (KKMC + psi(4260), omega -> pi+ pi- pi0, pi0 -> gamma gamma) ---
decay_card_omega = <<~DECAYCARD
    Decay psi(4260)
    1.000 p+ anti-p- omega PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# --- Exclusive MC: 100k events at every scan point for both modes ---
exMC_eta = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name   = "exmc_ppbar_eta"
  config.events        = 100_000
  config.decay_card    = decay_card_eta
  config.cross_section = :default
end

exMC_omega = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name   = "exmc_ppbar_omega"
  config.events        = 100_000
  config.decay_card    = decay_card_omega
  config.cross_section = :default
end

# ============================================================================
# Event selection (BOSS)
# ============================================================================
# Shared vocabulary for both channels:
#   charged tracks : |cos(theta)| < 0.93, |Vz| < 10 cm, Vr < 1 cm
#   photons        : E > 25 MeV (barrel) / 50 MeV (endcap), EMC timing 0-14,
#                    > 10 deg away from any track, at least two photons
#   PID            : probability method, prob_cut 0.001

# --------------------------- Channel 1: e+e- -> p pbar eta ---------------------------
alg_name_eta = "PpbarEta"
alg_eta = Algorithm.new(alg_name_eta)
alg_eta.set_header(["#{alg_name_eta}Alg/#{alg_name_eta}.h"])
       .set_constant({"ECMS" => [:double, 4.260]})   # nominal scan energy (varies point to point)
       .set_alias({"std::vector<double>" => "Vdouble"})

sel_eta = Selection.new
sel_eta.select_track {                                  # charged track selection
          cos_theta 0.93        # |cos(theta)| < 0.93
          Vz        10.0        # |Vz| < 10 cm
          Vr        1.0         # Vr < 1 cm
          nChrp     "==1"       # one positive track
          nChrn     "==1"       # one negative track
          nNet      "==0"       # net charge zero
        }
        .select_photon {                                # photon selection
          tdc_emc_start     0
          tdc_emc_end       14
          angle_to_track    10.0
          energyThreshold_b 0.025   # 25 MeV (barrel)
          energyThreshold_e 0.050   # 50 MeV (endcap)
          nGam              ">=2"   # at least two photons
        }
        .pid(method: :probability) {                    # particle identification
          prob_cut 0.001
          identify :proton, against: [:kaon, :pion]     # p / pbar, separated from kaons and pions
          nprp "==1"
          nprm "==1"
        }
        # Reconstruct eta from two photons (1-C mass constraint, chi2 < 25)
        .kalman_kinematic_fit([:gamma, :gamma]) {
          invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
          chi2_cut 25
          neta ">=1"
        }
        # 4C kinematic fit to p pbar eta (nominal masses)
        .kinematic_fit([:prp, :prm, :eta]) {
          nominal
          constrain_four_momentum
          chi2_cut 200          # loose BOSS-level cut; tight range applied in ROOT
        }

alg_eta.with_decay_card(decay_card_eta).apply(sel_eta)

# --------------------------- Channel 2: e+e- -> p pbar omega ---------------------------
alg_name_omega = "PpbarOmega"
alg_omega = Algorithm.new(alg_name_omega)
alg_omega.set_header(["#{alg_name_omega}Alg/#{alg_name_omega}.h"])
         .set_constant({"ECMS" => [:double, 4.260]})   # nominal scan energy (varies point to point)
         .set_alias({"std::vector<double>" => "Vdouble"})

sel_omega = Selection.new
sel_omega.select_track {                                # charged track selection
            cos_theta 0.93      # |cos(theta)| < 0.93
            Vz        10.0      # |Vz| < 10 cm
            Vr        1.0       # Vr < 1 cm
            nChrp     "==2"     # two positive tracks
            nChrn     "==2"     # two negative tracks
            nNet      "==0"     # net charge zero
          }
          .select_photon {                              # photon selection
            tdc_emc_start     0
            tdc_emc_end       14
            angle_to_track    10.0
            energyThreshold_b 0.025
            energyThreshold_e 0.050
            nGam              ">=2"
          }
          .pid(method: :probability) {                  # particle identification
            prob_cut 0.001
            identify :proton, against: [:kaon, :pion]   # p / pbar
            identify :pion,   against: [:kaon, :proton] # pi+ / pi-
            nprp "==1"
            nprm "==1"
            npip "==1"
            npim "==1"
          }
          # Reconstruct pi0 from two photons (1-C mass constraint, chi2 < 25)
          .kalman_kinematic_fit([:gamma, :gamma]) {
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
            chi2_cut 25
            npi0 ">=1"
          }
          # 4C kinematic fit to p pbar pi+ pi- pi0 (nominal masses)
          .kinematic_fit([:prp, :prm, :pip, :pim, :pi0]) {
            nominal
            constrain_four_momentum
            invariant_mass_of(:pip, :pim, :pi0).within(0.752, 0.813)  # |M(pi+pi-pi0) - m_omega| < 30 MeV/c^2
            chi2_cut 200
          }

alg_omega.note(:background_veto,
  "K_S0 -> pi+ pi- and Lambda -> p pi- peaking backgrounds vetoed in the p pbar omega channel; " \
  "the veto windows (M(pi+pi-) around the K_S0 mass and M(p pi-) around the Lambda mass) are applied " \
  "in the ROOT signal-extraction fit, together with the simultaneous fit to the M(p pbar) recoil-mass spectra")

alg_omega.with_decay_card(decay_card_omega).apply(sel_omega)

# ============================================================================
# Execute on real data, inclusive MC, and exclusive MC
# ============================================================================
root_files_eta   = alg_eta.execute_on(scan_data + scan_incMC + exMC_eta)
root_files_omega = alg_omega.execute_on(scan_data + scan_incMC + exMC_omega)