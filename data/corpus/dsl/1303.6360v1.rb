# arXiv:1303.6360v1 — Study of psi(3686) -> omega K K pi (BESIII, 1.06e8 psi(3686) events)
# Mode I   : psi(3686) -> omega K_S0 K+ pi-,  omega -> pi+ pi- pi0, K_S0 -> pi+ pi-, pi0 -> gamma gamma
# Mode II  : psi(3686) -> omega K+ K- pi0,    omega -> pi+ pi- pi0, pi0 -> gamma gamma
# Each independent decay mode gets its own Algorithm object (Rule T1).

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")     # psi(3686) data
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # 1.06e8 inclusive MC events

# ---------------- Mode I: psi(3686) -> omega K_S0 K+ pi- ----------------
decay_card_modeI = <<~DECAYCARD
    Decay psi(3686)
    1.0000 omega K_S0 K+ pi-    PHSP;
    Enddecay

    Decay omega
    1.0000 pi+ pi- pi0    PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi-    PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma    PHSP;
    Enddecay

    End
DECAYCARD

# ---------------- Mode II: psi(3686) -> omega K+ K- pi0 ----------------
decay_card_modeII = <<~DECAYCARD
    Decay psi(3686)
    1.0000 omega K+ K- pi0    PHSP;
    Enddecay

    Decay omega
    1.0000 pi+ pi- pi0    PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma    PHSP;
    Enddecay

    End
DECAYCARD

exMC_modeI = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_3686_omegaKSKpi"
    config.related_dataset = psip_data
    config.events          = 100000
    config.decay_card      = decay_card_modeI
    config.cross_section   = :default
end

exMC_modeII = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_3686_omegaKKpi0"
    config.related_dataset = psip_data
    config.events          = 100000
    config.decay_card      = decay_card_modeII
    config.cross_section   = :default
end

### Event selection (BOSS) ###

# ============================ Mode I ============================
alg_name_I = "OmegaKSKPi"
alg_modeI = Algorithm.new(alg_name_I)
alg_modeI.set_header(["#{alg_name_I}Alg/#{alg_name_I}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})
         .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeI = Selection.new
sel_modeI.select_track {
              cos_theta       0.93    # |cos(theta)| < 0.93 (MDC fiducial volume)
              Vz              20.0    # |Vz| < 20 cm along the beam direction
              Vr              2.0     # Vr < 2 cm in the radial direction
              nChrp           ">=2"   # charged tracks with net charge zero
              nChrn           ">=2"
              nNet            "==0"
          }
         .pid(method: :probability) {
              prob_cut   0.001
              # A track is a kaon if Prob(K) > Prob(pi) and > Prob(p);
              # at least one track must be positively identified as a kaon.
              identify :kaon, against: [:pion, :proton]
              nkp   ">=1"
         }
         .remove([:kp <= :chrgp, :km <= :chrgn])   # keep kaons out of the pion lists
         .assign({:chrgp => :pip, :chrgn => :pim}) # remaining tracks treated as pions
         # K_S0 from a secondary vertex fit to an oppositely charged track pair;
         # the combination closest to the nominal K_S0 mass is kept.
         .secondary_vertex_fit([:pip, :pim]) {
              build_virtual_particle(:K_S0).by_minimizing_mass_difference
              remove_used_particle_from_candidate_list
         }
         .select_photon {
              tdc_emc_start     0
              tdc_emc_end       14      # 0 <= t <= 700 ns coincidence with the collision
              energyThreshold_b 0.025   # E > 25 MeV in the barrel (|cos(theta)| < 0.80)
              energyThreshold_e 0.050   # E > 50 MeV in the endcap (0.86 < |cos(theta)| < 0.92)
              angle_to_track    10.0
              nGam              ">=2"
         }
         # pi0 reconstructed from a photon pair with the nominal pi0 mass constraint
         .kalman_kinematic_fit([:gamma, :gamma]) {
              invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
              chi2_cut 25
              npi0  ">=1"
         }
         # 4C kinematic fit; when more than two photons are present the combination
         # with the smallest chi^2 of the fit is retained.
         .kinematic_fit([:K_S0, :kp, :pip, :pim, :pim, :pi0]) {
              nominal
              constrain_four_momentum   # 4C fit to the psi(3686) four-momentum
              chi2_cut 40
         }

alg_modeI
  .note(:background_veto, "psi(3686) -> pi+pi- J/psi background vetoed by requiring
    |M_recoil(pi+pi-) - m_J/psi| > 0.007 GeV/c^2 (recoil mass not expressible in the
    BOSS selection; applied in the ROOT analysis)")
  .note(:background_veto, "K_S0 and omega sideband regions used for background
    estimation: 0.012 < |M(pi+pi-) - m_K_S0| < 0.020 GeV/c^2 and
    0.06 < |M(pi+pi-pi0) - m_omega| < 0.10 GeV/c^2; signal windows
    0.489 < M(pi+pi-) < 0.505 GeV/c^2 and 0.743 < M(pi+pi-pi0) < 0.823 GeV/c^2")
  .note(:efficiency_curve, "background estimated from the K_S0 and omega sidebands,
    normalized by the ratio of MC events in the sideband to the signal region")
  .with_decay_card(decay_card_modeI)
  .apply(sel_modeI)

# ============================ Mode II ============================
alg_name_II = "OmegaKKPi0"
alg_modeII = Algorithm.new(alg_name_II)
alg_modeII.set_header(["#{alg_name_II}Alg/#{alg_name_II}.h"])
          .set_constant({"ECMS" => [:double, 3.686]})
          .set_alias({"std::vector<double>" => "Vdouble"})

sel_modeII = Selection.new
sel_modeII.select_track {
              cos_theta       0.93    # |cos(theta)| < 0.93
              Vz              20.0
              Vr              2.0
              nChrp           "==2"   # four charged tracks with zero net charge
              nChrn           "==2"
              nNet            "==0"
          }
          .pid(method: :probability) {
              prob_cut   0.001
              # At least two charged tracks identified as kaons; when more than two
              # are identified the two oppositely charged tracks with the largest
              # Prob(K) are chosen.
              identify :kaon, against: [:pion, :proton]
              nkp   ">=1"
              nkm   ">=1"
          }
          .remove([:kp <= :chrgp, :km <= :chrgn])
          .assign({:chrgp => :pip, :chrgn => :pim})
          .select_photon {
              tdc_emc_start     0
              tdc_emc_end       14
              energyThreshold_b 0.025
              energyThreshold_e 0.050
              angle_to_track    10.0
              nGam              ">=4"   # two pi0 -> four photons
          }
          # Competing kinematic-fit hypotheses used to reject the 3-photon and
          # 5-photon backgrounds. No chi2_cut and no nominal: only the chi^2 values
          # are stored, and the comparison is applied in the ROOT analysis.
          .kinematic_fit([:kp, :km, :pip, :pim, :gamma, :gamma, :gamma]) {
              constrain_four_momentum
          }
          .kinematic_fit([:kp, :km, :pip, :pim, :gamma, :gamma, :gamma, :gamma, :gamma]) {
              constrain_four_momentum
          }
          # Two pi0 candidates reconstructed from photon pairs (mass constrained)
          .kalman_kinematic_fit([:gamma, :gamma]) {
              invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
              chi2_cut 25
              npi0  ">=2"
          }
          # 6C kinematic fit: four-momentum conservation plus the two pi0 mass
          # constraints; the combination with the least chi^2 is chosen.
          .kinematic_fit([:kp, :km, :pip, :pim, :pi0, :pi0]) {
              nominal
              constrain_four_momentum
              chi2_cut 100
          }

alg_modeII
  .note(:background_veto, "the nominal 4-gamma hypothesis is required to have a
    smaller chi^2_4C than the K+K-pi+pi-3gamma and K+K-pi+pi-5gamma hypotheses
    (the two competing fits above store their chi^2 values for the ROOT-level veto)")
  .note(:background_veto, "psi(3686) -> pi+pi- J/psi, pi0pi0 J/psi and gamma gamma J/psi
    backgrounds removed with |M_recoil(pi+pi-) - m_J/psi| > 0.007 GeV/c^2,
    |M_recoil(pi0pi0) - m_J/psi| > 0.06 GeV/c^2 and
    |M_recoil(gamma gamma) - m_J/psi| > 0.05 GeV/c^2 (applied in the ROOT analysis)")
  .note(:background_veto, "omega signal window 0.743 < M(pi+pi-pi0) < 0.823 GeV/c^2;
    non-omega background estimated from the omega sideband
    0.06 < |M(pi+pi-pi0) - m_omega| < 0.10 GeV/c^2")
  .with_decay_card(decay_card_modeII)
  .apply(sel_modeII)

### Job submission ###
root_files_I  = alg_modeI.execute_on([psip_data, psip_incMC, exMC_modeI])
root_files_II = alg_modeII.execute_on([psip_data, psip_incMC, exMC_modeII])
