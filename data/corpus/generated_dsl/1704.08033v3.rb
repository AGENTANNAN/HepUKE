# ============================================================================
# BOSS DSL  —  e+e- -> eta h_c,  h_c -> gamma eta_c,  eta -> gamma gamma,
#              eta_c -> X_i  (16 exclusive hadronic modes)
# 703 scan: 15 c.m. energy points, 4.085 - 4.600 GeV
# ============================================================================

### Dataset description ###
# 15 c.m. energy points of the 703 scan (4.085 --> 4.600 GeV)
energy_point_names = %w[
  703_4090 703_4180 703_4190 703_4200 703_4210
  703_4220 703_4230 703_4245 703_4260 703_4310
  703_4360 703_4420 703_4470 703_4530 703_4600
]
real_data_points = energy_point_names.map { |n| DatasetManager.real_data.find(n) }

# Inclusive MC samples at 4.23, 4.26 and 4.36 GeV
inclusive_mc_points = %w[703_4230 703_4260 703_4360].map do |n|
  DatasetManager.inclusive_mc.find(n)
end

# ---- auxiliary decay-card fragments ----
extra_pi0 = "Decay pi0\n1.0000 gamma gamma PHSP;\nEnddecay\n\n"
extra_ks  = "Decay K_S0\n1.0000 pi+ pi- PHSP;\nEnddecay\n\n"

# ---- the 16 eta_c -> X_i hadronic modes ----
modes = [
  # 1. eta_c -> p pbar
  { name: "ppbar", products: "p+ anti-p-", extra: "",
    nChrp: ">=1", nChrn: ">=1",
    has_proton: true, has_kaon: false, has_pion: false,
    nprp: ">=1", nprm: ">=1",
    nGam: ">=3", nEtaReq: ">=1", nPi0Req: nil, need_ks: false,
    participants: [:gamma, :eta, :prp, :prm], vtx: [2, 3] },

  # 2. eta_c -> 2(pi+ pi-)
  { name: "2pip2pim", products: "pi+ pi- pi+ pi-", extra: "",
    nChrp: ">=2", nChrn: ">=2",
    has_proton: false, has_kaon: false, has_pion: true,
    npip: ">=2", npim: ">=2",
    nGam: ">=3", nEtaReq: ">=1", nPi0Req: nil, need_ks: false,
    participants: [:gamma, :eta, :pip, :pip, :pim, :pim], vtx: [2, 3, 4, 5] },

  # 3. eta_c -> 2(K+ K-)
  { name: "2kp2km", products: "K+ K- K+ K-", extra: "",
    nChrp: ">=2", nChrn: ">=2",
    has_proton: false, has_kaon: true, has_pion: false,
    nkp: ">=2", nkm: ">=2",
    nGam: ">=3", nEtaReq: ">=1", nPi0Req: nil, need_ks: false,
    participants: [:gamma, :eta, :kp, :kp, :km, :km], vtx: [2, 3, 4, 5] },

  # 4. eta_c -> K+ K- pi+ pi-
  { name: "kpkm2pi", products: "K+ K- pi+ pi-", extra: "",
    nChrp: ">=2", nChrn: ">=2",
    has_proton: false, has_kaon: true, has_pion: true,
    nkp: ">=1", nkm: ">=1", npip: ">=1", npim: ">=1",
    nGam: ">=3", nEtaReq: ">=1", nPi0Req: nil, need_ks: false,
    participants: [:gamma, :eta, :kp, :km, :pip, :pim], vtx: [2, 3, 4, 5] },

  # 5. eta_c -> p pbar pi+ pi-
  { name: "ppbar2pi", products: "p+ anti-p- pi+ pi-", extra: "",
    nChrp: ">=2", nChrn: ">=2",
    has_proton: true, has_kaon: false, has_pion: true,
    nprp: ">=1", nprm: ">=1", npip: ">=1", npim: ">=1",
    nGam: ">=3", nEtaReq: ">=1", nPi0Req: nil, need_ks: false,
    participants: [:gamma, :eta, :prp, :prm, :pip, :pim], vtx: [2, 3, 4, 5] },

  # 6. eta_c -> 3(pi+ pi-)
  { name: "3pip3pim", products: "pi+ pi- pi+ pi- pi+ pi-", extra: "",
    nChrp: ">=3", nChrn: ">=3",
    has_proton: false, has_kaon: false, has_pion: true,
    npip: ">=3", npim: ">=3",
    nGam: ">=3", nEtaReq: ">=1", nPi0Req: nil, need_ks: false,
    participants: [:gamma, :eta, :pip, :pip, :pip, :pim, :pim, :pim], vtx: [2, 3, 4, 5, 6, 7] },

  # 7. eta_c -> K+ K- 2(pi+ pi-)
  { name: "kpkm2pip2pim", products: "K+ K- pi+ pi- pi+ pi-", extra: "",
    nChrp: ">=3", nChrn: ">=3",
    has_proton: false, has_kaon: true, has_pion: true,
    nkp: ">=1", nkm: ">=1", npip: ">=2", npim: ">=2",
    nGam: ">=3", nEtaReq: ">=1", nPi0Req: nil, need_ks: false,
    participants: [:gamma, :eta, :kp, :km, :pip, :pip, :pim, :pim], vtx: [2, 3, 4, 5, 6, 7] },

  # 8. eta_c -> K+ K- pi0
  { name: "kpkm_pi0", products: "K+ K- pi0", extra: extra_pi0,
    nChrp: ">=1", nChrn: ">=1",
    has_proton: false, has_kaon: true, has_pion: false,
    nkp: ">=1", nkm: ">=1",
    nGam: ">=5", nEtaReq: ">=1", nPi0Req: ">=1", need_ks: false,
    participants: [:gamma, :eta, :kp, :km, :pi0], vtx: [2, 3] },

  # 9. eta_c -> p pbar pi0
  { name: "ppbar_pi0", products: "p+ anti-p- pi0", extra: extra_pi0,
    nChrp: ">=1", nChrn: ">=1",
    has_proton: true, has_kaon: false, has_pion: false,
    nprp: ">=1", nprm: ">=1",
    nGam: ">=5", nEtaReq: ">=1", nPi0Req: ">=1", need_ks: false,
    participants: [:gamma, :eta, :prp, :prm, :pi0], vtx: [2, 3] },

  # 10. eta_c -> K_S0 K+ pi-
  { name: "Ks_kpi", products: "K_S0 K+ pi-", extra: extra_ks,
    nChrp: ">=2", nChrn: ">=2",
    has_proton: false, has_kaon: true, has_pion: true,
    nkp: ">=1", npip: ">=1", npim: ">=2",
    nGam: ">=3", nEtaReq: ">=1", nPi0Req: nil, need_ks: true,
    participants: [:gamma, :eta, :K_S0, :kp, :pim], vtx: [3, 4] },

  # 11. eta_c -> K_S0 K+ pi- pi+ pi-
  { name: "Ks_kpipipi", products: "K_S0 K+ pi- pi+ pi-", extra: extra_ks,
    nChrp: ">=3", nChrn: ">=3",
    has_proton: false, has_kaon: true, has_pion: true,
    nkp: ">=1", npip: ">=2", npim: ">=3",
    nGam: ">=3", nEtaReq: ">=1", nPi0Req: nil, need_ks: true,
    participants: [:gamma, :eta, :K_S0, :kp, :pip, :pim], vtx: [3, 4, 5] },

  # 12. eta_c -> pi+ pi- eta
  { name: "2pipi_eta", products: "pi+ pi- eta", extra: "",
    nChrp: ">=1", nChrn: ">=1",
    has_proton: false, has_kaon: false, has_pion: true,
    npip: ">=1", npim: ">=1",
    nGam: ">=5", nEtaReq: ">=2", nPi0Req: nil, need_ks: false,
    participants: [:gamma, :eta, :pip, :pim, :eta], vtx: [2, 3] },

  # 13. eta_c -> K+ K- eta
  { name: "kpkm_eta", products: "K+ K- eta", extra: "",
    nChrp: ">=1", nChrn: ">=1",
    has_proton: false, has_kaon: true, has_pion: false,
    nkp: ">=1", nkm: ">=1",
    nGam: ">=5", nEtaReq: ">=2", nPi0Req: nil, need_ks: false,
    participants: [:gamma, :eta, :kp, :km, :eta], vtx: [2, 3] },

  # 14. eta_c -> 2(pi+ pi-) eta
  { name: "2pipi2pim_eta", products: "pi+ pi- pi+ pi- eta", extra: "",
    nChrp: ">=2", nChrn: ">=2",
    has_proton: false, has_kaon: false, has_pion: true,
    npip: ">=2", npim: ">=2",
    nGam: ">=5", nEtaReq: ">=2", nPi0Req: nil, need_ks: false,
    participants: [:gamma, :eta, :pip, :pip, :pim, :pim, :eta], vtx: [2, 3, 4, 5] },

  # 15. eta_c -> pi+ pi- pi0 pi0
  { name: "2pipi_2pi0", products: "pi+ pi- pi0 pi0", extra: extra_pi0,
    nChrp: ">=1", nChrn: ">=1",
    has_proton: false, has_kaon: false, has_pion: true,
    npip: ">=1", npim: ">=1",
    nGam: ">=7", nEtaReq: ">=1", nPi0Req: ">=2", need_ks: false,
    participants: [:gamma, :eta, :pip, :pim, :pi0, :pi0], vtx: [2, 3] },

  # 16. eta_c -> 2(pi+ pi-) pi0 pi0
  { name: "2pipi2pim_2pi0", products: "pi+ pi- pi+ pi- pi0 pi0", extra: extra_pi0,
    nChrp: ">=2", nChrn: ">=2",
    has_proton: false, has_kaon: false, has_pion: true,
    npip: ">=2", npim: ">=2",
    nGam: ">=7", nEtaReq: ">=1", nPi0Req: ">=2", need_ks: false,
    participants: [:gamma, :eta, :pip, :pip, :pim, :pim, :pi0, :pi0], vtx: [2, 3, 4, 5] }
]

### Signal decay cards / exclusive MC / event selection (one Algorithm per mode) ###
modes.each do |mode|
  # ---- decay card for e+e- -> eta h_c, h_c -> gamma eta_c, eta -> gamma gamma, eta_c -> X_i ----
  decay_card = <<~DECAYCARD
    Decay psi(4260)
    1.0000 eta h_c PHSP;
    Enddecay

    Decay h_c
    1.0000 gamma eta_c PHSP;
    Enddecay

    Decay eta
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay eta_c
    1.0000 #{mode[:products]} PHSP;
    Enddecay

#{mode[:extra]}End
  DECAYCARD

  # ---- 40k-event exclusive MC at each of the 15 energy points ----
  exMC_mode = DatasetManager.create_exclusive_mc_for(real_data_points) do |config|
    config.sample_name   = "exmc_eta_hc_etac_#{mode[:name]}"
    config.events        = 40_000
    config.decay_card    = decay_card
    config.cross_section = :default
  end

  # ---- Algorithm ----
  alg_name = "EtaHc_#{mode[:name]}"
  alg = Algorithm.new(alg_name)
  alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
     .set_constant({ "ECMS" => [:double, 4.260] })
     .note(:ecms_per_energy_point,
           "the 15 scan points span 4.085-4.600 GeV; the 4C fit must use the measured " \
           "beam energy of the current run (per-run CMS energy). ECMS=4.260 GeV is only " \
           "a nominal default for this spec.")

  # ---- event selection ----
  sel = Selection.new
  sel.select_track {          # charged tracks (mode-specific minimum multiplicities)
        cos_theta 0.93        # |cos(theta)| < 0.93
        Vz 10.0               # |Vz| < 10 cm
        Vr 1.0                # Vr < 1 cm
        nChrp mode[:nChrp]
        nChrn mode[:nChrn]
      }
     .select_photon {         # photons
        tdc_emc_start 0
        tdc_emc_end 14
        angle_to_track 10.0
        energyThreshold_b 0.025   # > 25 MeV (barrel)
        energyThreshold_e 0.050   # > 50 MeV (endcap)
        nGam mode[:nGam]
      }
     .pid(method: :probability) {
        prob_cut 0.001
        # identify the hadron species present in this mode against the other hypotheses
        identify :proton, against: [:kaon, :pion] if mode[:has_proton]
        identify :kaon,   against: [:pion, :proton] if mode[:has_kaon]
        identify :pion,   against: [:kaon, :proton] if mode[:has_pion]
        # corresponding minimum counts
        nprp mode[:nprp] if mode[:nprp]
        nprm mode[:nprm] if mode[:nprm]
        nkp  mode[:nkp]  if mode[:nkp]
        nkm  mode[:nkm]  if mode[:nkm]
        npip mode[:npip] if mode[:npip]
        npim mode[:npim] if mode[:npim]
      }

  # ---- 1C mass-constrained fit of gamma-gamma pairs to the eta mass (top eta), chi2 < 25 ----
  sel = sel.kalman_kinematic_fit([:gamma, :gamma]) {
          invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
          chi2_cut 25
          neta mode[:nEtaReq]
        }

  # additional 1C fit of gamma-gamma pairs to the pi0 mass where a pi0 appears in X_i
  if mode[:nPi0Req]
    sel = sel.kalman_kinematic_fit([:gamma, :gamma]) {
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
            chi2_cut 25
            npi0 mode[:nPi0Req]
          }
  end

  # ---- K_S0 -> pi+ pi- secondary-vertex fit for modes containing a K_S0 ----
  if mode[:need_ks]
    sel = sel.secondary_vertex_fit([:pip, :pim]) {
            build_virtual_particle(:K_S0).by_minimizing_mass_difference
            remove_used_particle_from_candidate_list
          }
  end

  # ---- 4C kinematic fit to gamma eta + X_i ----
  # common production-vertex fit on the prompt charged tracks, 4-momentum
  # constraint, nominal, loose chi2 < 200 (tight chi2 < 25, the E_gamma window
  # and the recoil / hadronic mass windows are applied at the ROOT level).
  sel = sel.kinematic_fit(mode[:participants]) {
          nominal
          vertex_fit(mode[:vtx])
          constrain_four_momentum
          chi2_cut 200
        }

  alg.with_decay_card(decay_card).apply(sel)

  # ---- run on the 15 real-data points, the 3 inclusive-MC samples and this mode's exclusive MC ----
  root_files = alg.execute_on(real_data_points + inclusive_mc_points + exMC_mode)
end