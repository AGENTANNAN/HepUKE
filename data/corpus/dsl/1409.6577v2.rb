# BESIII: e+e- -> pi0 pi0 hc, hc -> gamma eta_c, eta_c -> 16 hadronic modes Xi
# Data at sqrt(s) = 4.23, 4.26, 4.36 GeV. Neutral Zc(4020)^0 search in pi0 hc.

### Dataset preparation ###
data_4230 = DatasetManager.real_data.find("703_4230")
data_4260 = DatasetManager.real_data.find("703_4260")
data_4360 = DatasetManager.real_data.find("703_4360")
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230")
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")
incMC_4360 = DatasetManager.inclusive_mc.find("703_4360")

# Common hc -> gamma eta_c template; eta_c decay to Xi varies by mode
def decay_card_for_mode(eta_c_decay_line)
  <<~DECAYCARD
    Decay psi(4260)
    1.0000 pi0 pi0 h_c            PHSP;
    Enddecay

    Decay h_c
    1.0000 gamma eta_c             PHSP;
    Enddecay

    Decay eta_c
    #{eta_c_decay_line}
    Enddecay

    Decay pi0
    1.000 gamma gamma              PHSP;
    Enddecay

    End
  DECAYCARD
end

# 16 eta_c hadronic final states
eta_c_modes = {
  "ppbar"             => "1.000 p+ anti-p-                                     PHSP;",
  "pipiKK"            => "1.000 pi+ pi- K+ K-                                  PHSP;",
  "pipippbar"         => "1.000 pi+ pi- p+ anti-p-                             PHSP;",
  "2KK"               => "1.000 K+ K- K+ K-                                    PHSP;",
  "2pipi"             => "1.000 pi+ pi- pi+ pi-                                PHSP;",
  "3pipi"             => "1.000 pi+ pi- pi+ pi- pi+ pi-                        PHSP;",
  "2pipiKK"           => "1.000 pi+ pi- pi+ pi- K+ K-                          PHSP;",
  "KsKpi"             => "1.000 K_S0 K+ pi-                                    PHSP;",
  "KsKpipipi"         => "1.000 K_S0 K+ pi- pi+ pi-                            PHSP;",
  "KKpi0"             => "1.000 K+ K- pi0                                      PHSP;",
  "KKeta"             => "1.000 K+ K- eta                                      PHSP;",
  "ppbarpi0"          => "1.000 p+ anti-p- pi0                                 PHSP;",
  "pipieta"           => "1.000 pi+ pi- eta                                    PHSP;",
  "pipipi0pi0"        => "1.000 pi+ pi- pi0 pi0                                PHSP;",
  "2pipieta"          => "1.000 pi+ pi- pi+ pi- eta                            PHSP;",
  "2pipipi0"          => "1.000 pi+ pi- pi0 pi+ pi- pi0                        PHSP;"
}

# Build exclusive MC for each mode at each energy point
exclusive_mc = {}
[[data_4230, "4230"], [data_4260, "4260"], [data_4360, "4360"]].each do |ds, tag|
  eta_c_modes.each do |mode_name, decay_line|
    exclusive_mc["#{mode_name}_#{tag}"] = DatasetManager.create_exclusive_mc do |config|
      config.sample_name     = "pi0pi0hc_gam_etac_#{mode_name}_#{tag}"
      config.related_dataset = ds
      config.events          = 200000
      config.decay_card      = decay_card_for_mode(decay_line)
      config.cross_section   = :default
    end
  end
end

### Event selection (BOSS) ###
# Mode-dependent track/photon multiplicity and kinematic-fit participants.
# Track counts follow eta_c decay + gamma + 4 photons (2*pi0 -> gamma gamma).
# Photon counts: 5 gammas (1 from hc E1 + 2 for pi0 + 2 for pi0). Extra gammas for pi0/eta in eta_c.
# chi2_4C cut: <30 for eta_c modes with only charged tracks or K_S0; <25 otherwise.

# ---- Helper to build a Selection chain for a given eta_c mode ----
def build_selection_for(track_cfg, photon_min, kf_participants, chi2_cut_val)
  sel = Selection.new
  sel.select_track {
        cos_theta 0.93
        Vz        10.0
        Vr        1.0
        nChrp track_cfg[:nChrp]
        nChrn track_cfg[:nChrn]
        nNet  "==0"
      }
     .select_photon {
        nGam            ">=#{photon_min}"
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        tdc_emc_start   0
        tdc_emc_end     14
      }
     .pid(method: :probability) {
        prob_cut 0.001
        identify :pion,   against: [:kaon, :proton]
        identify :kaon,   against: [:pion, :proton]
        identify :proton, against: [:pion, :kaon]
      }
     .kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 200
        npi0 ">=2"
      }
     .kinematic_fit(kf_participants) {
        nominal
        constrain_four_momentum
        chi2_cut chi2_cut_val
        # gamma+pi0+pi0 recoil (eta_c region) and pi0 pi0 recoil (hc region)
        # are stored automatically for downstream ROOT analysis.
      }
  sel
end

# Mode-specific configurations
# key: [nChrp, nChrn, extra_photon_min_for_etac_modes, kf_participants_after_pi0pi0_hc, chi2_cut]
# The pi0 pi0 hc final state contributes 1 gamma (hc E1) + eta_c products.
# For pi0 pi0: reconstruct as two :pi0 via Kalman.
mode_configs = {
  "ppbar"       => { nChrp:"==1", nChrn:"==1", nGam:">=5",
                     kf:[:pi0, :pi0, :gamma, :prp, :prm],                                   chi2:30 },
  "pipiKK"      => { nChrp:"==2", nChrn:"==2", nGam:">=5",
                     kf:[:pi0, :pi0, :gamma, :pip, :pim, :kp, :km],                          chi2:30 },
  "pipippbar"   => { nChrp:"==2", nChrn:"==2", nGam:">=5",
                     kf:[:pi0, :pi0, :gamma, :pip, :pim, :prp, :prm],                        chi2:30 },
  "2KK"         => { nChrp:"==2", nChrn:"==2", nGam:">=5",
                     kf:[:pi0, :pi0, :gamma, :kp, :km, :kp, :km],                            chi2:30 },
  "2pipi"       => { nChrp:"==2", nChrn:"==2", nGam:">=5",
                     kf:[:pi0, :pi0, :gamma, :pip, :pim, :pip, :pim],                        chi2:30 },
  "3pipi"       => { nChrp:"==3", nChrn:"==3", nGam:">=5",
                     kf:[:pi0, :pi0, :gamma, :pip, :pim, :pip, :pim, :pip, :pim],            chi2:30 },
  "2pipiKK"     => { nChrp:"==3", nChrn:"==3", nGam:">=5",
                     kf:[:pi0, :pi0, :gamma, :pip, :pim, :pip, :pim, :kp, :km],              chi2:30 },
  "KsKpi"       => { nChrp:"==2", nChrn:"==2", nGam:">=5",
                     kf:[:pi0, :pi0, :gamma, :K_S0, :kp, :pim],                              chi2:30 },
  "KsKpipipi"   => { nChrp:"==3", nChrn:"==3", nGam:">=5",
                     kf:[:pi0, :pi0, :gamma, :K_S0, :kp, :pim, :pip, :pim],                  chi2:30 },
  "KKpi0"       => { nChrp:"==1", nChrn:"==1", nGam:">=7",
                     kf:[:pi0, :pi0, :pi0, :gamma, :kp, :km],                                chi2:25 },
  "KKeta"       => { nChrp:"==1", nChrn:"==1", nGam:">=7",
                     kf:[:pi0, :pi0, :eta, :gamma, :kp, :km],                                chi2:25 },
  "ppbarpi0"    => { nChrp:"==1", nChrn:"==1", nGam:">=7",
                     kf:[:pi0, :pi0, :pi0, :gamma, :prp, :prm],                              chi2:25 },
  "pipieta"     => { nChrp:"==1", nChrn:"==1", nGam:">=7",
                     kf:[:pi0, :pi0, :eta, :gamma, :pip, :pim],                              chi2:25 },
  "pipipi0pi0"  => { nChrp:"==1", nChrn:"==1", nGam:">=9",
                     kf:[:pi0, :pi0, :pi0, :pi0, :gamma, :pip, :pim],                        chi2:25 },
  "2pipieta"    => { nChrp:"==2", nChrn:"==2", nGam:">=7",
                     kf:[:pi0, :pi0, :eta, :gamma, :pip, :pim, :pip, :pim],                  chi2:25 },
  "2pipipi0"    => { nChrp:"==2", nChrn:"==2", nGam:">=9",
                     kf:[:pi0, :pi0, :pi0, :pi0, :gamma, :pip, :pim, :pip, :pim],            chi2:25 }
}

algorithms = []
mode_configs.each do |mode_name, cfg|
  alg_name = "PPHc_#{mode_name}"
  alg = Algorithm.new(alg_name)
  alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
     .set_constant({"ECMS" => [:double, 4.230]})  # per-energy ECMS set at execute_on time

  # For K_S0-containing modes, secondary vertex fit of pi+pi- pair
  sel = Selection.new
  sel.select_track {
        cos_theta 0.93
        Vz        10.0
        Vr        1.0
        nChrp cfg[:nChrp]
        nChrn cfg[:nChrn]
        nNet  "==0"
      }
     .select_photon {
        nGam            cfg[:nGam]
        energyThreshold_b 0.025
        energyThreshold_e 0.050
        tdc_emc_start   0
        tdc_emc_end     14
      }
     .pid(method: :probability) {
        prob_cut 0.001
        identify :pion,   against: [:kaon, :proton]
        identify :kaon,   against: [:pion, :proton]
        identify :proton, against: [:pion, :kaon]
      }

  # Reconstruct pi0 candidates (at least 2 pi0 for pi0 pi0 hc; some modes need more)
  n_pi0_expected = cfg[:kf].count(:pi0)
  sel.kalman_kinematic_fit([:gamma, :gamma]) {
       invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
       chi2_cut 200
       npi0 ">=#{n_pi0_expected}"
     }

  # Reconstruct eta candidates for modes with eta in eta_c decay
  if cfg[:kf].include?(:eta)
    sel.kalman_kinematic_fit([:gamma, :gamma]) {
         invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
         chi2_cut 200
         neta ">=1"
       }
  end

  # Reconstruct K_S0 for K_S0-containing modes via secondary vertex fit
  if cfg[:kf].include?(:K_S0)
    sel.assign({:chrgp => :pip, :chrgn => :pim})
       .secondary_vertex_fit([:pip, :pim]) {
          build_virtual_particle(:K_S0).by_minimizing_mass_difference
          remove_used_particle_from_candidate_list
        }
  end

  # Main 4C kinematic fit: total 4-momentum constrained to CMS; chi2_4C cut mode-dependent
  sel.kinematic_fit(cfg[:kf]) {
       nominal
       constrain_four_momentum
       chi2_cut cfg[:chi2]
     }

  alg.note(:pid_chi2_sum,
    "The best photon/pi0/eta combination is chosen by minimising " \
    "chi^2_total = chi^2_4C + sum_i chi^2_PID(i) + chi^2_1C(pi0/eta). " \
    "This composite chi^2 is not directly expressed in the DSL; the loose " \
    "chi2_4C cut here matches the paper's requirement (<30 for charged/KS0-only " \
    "modes, <25 for pi0/eta-containing modes).")
     .note(:etac_mass_window,
    "eta_c candidate mass window: |M(reconstructed) - m(eta_c)| < 35 MeV/c^2, " \
    "applied via M_recoil(gamma pi0 pi0) in [2.945, 3.015] GeV/c^2 (signal region). " \
    "Sideband regions [2.865, 2.900] and [3.050, 3.085] GeV/c^2 for background study.")
     .note(:hc_signal_window,
    "hc signal window: M_recoil(pi0 pi0) in [3.51, 3.55] GeV/c^2; " \
    "sidebands [3.45, 3.49] and [3.57, 3.61] GeV/c^2 for background study.")
     .note(:pi0_recoil_max,
    "Two pi0 recoil-mass combinations per event; retain the one with the LARGER " \
    "pi0 recoil-mass value (M_recoil(pi0)|max) for the Zc(4020)^0 fit.")
     .note(:preselection_windows,
    "Preselection: at least one gamma pi0 pi0 combination with " \
    "M_recoil(pi0 pi0) in [3.3, 3.7] GeV/c^2 (hc region) and " \
    "M_recoil(gamma pi0 pi0) in [2.8, 3.2] GeV/c^2 (eta_c region).")

  alg.with_decay_card(decay_card_for_mode(eta_c_modes[mode_name])).apply(sel)

  # Execute across all three CME points and their exclusive MC samples
  alg.execute_on([
    data_4230, incMC_4230, exclusive_mc["#{mode_name}_4230"],
    data_4260, incMC_4260, exclusive_mc["#{mode_name}_4260"],
    data_4360, incMC_4360, exclusive_mc["#{mode_name}_4360"]
  ])
  algorithms << alg
end
