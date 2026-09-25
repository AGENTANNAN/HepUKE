# Core DSL classes and dependencies will be loaded automatically at execution
# ============================================================================
# Search for the neutral Zc(4020)^0 in e+e- -> pi0 pi0 h_c, h_c -> gamma eta_c
# with eta_c decaying into 16 hadronic final states.
# Real data + inclusive MC at 4.23, 4.26 and 4.36 GeV.
# 200k-event exclusive MC for every eta_c mode at every energy point.
# ============================================================================

### Dataset preparation ###
data_4230  = DatasetManager.real_data.find("703_4230")    # e+e- data at sqrt(s) = 4.23 GeV
data_4260  = DatasetManager.real_data.find("703_4260")    # e+e- data at sqrt(s) = 4.26 GeV
data_4360  = DatasetManager.real_data.find("703_4360")    # e+e- data at sqrt(s) = 4.36 GeV
incMC_4230 = DatasetManager.inclusive_mc.find("703_4230") # inclusive MC at 4.23 GeV
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260") # inclusive MC at 4.26 GeV
incMC_4360 = DatasetManager.inclusive_mc.find("703_4360") # inclusive MC at 4.36 GeV

data_points = [data_4230, data_4260, data_4360]           # the three energy points
bkg_samples = [incMC_4230, incMC_4260, incMC_4360]        # the three inclusive MC samples

### Decay-card template (common signal chain, EvtGen syntax) ###
# Top mother is psi(4260) (BESIII / KKMC convention). Only the eta_c line varies per mode.
signal_decay_template = <<~DECAYCARD
  Decay psi(4260)
  1.0000 pi0 pi0 h_c PHSP;
  Enddecay

  Decay h_c
  1.000 gamma eta_c PHSP;
  Enddecay

  Decay eta_c
  __ETAC_MODE__
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  Decay eta
  1.000 gamma gamma PHSP;
  Enddecay

  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

### Per-mode specification table (16 eta_c hadronic modes) ###
#  nchr   : charged multiplicity per charge (1+1, 2+2 or 3+3)
#  ngam   : minimum photon multiplicity (5 = base, 7 = +one pi0/eta, 9 = +two pi0)
#  npi0   : minimum number of pi0 candidates (2 signal pi0 + mode pi0)
#  chi2   : 4C-fit chi2 cut (30 for charged-only / K_S0 modes, 25 for pi0/eta modes)
modes = [
  { key: "ppbar", alg: "Zc4020EtacPPbar",
    decay: "1.0000 anti-p- p+ PHSP;",
    nchr: "==1", ngam: ">=5",
    identify: [{ particles: [:proton], against: [:kaon, :pion] }],
    counts: { nprp: "==1", nprm: "==1" },
    npi0: ">=2", has_eta: false, has_ks: false,
    fit: [:pi0, :pi0, :gamma, :prp, :prm], chi2: 30 },

  { key: "pip_pim_Kp_Km", alg: "Zc4020EtacPiPiKK",
    decay: "1.0000 pi+ pi- K+ K- PHSP;",
    nchr: "==2", ngam: ">=5",
    identify: [{ particles: [:pion], against: [:kaon, :proton] },
               { particles: [:kaon], against: [:pion, :proton] }],
    counts: { npip: "==1", npim: "==1", nkp: "==1", nkm: "==1" },
    npi0: ">=2", has_eta: false, has_ks: false,
    fit: [:pi0, :pi0, :gamma, :pip, :pim, :kp, :km], chi2: 30 },

  { key: "pip_pim_ppbar", alg: "Zc4020EtacPiPiPPbar",
    decay: "1.0000 pi+ pi- anti-p- p+ PHSP;",
    nchr: "==2", ngam: ">=5",
    identify: [{ particles: [:pion], against: [:kaon, :proton] },
               { particles: [:proton], against: [:kaon, :pion] }],
    counts: { npip: "==1", npim: "==1", nprp: "==1", nprm: "==1" },
    npi0: ">=2", has_eta: false, has_ks: false,
    fit: [:pi0, :pi0, :gamma, :pip, :pim, :prp, :prm], chi2: 30 },

  { key: "Kp_Km_Kp_Km", alg: "Zc4020EtacKKKK",
    decay: "1.0000 K+ K- K+ K- PHSP;",
    nchr: "==2", ngam: ">=5",
    identify: [{ particles: [:kaon], against: [:pion, :proton] }],
    counts: { nkp: "==2", nkm: "==2" },
    npi0: ">=2", has_eta: false, has_ks: false,
    fit: [:pi0, :pi0, :gamma, :kp, :kp, :km, :km], chi2: 30 },

  { key: "2pip_2pim", alg: "Zc4020Etac2Pi2Pi",
    decay: "1.0000 pi+ pi- pi+ pi- PHSP;",
    nchr: "==2", ngam: ">=5",
    identify: [{ particles: [:pion], against: [:kaon, :proton] }],
    counts: { npip: "==2", npim: "==2" },
    npi0: ">=2", has_eta: false, has_ks: false,
    fit: [:pi0, :pi0, :gamma, :pip, :pip, :pim, :pim], chi2: 30 },

  { key: "3pip_3pim", alg: "Zc4020Etac3Pi3Pi",
    decay: "1.0000 pi+ pi- pi+ pi- pi+ pi- PHSP;",
    nchr: "==3", ngam: ">=5",
    identify: [{ particles: [:pion], against: [:kaon, :proton] }],
    counts: { npip: "==3", npim: "==3" },
    npi0: ">=2", has_eta: false, has_ks: false,
    fit: [:pi0, :pi0, :gamma, :pip, :pip, :pip, :pim, :pim, :pim], chi2: 30 },

  { key: "2pip_2pim_Kp_Km", alg: "Zc4020Etac2Pi2PiKK",
    decay: "1.0000 pi+ pi- pi+ pi- K+ K- PHSP;",
    nchr: "==3", ngam: ">=5",
    identify: [{ particles: [:pion], against: [:kaon, :proton] },
               { particles: [:kaon], against: [:pion, :proton] }],
    counts: { npip: "==2", npim: "==2", nkp: "==1", nkm: "==1" },
    npi0: ">=2", has_eta: false, has_ks: false,
    fit: [:pi0, :pi0, :gamma, :pip, :pip, :pim, :pim, :kp, :km], chi2: 30 },

  { key: "KS_K_pi", alg: "Zc4020EtacKSKPi",
    decay: "1.0000 K_S0 K+ pi- PHSP;",
    nchr: "==2", ngam: ">=5",
    identify: [{ particles: [:pion], against: [:kaon, :proton] },
               { particles: [:kaon], against: [:pion, :proton] }],
    counts: { npip: ">=1", npim: ">=1", nkp: ">=1" },
    npi0: ">=2", has_eta: false, has_ks: true,
    fit: [:pi0, :pi0, :gamma, :K_S0, :kp, :pim], chi2: 30 },

  { key: "KS_K_3pi", alg: "Zc4020EtacKSK3Pi",
    decay: "1.0000 K_S0 K- pi+ pi+ pi- PHSP;",
    nchr: "==3", ngam: ">=5",
    identify: [{ particles: [:pion], against: [:kaon, :proton] },
               { particles: [:kaon], against: [:pion, :proton] }],
    counts: { npip: ">=1", npim: ">=1", nkm: ">=1" },
    npi0: ">=2", has_eta: false, has_ks: true,
    fit: [:pi0, :pi0, :gamma, :K_S0, :km, :pip, :pip, :pim], chi2: 30 },

  { key: "Kp_Km_pi0", alg: "Zc4020EtacKKPi0",
    decay: "1.0000 K+ K- pi0 PHSP;",
    nchr: "==1", ngam: ">=7",
    identify: [{ particles: [:kaon], against: [:pion, :proton] }],
    counts: { nkp: "==1", nkm: "==1" },
    npi0: ">=3", has_eta: false, has_ks: false,
    fit: [:pi0, :pi0, :gamma, :pi0, :kp, :km], chi2: 25 },

  { key: "Kp_Km_eta", alg: "Zc4020EtacKKEta",
    decay: "1.0000 K+ K- eta PHSP;",
    nchr: "==1", ngam: ">=7",
    identify: [{ particles: [:kaon], against: [:pion, :proton] }],
    counts: { nkp: "==1", nkm: "==1" },
    npi0: ">=2", has_eta: true, has_ks: false,
    fit: [:pi0, :pi0, :gamma, :eta, :kp, :km], chi2: 25 },

  { key: "ppbar_pi0", alg: "Zc4020EtacPPbarPi0",
    decay: "1.0000 anti-p- p+ pi0 PHSP;",
    nchr: "==1", ngam: ">=7",
    identify: [{ particles: [:proton], against: [:kaon, :pion] }],
    counts: { nprp: "==1", nprm: "==1" },
    npi0: ">=3", has_eta: false, has_ks: false,
    fit: [:pi0, :pi0, :gamma, :pi0, :prp, :prm], chi2: 25 },

  { key: "pip_pim_eta", alg: "Zc4020EtacPiPiEta",
    decay: "1.0000 pi+ pi- eta PHSP;",
    nchr: "==1", ngam: ">=7",
    identify: [{ particles: [:pion], against: [:kaon, :proton] }],
    counts: { npip: "==1", npim: "==1" },
    npi0: ">=2", has_eta: true, has_ks: false,
    fit: [:pi0, :pi0, :gamma, :eta, :pip, :pim], chi2: 25 },

  { key: "pip_pim_2pi0", alg: "Zc4020EtacPiPi2Pi0",
    decay: "1.0000 pi+ pi- pi0 pi0 PHSP;",
    nchr: "==1", ngam: ">=9",
    identify: [{ particles: [:pion], against: [:kaon, :proton] }],
    counts: { npip: "==1", npim: "==1" },
    npi0: ">=4", has_eta: false, has_ks: false,
    fit: [:pi0, :pi0, :gamma, :pi0, :pi0, :pip, :pim], chi2: 25 },

  { key: "2pip_2pim_eta", alg: "Zc4020Etac2Pi2PiEta",
    decay: "1.0000 pi+ pi- pi+ pi- eta PHSP;",
    nchr: "==2", ngam: ">=7",
    identify: [{ particles: [:pion], against: [:kaon, :proton] }],
    counts: { npip: "==2", npim: "==2" },
    npi0: ">=2", has_eta: true, has_ks: false,
    fit: [:pi0, :pi0, :gamma, :eta, :pip, :pip, :pim, :pim], chi2: 25 },

  { key: "2pip_2pim_pi0", alg: "Zc4020Etac2Pi2PiPi0",
    decay: "1.0000 pi+ pi- pi+ pi- pi0 PHSP;",
    nchr: "==2", ngam: ">=7",
    identify: [{ particles: [:pion], against: [:kaon, :proton] }],
    counts: { npip: "==2", npim: "==2" },
    npi0: ">=3", has_eta: false, has_ks: false,
    fit: [:pi0, :pi0, :gamma, :pi0, :pip, :pip, :pim, :pim], chi2: 25 }
]

### Decay cards + exclusive MC (200k events per mode per energy point) ###
decay_cards   = {}
exclusive_mcs = {}
modes.each do |m|
  decay_cards[m[:key]] = signal_decay_template.sub("__ETAC_MODE__", m[:decay])

  # One exclusive MC per energy point, sharing the same decay card / cross section / statistics
  exclusive_mcs[m[:key]] = DatasetManager.create_exclusive_mc_for(data_points) do |config|
    config.sample_name   = "exmc_zc4020_#{m[:key]}"
    config.events        = 200_000
    config.decay_card    = decay_cards[m[:key]]
    config.cross_section = :default
  end
end

### Event selection (BOSS) + execution, one Algorithm per eta_c decay mode ###
algorithms = {}
modes.each do |m|
  alg = Algorithm.new(m[:alg])
  alg.set_header(["#{m[:alg]}Alg/#{m[:alg]}.h"])
     .set_constant({ "ECMS" => [:double, 4.26] })   # central energy point (see note below)

  sel = Selection.new
  sel.select_track {                                # charged-track selection
        cos_theta 0.93                              # |cos(theta)| < 0.93
        Vz        10.0                              # |Vz| < 10 cm
        Vr        1.0                               # Vr < 1 cm
        nChrp     m[:nchr]                          # mode-dependent positive multiplicity
        nChrn     m[:nchr]                          # mode-dependent negative multiplicity
        nNet      "==0"                             # net charge zero
      }
     .select_photon {                               # photon selection
        tdc_emc_start     0                         # TDC window 0-14
        tdc_emc_end       14
        angle_to_track    10.0
        energyThreshold_b 0.025                     # E > 25 MeV (barrel)
        energyThreshold_e 0.050                     # E > 50 MeV (endcap)
        nGam              m[:ngam]                  # >= 5 / 7 / 9 photons
      }
     .pid(method: :probability) {                   # PID: probability method, 0.001 cut
        prob_cut 0.001
        m[:identify].each { |spec| identify(*spec[:particles], against: spec[:against]) }
        m[:counts].each   { |name, expr| public_send(name, expr) }  # counts matched to final state
      }

  # pi0 -> gamma gamma mass-constrained reconstruction (>= 2 signal pi0)
  sel.kalman_kinematic_fit([:gamma, :gamma]) {
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 200
        npi0 m[:npi0]
      }

  # eta -> gamma gamma mass-constrained reconstruction (only for eta modes)
  if m[:has_eta]
    sel.kalman_kinematic_fit([:gamma, :gamma]) {
          invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
          chi2_cut 200
          neta ">=1"
        }
  end

  # K_S0 -> pi+ pi- through a secondary vertex fit (only for K_S0 modes)
  if m[:has_ks]
    sel.secondary_vertex_fit([:pip, :pim]) {
          build_virtual_particle(:K_S0).by_minimizing_mass_difference
          remove_used_particle_from_candidate_list
        }
  end

  # 4C kinematic fit to the nominal four-momentum; best photon/pi0/eta combination
  # is chosen automatically by minimizing the composite chi2.
  sel.kinematic_fit(m[:fit]) {
        nominal
        constrain_four_momentum
        chi2_cut m[:chi2]                           # 30 (charged/K_S0) or 25 (pi0/eta) modes
      }

  alg.note(:recoil_mass_preselection,
           "preselection window on the gamma pi0 pi0 system (h_c region [3.3, 3.7] GeV/c^2, " \
           "eta_c region [2.8, 3.2] GeV/c^2); recoil masses and the Zc-side mass-window / " \
           "sideband / larger-of-two-pi0-recoil-mass selection are applied after the 4C fit " \
           "in ROOT and are not expressible in the BOSS DSL selection chain")
     .note(:multi_energy_ecms,
           "this algorithm runs on three energy points (4.23/4.26/4.36 GeV); the 4C fit " \
           "requires the per-point CMS energy, here fixed through the ECMS constant at the " \
           "central 4.26 GeV value and to be overridden per energy point")
     .with_decay_card(decay_cards[m[:key]])
     .apply(sel)

  # Execute on real data + inclusive MC (all three energies) + the mode's exclusive MC
  alg.execute_on(data_points + bkg_samples + exclusive_mcs[m[:key]])
  algorithms[m[:key]] = alg
end