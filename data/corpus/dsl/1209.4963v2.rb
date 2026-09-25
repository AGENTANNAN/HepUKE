# =============================================================================
# BESIII paper 1209.4963v2 — Measurement of the mass and width of h_c and
#   eta_c, and of eta_c exclusive branching fractions, via
#   psi(3686) -> pi0 h_c, h_c -> gamma eta_c, eta_c -> X_i
# arXiv:1209.4963v2
# Data: (106 +- 4) x 10^6 psi(3686) events at sqrt(s) = 3.686 GeV,
#       plus 42 pb^-1 continuum data at sqrt(s) = 3.65 GeV.
#
# The eta_c is reconstructed in 16 exclusive hadronic decay modes X_i:
#   p pbar, 2(pi+pi-), 2(K+K-), K+K-pi+pi-, p pbar pi+pi-, 3(pi+pi-),
#   K+K-2(pi+pi-), K+K-pi0, p pbar pi0, K_S0 K+-pi-+, K_S0 K+-pi-+pi+pi-,
#   pi+pi-eta, K+K-eta, 2(pi+pi-)eta, pi+pi-pi0pi0, 2(pi+pi-)pi0pi0.
# Each mode has its own Algorithm and Selection chain (different track/photon
# multiplicities and different 4C-fit participant lists).
# =============================================================================

### Dataset preparation ###
psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")
cont_data  = DatasetManager.real_data.find("709_3650")     # continuum at 3.65 GeV
cont_incMC = DatasetManager.inclusive_mc.find("709_3650")

# ---- Decay cards -------------------------------------------------------------
# Common template: psi(2S) -> pi0 h_c, h_c -> gamma eta_c, eta_c -> X_i.
# The eta_c decay line varies by mode; the soft pi0 always decays to gamma gamma.
def decay_card_for_mode(eta_c_decay_line)
  <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi0 h_c                  PHSP;
    Enddecay

    Decay h_c
    1.0000 gamma eta_c              PHSP;
    Enddecay

    Decay eta_c
    #{eta_c_decay_line}
    Enddecay

    Decay pi0
    1.0000 gamma gamma              PHSP;
    Enddecay

    End
  DECAYCARD
end

# The 16 exclusive eta_c hadronic final states
eta_c_modes = {
  "ppbar"            => "1.0000 p+ anti-p-                             PHSP;",
  "2pipi"            => "1.0000 pi+ pi- pi+ pi-                         PHSP;",
  "2KK"              => "1.0000 K+ K- K+ K-                             PHSP;",
  "KKpipi"           => "1.0000 K+ K- pi+ pi-                          PHSP;",
  "ppbarpipi"        => "1.0000 p+ anti-p- pi+ pi-                      PHSP;",
  "3pipi"            => "1.0000 pi+ pi- pi+ pi- pi+ pi-                 PHSP;",
  "KK2pipi"          => "1.0000 K+ K- pi+ pi- pi+ pi-                  PHSP;",
  "KKpi0"            => "1.0000 K+ K- pi0                              PHSP;",
  "ppbarpi0"         => "1.0000 p+ anti-p- pi0                          PHSP;",
  "KsKpi"            => "1.0000 K_S0 K+ pi-                            PHSP;",
  "KsKpipipi"        => "1.0000 K_S0 K+ pi- pi+ pi-                    PHSP;",
  "pipieta"          => "1.0000 pi+ pi- eta                             PHSP;",
  "KKeta"            => "1.0000 K+ K- eta                              PHSP;",
  "2pipieta"         => "1.0000 pi+ pi- pi+ pi- eta                    PHSP;",
  "pipipi0pi0"       => "1.0000 pi+ pi- pi0 pi0                        PHSP;",
  "2pipipi0pi0"      => "1.0000 pi+ pi- pi+ pi- pi0 pi0                PHSP;"
}

# ---- Exclusive MC samples (one per eta_c decay mode) -------------------------
exclusive_mc = {}
eta_c_modes.each do |mode_name, decay_line|
  exclusive_mc[mode_name] = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "psip_pi0_hc_gamma_etac_#{mode_name}"
    config.related_dataset = psip_data
    config.events          = 200000
    config.decay_card      = decay_card_for_mode(decay_line)
    config.cross_section   = :default
  end
end

### Event selection (BOSS) ###
# The final state psi(3686) -> pi0 h_c, h_c -> gamma eta_c -> gamma X_i contains
#   1 soft pi0 (2 photons) + 1 E1 photon + the eta_c decay products X_i.
# Track and photon multiplicities therefore follow the X_i final state; the
# per-channel chi2_4C and PID requirements are the optimised ones of Table I.
# The paper quotes M(pi0 particle momentum) ~ 84 MeV/c and the E1 photon
# energy ~ 503 MeV in the h_c rest frame.

# Helper: build the Selection chain for one eta_c mode.
# cfg keys: :nChrp, :nChrn, :nGam, :npi0, :neta (for kalman reconstruction),
#           :ks0 (reconstruct K_S0), :kf (4C-fit participants),
#           :pid_tag (Table I PID requirement), :chi2 (paper's chi2_4C value)
def build_selection_for(cfg)
  sel = Selection.new
  sel.select_track {
        cos_theta 0.93    # |cos(theta)| < 0.93
        Vz        10.0    # within +-10 cm of the IP along the beam direction
        Vr         1.0    # within +-1 cm in the plane perpendicular to the beam
        nChrp cfg[:nChrp]
        nChrn cfg[:nChrn]
        nNet  "==0"
      }
     .select_photon {
        nGam              cfg[:nGam]
        energyThreshold_b 0.025   # barrel (|cos theta| < 0.8): E > 25 MeV
        energyThreshold_e 0.050   # endcap (0.86 < |cos theta| < 0.92): E > 50 MeV
        angle_to_track    10.0
        tdc_emc_start     0
        tdc_emc_end       14
      }

  # Particle identification (channel dependent; modes with no PID requirement
  # in Table I of the paper get no identify() block)
  case cfg[:pid_tag]
  when :p
    sel.pid(method: :probability) {
      prob_cut 0.001
      identify :proton, against: [:pion, :kaon]
      nprp ">=1"
    }
  when :p2
    sel.pid(method: :probability) {
      prob_cut 0.001
      identify :proton, against: [:pion, :kaon]
      identify :pion,   against: [:kaon, :proton]
      nprp ">=2"
    }
  when :pi3
    sel.pid(method: :probability) {
      prob_cut 0.001
      identify :pion, against: [:kaon, :proton]
      npip ">=3"
    }
  when :pi4
    sel.pid(method: :probability) {
      prob_cut 0.001
      identify :pion, against: [:kaon, :proton]
      npip ">=4"
    }
  when :k3
    sel.pid(method: :probability) {
      prob_cut 0.001
      identify :kaon, against: [:pion, :proton]
      nkp ">=3"
    }
  when :k1
    sel.pid(method: :probability) {
      prob_cut 0.001
      identify :kaon, against: [:pion, :proton]
      nkp ">=1"
    }
  when :k2pi0
    sel.pid(method: :probability) {
      prob_cut 0.001
      identify :kaon, against: [:pion, :proton]
      nkp ">=2"
    }
  when :k2pi2
    sel.pid(method: :probability) {
      prob_cut 0.001
      identify :kaon, against: [:pion, :proton]
      identify :pion, against: [:kaon, :proton]
      nkp  ">=2"
      npip ">=2"
    }
  end

  # Reconstruct the soft pi0 from psi(3686) -> pi0 h_c (and any pi0 in X_i):
  # |M(gamma gamma) - m_pi0| < 15 MeV/c^2 with a 1C mass-constrained fit.
  if cfg[:npi0].to_i > 0
    sel.kalman_kinematic_fit([:gamma, :gamma]) {
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 25
      npi0 ">=#{cfg[:npi0]}"
    }
  end

  # Reconstruct eta from gamma gamma for the modes containing eta
  if cfg[:neta].to_i > 0
    sel.kalman_kinematic_fit([:gamma, :gamma]) {
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
      chi2_cut 25
      neta ">=#{cfg[:neta]}"
    }
  end

  # Reconstruct K_S0 -> pi+ pi- through a secondary vertex fit for K_S0 modes
  if cfg[:ks0]
    sel.assign({chrgp: :pip, chrgn: :pim})
       .secondary_vertex_fit([:pip, :pim]) {
         build_virtual_particle(:K_S0).by_minimizing_mass_difference
         remove_used_particle_from_candidate_list
       }
  end

  # Main 4C kinematic fit under the psi(3686) -> pi0 gamma X_i hypothesis
  sel.kinematic_fit(cfg[:kf]) {
    nominal
    constrain_four_momentum
    chi2_cut 200   # loose BOSS cut; the paper's optimised channel-dependent
                   # chi2_4C value (Table I) is applied in ROOT
  }
  sel
end

# Mode configurations.
#   npi0 counts the soft pi0 from psi(3686) -> pi0 h_c plus any pi0 in X_i;
#   neta counts eta in X_i; nGam = 2 (soft pi0) + 1 (E1 photon) + photons in X_i.
#   pid_tag encodes the Table I PID requirement; chi2 the Table I chi2_4C value.
mode_configs = {
  "ppbar" => {
    nChrp: "==1", nChrn: "==1", nGam: ">=3", npi0: 1, neta: 0, ks0: false,
    pid_tag: :p,
    kf: [:pi0, :gamma, :prp, :prm], chi2: 30 },
  "2pipi" => {
    nChrp: "==2", nChrn: "==2", nGam: ">=3", npi0: 1, neta: 0, ks0: false,
    pid_tag: :pi3,
    kf: [:pi0, :gamma, :pip, :pim, :pip, :pim], chi2: 60 },
  "2KK" => {
    nChrp: "==2", nChrn: "==2", nGam: ">=3", npi0: 1, neta: 0, ks0: false,
    pid_tag: :k3,
    kf: [:pi0, :gamma, :kp, :km, :kp, :km], chi2: 60 },
  "KKpipi" => {
    nChrp: "==2", nChrn: "==2", nGam: ">=3", npi0: 1, neta: 0, ks0: false,
    pid_tag: :k2pi0,
    kf: [:pi0, :gamma, :kp, :km, :pip, :pim], chi2: 40 },
  "ppbarpipi" => {
    nChrp: "==2", nChrn: "==2", nGam: ">=3", npi0: 1, neta: 0, ks0: false,
    pid_tag: :p2,
    kf: [:pi0, :gamma, :prp, :prm, :pip, :pim], chi2: 30 },
  "3pipi" => {
    nChrp: "==3", nChrn: "==3", nGam: ">=3", npi0: 1, neta: 0, ks0: false,
    pid_tag: :pi4,
    kf: [:pi0, :gamma, :pip, :pim, :pip, :pim, :pip, :pim], chi2: 50 },
  "KK2pipi" => {
    nChrp: "==3", nChrn: "==3", nGam: ">=3", npi0: 1, neta: 0, ks0: false,
    pid_tag: :k2pi2,
    kf: [:pi0, :gamma, :kp, :km, :pip, :pim, :pip, :pim], chi2: 70 },
  "KKpi0" => {
    nChrp: "==1", nChrn: "==1", nGam: ">=5", npi0: 2, neta: 0, ks0: false,
    pid_tag: :k1,
    kf: [:pi0, :gamma, :pi0, :kp, :km], chi2: 50 },
  "ppbarpi0" => {
    nChrp: "==1", nChrn: "==1", nGam: ">=5", npi0: 2, neta: 0, ks0: false,
    pid_tag: :p,
    kf: [:pi0, :gamma, :pi0, :prp, :prm], chi2: 40 },
  "KsKpi" => {
    nChrp: "==2", nChrn: "==2", nGam: ">=3", npi0: 1, neta: 0, ks0: true,
    pid_tag: nil,
    kf: [:pi0, :gamma, :K_S0, :kp, :pim], chi2: 70 },
  "KsKpipipi" => {
    nChrp: "==3", nChrn: "==3", nGam: ">=3", npi0: 1, neta: 0, ks0: true,
    pid_tag: nil,
    kf: [:pi0, :gamma, :K_S0, :kp, :pim, :pip, :pim], chi2: 50 },
  "pipieta" => {
    nChrp: "==1", nChrn: "==1", nGam: ">=5", npi0: 1, neta: 1, ks0: false,
    pid_tag: nil,
    kf: [:pi0, :gamma, :eta, :pip, :pim], chi2: 50 },
  "KKeta" => {
    nChrp: "==1", nChrn: "==1", nGam: ">=5", npi0: 1, neta: 1, ks0: false,
    pid_tag: :k1,
    kf: [:pi0, :gamma, :eta, :kp, :km], chi2: 70 },
  "2pipieta" => {
    nChrp: "==2", nChrn: "==2", nGam: ">=5", npi0: 1, neta: 1, ks0: false,
    pid_tag: nil,
    kf: [:pi0, :gamma, :eta, :pip, :pim, :pip, :pim], chi2: 30 },
  "pipipi0pi0" => {
    nChrp: "==1", nChrn: "==1", nGam: ">=7", npi0: 3, neta: 0, ks0: false,
    pid_tag: nil,
    kf: [:pi0, :gamma, :pi0, :pi0, :pip, :pim], chi2: 40 },
  "2pipipi0pi0" => {
    nChrp: "==2", nChrn: "==2", nGam: ">=7", npi0: 3, neta: 0, ks0: false,
    pid_tag: nil,
    kf: [:pi0, :gamma, :pi0, :pi0, :pip, :pim, :pip, :pim], chi2: 60 }
}

algorithms = []
mode_configs.each do |mode_name, cfg|
  alg_name = "PsipPi0Hc_#{mode_name}"
  alg = Algorithm.new(alg_name)
  alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
     .set_constant({"ECMS" => [:double, 3.686]})

  sel = build_selection_for(cfg)

  alg.note(:pi0_eta_reconstruction,
           "pi0 (eta) candidates are reconstructed from photon pairs with " \
           "|M(gamma gamma) - m_pi0(eta)| < 15 MeV/c^2 and a 1C kinematic fit constraining " \
           "M(gamma gamma) to the known pi0 (eta) mass; the pi0 (eta) mass window and the " \
           "1C-fit chi2 enter the composite candidate-quality chi2 used to select the best " \
           "combination per event.")
     .note(:ks0_reconstruction,
           "K_S0 -> pi+ pi- candidates require |M(pi+ pi-) - m_K_S0| < 20 MeV/c^2; a " \
           "secondary-vertex fit imposing the kinematic constraint between the production and " \
           "decay vertices is used to reject random pi+ pi- combinations and the accepted K_S0 " \
           "candidates must have a decay length of at least twice the vertex resolution. Tracks " \
           "used in reconstructing K_S0 decays are exempt from the Vz and Vr requirements on " \
           "the primary track selection.")
     .note(:photon_selection,
           "photon candidates are EMC showers with E > 25 MeV in the barrel (|cos theta| < 0.8) " \
           "and E > 50 MeV in the end-cap (0.86 < |cos theta| < 0.92); showers in the " \
           "intermediate barrel/end-cap region are not well measured and are rejected, and an " \
           "EMC timing requirement suppresses electronic noise and energy deposits unrelated to " \
           "the event.")
     .note(:etac_candidate_selection,
           "the reconstructed eta_c mass is required to satisfy 2.900 < M(eta_c) < 3.050 GeV/c^2 " \
           "and the E1 transition photon energy 0.450 < E(gamma_E1) < 0.550 GeV. Applied in ROOT.")
     .note(:best_candidate_selection,
           "when several h_c candidates are found in an event, the one with the smallest value " \
           "of chi2 = chi2_4C + chi2_1C + chi2_pid + chi2_vertex is accepted (chi2_1C is the " \
           "1C-fit chi2 of the pi0/eta, chi2_pid the PID chi2 summed over the charged tracks of " \
           "the h_c candidate, and chi2_vertex the K_S0 vertex-fit chi2). If an event has no " \
           "pi0/eta (K_S0), the corresponding chi2_1C (chi2_vertex) is set to zero. The composite " \
           "chi2 is not directly expressible in the DSL and is applied in ROOT. The paper's " \
           "optimised channel-dependent chi2_4C cut (this mode: <#{cfg[:chi2]}) is likewise " \
           "applied in ROOT instead of the loose BOSS cut.")
     .note(:background_veto,
           "five background-suppression requirements are applied channel by channel on the basis " \
           "of the 100-million-event inclusive MC study: (1) reject events with the pi+ pi- " \
           "recoil mass M_X within +-12 MeV/c^2 of the J/psi mass; (2) reject events with the " \
           "pi0 pi0 recoil mass within +-15 MeV/c^2 of the J/psi mass (+-10 MeV/c^2 for the " \
           "pi+pi-pi0pi0 mode, where the lower pi0 momentum moves the recoil mass near " \
           "3.1 GeV/c^2); (3) reject candidates containing a pi0 one of whose daughter photons " \
           "has an energy within +-5 MeV of that expected for the psi(3686) -> gamma chi_c2 " \
           "radiative transition (128 MeV); (4) reject an E1 photon candidate that can be " \
           "combined with another photon in the event to form a pi0 within +-10 MeV/c^2; " \
           "(5) reject a pi0 candidate whose M(pi+ pi- pi0) combination (for any combination in " \
           "the event) is within +-15 MeV/c^2 of the known eta mass.")
     .note(:signal_extraction,
           "the h_c mass, width and per-channel yields are extracted from a simultaneous binned " \
           "maximum-likelihood fit of the pi0 recoil mass distributions of the 16 channels, " \
           "using a Breit-Wigner signal shape convolved with the channel-dependent MC resolution " \
           "and a background shape taken from the eta_c mass side bands (2300-2700 and " \
           "3070-3200 MeV/c^2). Only 1C pi0-mass-constrained fits are used to reconstruct the " \
           "pi0 recoil mass (not the 4C-fit four-momenta). The eta_c resonant parameters are " \
           "extracted from the 16 hadronic invariant-mass spectra accompanying the transition " \
           "photon, fitted with (E_gamma^3 x BW(m) x f_d(E_gamma)) convolved with a " \
           "double-Gaussian resolution, with the KEDR damping function and background from the " \
           "h_c side bands (3500-3515, 3535-3550 MeV/c^2). ROOT level.")
     .note(:systematics,
           "h_c mass/width: energy calibration (0.13 MeV/c^2, 0.07 MeV), signal shape, fitting " \
           "range 0.04 MeV/c^2, binning 0.02 MeV/c^2, background shape 0.16 MeV/c^2, background " \
           "veto, kinematic fit 0.01 MeV/c^2 and M(psi(3686)) 0.03 MeV/c^2 (total 0.14 MeV/c^2 " \
           "and 0.22 MeV). eta_c branching ratios: N(psi(3686)) 4.0%, tracking 2% per track, " \
           "kaon/pion PID 2%, photon detection 1% per photon, fit range, background shape, " \
           "signal shape 2.3%, kinematic-fit efficiency, background veto, cross-feed (2.5%, 1.4%, " \
           "1.3% for 2(pi+pi-), K+K-pi0, pi+pi-pi0pi0 from K_S0 K+-pi-+), eta_c decay models and " \
           "eta_c line shape; total per mode 10.6%-19.4%. K_S0 reconstruction 1%. ROOT level.")

  alg.with_decay_card(decay_card_for_mode(eta_c_modes[mode_name])).apply(sel)
  alg.execute_on([psip_data, psip_incMC, cont_data, cont_incMC,
                  exclusive_mc[mode_name]])
  algorithms << alg
end
