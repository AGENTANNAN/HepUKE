# BESIII DSL: Observation of η_c → ω ω in J/ψ → γ ω ω
# ArXiv: 1905.10318v1
# J/ψ radiative decay, ω → π+π-π0, π0 → γγ

### Dataset preparation ###

# J/psi data at 3.097 GeV
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Decay card for J/ψ → γ η_c → γ ω ω → γ (π+π-π0)(π+π-π0)
decay_card_etac_ww = <<~DECAYCARD
  Decay J/psi
  1.000 gamma eta_c PHSP;
  Enddecay

  Decay eta_c
  1.000 omega omega PHSP;
  Enddecay

  Alias another_omega omega

  Decay omega
  1.000 pi+ pi- pi0 OMEGA_DALITZ;
  Enddecay

  Decay another_omega
  1.000 pi+ pi- pi0 OMEGA_DALITZ;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# Exclusive MC for signal
sig_etac_ww = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_etac_ww_exclusive_mc"
  config.related_dataset = jpsi_data
  config.events          = 2_000_000
  config.decay_card      = decay_card_etac_ww
  config.cross_section   = :default
end

### Event selection (BOSS) ###

alg_etac_ww = Algorithm.new("EtacToOmegaOmega")
alg_etac_ww.set_header(["EtacToOmegaOmegaAlg/EtacToOmegaOmega.h"])
            .set_constant({ "ECMS" => [:double, 3.097] })

sel_etac_ww = Selection.new
  # Select charged tracks: exactly 4 charged (2π+ 2π-), net zero charge
  .select_track {
    nChrp "==2"
    nChrn "==2"
    nNet "==0"
    cos_theta 0.93
    Vz 100.0
    Vr 10.0
  }
  # Select photons: 5 photons (1 radiative + 4 from two π0 → γγ)
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 10.0
    nGam ">=5"
  }
  # PID: identify pions
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip "==2"
    npim "==2"
  }
  # Remove identified pions from charged track lists
  .remove([:pip <= :chrgp, :pim <= :chrgn])
  # Reconstruct π0 from photon pairs with 1C Kalman mass constraint
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=2"
  }
  # 6C kinematic fit: J/ψ → γ ω ω → γ π+π-π0 π+π-π0 (π0 mass constrained already)
  .kinematic_fit([:gamma, :pip, :pim, :pip, :pim, :pi0, :pi0]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_etac_ww
  .note(:q_factor_background_subtraction, "Q-factor method used for background subtraction: each event assigned weight Q between 0 (pure background) and 1 (pure signal); ~12.5% of selected events from background, weighted out")
  .note(:omega_mass_window, "omega candidates required within 26 MeV/c^2 of nominal omega mass; signal region m(ωω) >= 2.65 GeV/c^2")
  .note(:pwa_fit, "Partial Wave Analysis (PWA) with unbinned maximum likelihood fit over complete phase space using helicity formalism; best hypothesis H0 = {η_c, 0-, 1+, 2+}; η_c yield = 1705 ± 58 events")
  .note(:etac_line_shape_modification, "η_c line shape modified relativistic Breit-Wigner with Eγ^(3/2) factor from M1 transition plus empirical damping factor exp(-Eγ^2/(16β^2)) with β = 0.065 GeV")
  .note(:track_helix_correction, "track helix parameters in MC smeared with Gaussian to match data resolution for kinematic fit efficiency systematic uncertainty")
  .note(:kinematic_fit_chi2_cut, "χ2_6C < 60 applied as selection criterion on the 6C kinematic fit")
  .with_decay_card(decay_card_etac_ww)
  .apply(sel_etac_ww)

alg_etac_ww.execute_on([jpsi_data, jpsi_incMC, sig_etac_ww])