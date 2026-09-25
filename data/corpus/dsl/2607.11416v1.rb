# psi(3686) -> gamma eta_c -> gamma p pbar eta, with eta -> gamma gamma
# (2712.4 +/- 14.3) x 10^6 psi(3686) events

psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

decay_card_sig = <<~DECAYCARD
  Decay psi(2S)
  1.0000 gamma eta_c                                PHSP;
  Enddecay

  Decay eta_c
  1.0000 p+ anti-p- eta                             PHSP;
  Enddecay

  Decay eta
  1.0000 gamma gamma                                PHSP;
  Enddecay

  End
DECAYCARD

exMC_sig = DatasetManager.create_exclusive_mc do |c|
  c.sample_name     = "psip_gamma_etac_ppbar_eta"
  c.related_dataset = psip_data
  c.events          = 500000
  c.decay_card      = decay_card_sig
  c.cross_section   = :default
end

alg_name  = "PsipGammaPPbarEta"
algorithm = Algorithm.new(alg_name)
algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 3.686]})

event_selection = Selection.new
event_selection
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     "==1"
    nChrn     "==1"
    nNet      "==0"
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam              ">=3"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp "==1"
    nprm "==1"
  }
  .select_isolated_photon {
    angle_to_prm_track 20.0
    nGam ">=3"
  }
  # 4C kinematic fit for psi(3686) -> p pbar gamma gamma gamma
  .kinematic_fit([:prp, :prm, :gamma, :gamma, :gamma]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
  # Alternative 2-gamma and 4-gamma hypotheses to store chi2_4C for background suppression
  .kinematic_fit([:prp, :prm, :gamma, :gamma]) {
    constrain_four_momentum
  }
  .kinematic_fit([:prp, :prm, :gamma, :gamma, :gamma, :gamma]) {
    constrain_four_momentum
  }

algorithm
  .note(:chi2_4c_tight_in_root, "chi2_4C(3gamma p pbar) < 65 requirement applied in ROOT after loose BOSS-level chi2_cut 200")
  .note(:chi2_hypothesis_veto,  "require chi2_4C(3gamma p pbar) < chi2_4C(2gamma p pbar) AND chi2_4C(3gamma p pbar) < chi2_4C(4gamma p pbar) to suppress events better described by 2gamma or 4gamma hypotheses")
  .note(:eta_photon_pair_choice, "photons ordered by energy g1>g2>g3; select (g2,g3) for eta if M(g2 g3) in (0.4, 0.7) GeV/c^2, else try (g1,g3); else reject event")
  .note(:pi0_veto,              "reject events where any gamma gamma invariant mass falls in the pi0 window ~(0.12, 0.15) GeV/c^2 (FOM-optimized)")
  .note(:jpsi_veto,             "reject events where recoil mass of (g2 g3) falls in J/psi window (3.072, 3.125) GeV/c^2 to suppress psi(3686) -> eta J/psi")
  .note(:eta_mass_signal_region,"eta signal window (0.520, 0.570) GeV/c^2; sidebands (0.475, 0.500) U (0.600, 0.625) GeV/c^2 for background estimation")

algorithm.with_decay_card(decay_card_sig).apply(event_selection)
algorithm.execute_on([psip_data, psip_incMC, exMC_sig])
