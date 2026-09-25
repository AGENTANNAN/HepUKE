# Paper: arXiv:2404.03217v1
# Evidence of h_c -> K_S0 K+ pi- + c.c. in psi(3686) -> pi0 h_c
# Data: (2.712+-0.014)*10^9 psi(3686) at sqrt(s) = 3.686 GeV
# h_c reconstructed from K_S0 K+ pi-; pi0 from gamma gamma; K_S0 from pi+ pi-

### Dataset preparation ###
psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

decay_card = <<~DECAYCARD
    Decay psi(3686)
    1.000 pi0 h_c PHSP;
    Enddecay
    Decay h_c
    1.000 K_S0 K+ pi- PHSP;
    Enddecay
    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay
    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

exMC = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_psip_pi0_hc_KsKpi"
  config.related_dataset = psip_data
  config.events          = 1_000_000
  config.decay_card      = decay_card
  config.cross_section   = :default
end

### Event selection ###
alg = Algorithm.new("PsipPi0HcKsKPi")
alg.set_header(["PsipPi0HcKsKPiAlg/PsipPi0HcKsKPi.h"])
    .set_constant({ "ECMS" => [:double, 3.686] })

event_selection = Selection.new
  .select_track do
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=2"   # K+, pi+ (from K_S0)
    nChrn     ">=2"   # pi- (from K_S0), pi- (from h_c)
  end
  .select_photon do
    tdc_emc_start       0
    tdc_emc_end       700
    energyThreshold_b   0.025
    energyThreshold_e   0.050
    angle_to_track     10.0
    nGam               ">=2"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :kaon, against: [:pion, :proton]
    identify :pion, against: [:kaon, :proton]
    nkp  ">=1"
    npim ">=2"
  end
  .secondary_vertex_fit([:pip, :pim]) do
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  end
  .kinematic_fit([:K_S0, :kp, :pim, :pi0]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end

alg
  .note(:ks_fit_quality, "K_S0 secondary vertex fit requires L/deltaL > 2 for decay length significance. K_S0 mass window (0.487, 0.511) GeV/c^2 applied. Best K_S0 pair selected by mass difference. Not fully expressible in DSL.")
  .note(:pi0_mass_window, "pi0 mass window (0.12, 0.15) GeV/c^2 — applied in ROOT analysis after 5C kinematic fit.")
  .note(:"5c_kinematic_fit", "5C kinematic fit (4C energy-momentum + 1C pi0 mass). Chosen candidate has minimum chi2.")
  .note(:signal_yield, "Signal yield from unbinned ML fit to M(K_S0 K+ pi-) distribution. h_c signal shape from signal MC convolved with Gaussian (resolution difference). Peaking background from psi(3686)->gamma chi_c2, chi_c2->K_S0 K+ pi- + fake photon included in fit.")
  .note(:background_shape, "Smooth background modeled by ARGUS function with threshold at 3.551 GeV/c^2. Alternative: second-order Chebyshev polynomial tested for systematics.")
  .note(:peaking_bkg, "Peaking background: psi(3686) -> gamma chi_c2, chi_c2 -> K_S0 K+ pi- with a fake photon. Normalization fixed from PDG BFs. M(K+ pi0) window used to suppress chi_c2 in ROOT analysis.")
  .note(:charge_conjugate, "Charge conjugate mode (anti-h_c -> K_S0 K- pi+) always implied.")
  .note(:photon_quality, "Photon EMC time [0,700] ns; angle to nearest charged track > 10 deg. K+ PID vs pion; remaining tracks assigned as pions.")
  .with_decay_card(decay_card)
  .apply(event_selection)
  .execute_on([psip_data, psip_incMC, exMC])