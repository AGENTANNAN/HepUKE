# Paper: 2209.11175v2 — PWA of J/ψ → γK_S⁰K_S⁰π⁰
# Ordinary analysis at J/ψ (3.097 GeV), 10.09×10⁹ J/ψ events

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Exclusive MC: J/ψ → γ K_S⁰ K_S⁰ π⁰ phase space
decay_card = <<~DECAYCARD
Decay J/psi
1 gamma K_S0 K_S0 pi0 PHSP;
Enddecay
Decay K_S0
1 pi+ pi- PHSP;
Enddecay
Decay pi0
1 gamma gamma PHSP;
Enddecay
End
DECAYCARD

sig_mc = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_jpsi_gamma_ks_ks_pi0"
  config.related_dataset = jpsi_data
  config.events          = 500_000
  config.decay_card      = decay_card
  config.cross_section   = :default
end

alg = Algorithm.new("JpsiToGammaKsKsPi0")
alg.set_header(["JpsiToGammaKsKsPi0Alg/JpsiToGammaKsKsPi0.h"])
    .set_constant({ "ECMS" => [:double, 3.097] })
    .with_decay_card(decay_card)

event_selection = Selection.new
  .select_track do
    cos_theta  0.93
    Vz         20.0
    Vr         1.0
    nChrp      "==2"
    nChrn      "==2"
    nNet       "==0"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip ">=2"
    npim ">=2"
  end
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=3"
  end
  .secondary_vertex_fit([:pip, :pim]) do
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
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
  .kinematic_fit([:gamma, :K_S0, :K_S0, :pi0]) do
    constrain_four_momentum
    invariant_mass_of(:pi0).constrain_to_nominal_mass_of(:pi0)
    invariant_mass_of(:K_S0).constrain_to_nominal_mass_of(:K_S0)
    invariant_mass_of(:K_S0).constrain_to_nominal_mass_of(:K_S0)
    chi2_cut 200
    nominal
  end

alg.note(:track_quality, "Paper: |cosθ|<0.93, Vz<20cm (looser due to K_S⁰ lifetime), Vr<1cm. All tracks under pion hypothesis. DSL uses Vz:20.0.")
alg.note(:ks_secondary_vertex, "Paper: K_S⁰ secondary vertex fit + L/σ_L>2 on π⁺π⁻. DSL secondary_vertex_fit handles vertex+mass; L/σ_L cut is ROOT-level.")
alg.note(:photon_quality, "Paper: E>25MeV barrel, E>50MeV endcap, angle to track>10°, EMC time [0,700]ns. DSL defaults match.")
alg.note(:pi0_window, "Paper: |M(γγ)-M(π⁰)|<20MeV. ROOT-level cut; kalman_kinematic_fit chi2_cut 25 for loose pre-selection.")
alg.note(:ks_window, "Paper: both K_S⁰ |M(π⁺π⁻)-M(K_S⁰)|<12MeV. ROOT-level cut.")
alg.note(:eta_veto, "Paper: η veto — reject events with any |M(γγ)-M(η)|<30MeV. ROOT-level cut.")
alg.note(:omega_veto, "Paper: ω veto — reject events with |M(γπ⁰)-M(ω)|<40MeV. ROOT-level cut.")
alg.note(:pwa_mass_cut, "Paper: M(K_S⁰K_S⁰π⁰)<1.6GeV/c² for PWA. ROOT-level cut.")
alg.note(:four_c_combination, "Paper: 4C kinematic fit χ²<40 selects best γγπ⁰K_S⁰K_S⁰ combination. DSL nominal 7C fit iterates all combinations; tight 4C cut applied in ROOT.")
alg.note(:seven_c_fit, "Paper: 7C = 4C + π⁰ mass + K_S⁰ mass (×2). DSL emits 7C (constrain_four_momentum + 3 invariant_mass_of constraints).")

alg.with_decay_card(decay_card).apply(event_selection)
alg.execute_on([jpsi_data, jpsi_incMC, sig_mc])