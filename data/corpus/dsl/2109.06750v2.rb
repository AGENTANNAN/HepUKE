# BESIII DSL: First observation of psi(3686) -> omega K_S0 K_S0
# arXiv: 2109.06750v2
# omega -> pi+ pi- pi0, K_S0 -> pi+ pi-, pi0 -> gamma gamma
# Uses psi(3686) data at 3.686 GeV

### Dataset preparation ###
psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

decay_card = <<~DECAYCARD
    Decay psi(3686)
    1.000 omega K_S0 K_S0        PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0            OMEGA_DALITZ;
    Enddecay

    Decay K_S0
    1.000 pi+ pi-                PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma            PHSP;
    Enddecay

    End
DECAYCARD

exMC = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "psi3686_omega_ks_ks_exclusive_mc"
  config.related_dataset = psip_data
  config.events = 500_000
  config.decay_card = decay_card
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg = Algorithm.new("Psi3686OmegaKsKs")
alg.set_header(["Psi3686OmegaKsKsAlg/Psi3686OmegaKsKs.h"])
alg.set_constant({"ECMS" => [:double, 3.686]})

event_selection = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=3"
    nChrn ">=3"
    nNet "==0"
  end
  .select_photon do
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=2"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
    npip ">=3"
    npim ">=3"
  end
  # Reconstruct pi0 from photon pairs
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  end
  # Reconstruct K_S0 -> pi+ pi- (first)
  .secondary_vertex_fit([:pip, :pim]) do
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  # Reconstruct K_S0 -> pi+ pi- (second)
  .secondary_vertex_fit([:pip, :pim]) do
    build_virtual_particle(:K_S0_2).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  # 5C kinematic fit: 4C + pi0 mass constraint
  .kinematic_fit([:pip, :pim, :K_S0, :K_S0_2, :pi0]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end

# Note: the DSL doesn't support building K_S0_2 (second virtual particle of same type),
# so we use :K_S0_2 and note the actual implementation
alg.note(:ks_window, "K_S0 mass window (0.486, 0.510) GeV/c^2 (3 sigma); decay length > 2 sigma from IP")
  .note(:omega_window, "omega mass window (0.777, 0.807) GeV/c^2 used in final selection")
  .note(:pi0_mass_window, "pi0 mass window (0.115, 0.150) GeV/c^2")
  .note(:ks_sideband, "2D K_S0 K_S0 sideband method to subtract combinatorial background; sideband regions defined in M(pi+pi-)(1) vs M(pi+pi-)(2) plane")
  .note(:body3_generator, "BODY3 data-driven generator used for signal MC modeling; Dalitz plot of M(omega K_S0) vs M(omega K_S0) from efficiency-corrected data as input")
  .note(:kinematic_fit_5c, "5C kinematic fit (4C + pi0 mass constraint) with chi2_5C < 80; candidate selection by minimizing chi2_5C + chi2(K_S0_1) + chi2(K_S0_2)")
  .note(:continuum_subtraction, "Continuum contribution from sqrt(s)=3.65 GeV data (43.88 pb^-1) subtracted using scale factor f_c; cross section ~ 1/s dependence assumed")
  .note(:branching_fraction, "B(psi(3686)->omega K_S0 K_S0) = (7.04 +/- 0.39 +/- 0.37)e-5; first observation")
  .with_decay_card(decay_card)
  .apply(event_selection)

root_files = alg.execute_on([psip_data, psip_incMC, exMC])