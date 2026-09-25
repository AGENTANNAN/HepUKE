# psi(2S) -> gamma chi_cJ/eta_c(2S), chi_cJ/eta_c(2S) -> p pbar
# Searches for eta_c(2S) -> p pbar and measures chi_cJ -> p pbar BF

psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

decay_card = <<~DECAYCARD
    Decay psi(2S)
    0.2500 gamma chi_c0 PHSP;
    0.2500 gamma chi_c1 PHSP;
    0.2500 gamma chi_c2 PHSP;
    0.2500 gamma eta_c(2S) PHSP;
    Enddecay

    Decay chi_c0
    1.0000 p+ anti-p- PHSP;
    Enddecay

    Decay chi_c1
    1.0000 p+ anti-p- PHSP;
    Enddecay

    Decay chi_c2
    1.0000 p+ anti-p- PHSP;
    Enddecay

    Decay eta_c(2S)
    1.0000 p+ anti-p- PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "psip_gamma_ppbar_signal"
  config.related_dataset = psip_data
  config.events = 100000
  config.decay_card = decay_card
  config.cross_section = :default
end

algorithm = Algorithm.new("PsipToGammaPPbar")
algorithm.set_header(["PsipToGammaPPbarAlg/PsipToGammaPPbar.h"])
algorithm.set_constant({"ECMS" => [:double, 3.686]})

event_selection = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp "==1"
    nChrn "==1"
  end
  .select_photon do
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=1"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp ">=1"
    nprm ">=1"
  end
  .select_isolated_photon do
    angle_to_prp_track 20.0
    angle_to_prm_track 30.0
    nGam ">=1"
  end
  .kinematic_fit([:gamma, :prp, :prm]) do
    nominal
    vertex_fit([1, 2])
    constrain_four_momentum
    chi2_cut 60
  end
  .kinematic_fit([:prp, :prm]) do
    constrain_four_momentum
  end

algorithm
  .note(:chi2_4c_comparison, "chi2_4c(gamma p pbar) < chi2_4c(p pbar) required to suppress psi(2S)->p pbar with fake/FSR photon; comparison applied in ROOT using stored chi2_4c values from nominal and competing kinematic fits")
  .note(:three_c_fit, "3C kinematic fit (photon momentum magnitude allowed to float) used to obtain M_ppbar_3C for background separation in ROOT fit; not expressible in BOSS DSL")
  .with_decay_card(decay_card)
  .apply(event_selection)

algorithm.execute_on([psip_data, psip_incMC, exMC_signal])