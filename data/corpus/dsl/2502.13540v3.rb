# ψ(3686) → γ K_S^0 K_S^0 amplitude analysis
# Paper: 2502.13540v3

DatasetManager.load_real_data('assets/BES3_dataset.md')
DatasetManager.load_inclusive_mc('assets/BES3_incMC.md')

psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

decay_card_signal = <<~DECAYCARD
  Decay psi(2S)
  1.0000 gamma K_S0 anti-K_S0 PHSP;
  Enddecay
  Decay K_S0
  1.0000 pi+ pi- PHSP;
  Enddecay
  Decay anti-K_S0
  1.0000 pi+ pi- PHSP;
  Enddecay
  End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_psip_gamma_ks_ks"
  config.related_dataset = psip_data
  config.events          = 500_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

algorithm = Algorithm.new("PsipGammaKsKsAnalysis")
algorithm
  .set_header(["PsipGammaKsKsAnalysisAlg/PsipGammaKsKsAnalysis.h"])
  .set_constant(ECMS: 3.686)

event_selection = Selection.new
  .select_track do
    cos_theta   0.93
    Vz          20.0
    Vr          1.0
    nChrp       "==2"
    nChrn       "==2"
    nNet        "==0"
  end
  .assign({chrgp: :pip, chrgn: :pim})
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=1"
  end
  .secondary_vertex_fit([:pip, :pim]) do
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .secondary_vertex_fit([:pip, :pim]) do
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .kinematic_fit([:gamma, :pip, :pim, :pip, :pim]) do
    constrain_four_momentum
    chi2_cut 50
    nominal
  end

algorithm
  .note(:ks_mass_window, "K_S^0 mass window: |M(π+π-) - M(K_S^0)| < 12 MeV/c^2, corresponding to 2.5σ")
  .note(:ks_decay_length, "Decay length of each K_S^0 candidate > 2σ of its vertex fit resolution")
  .note(:ks_count, "Exactly two K_S^0 candidates required per event")
  .note(:mksks_cut, "M(K_S^0 K_S^0) < 2.8 GeV/c^2 to suppress χc0,2 backgrounds")
  .note(:continuum_suppression, "Continuum background estimated from 3650 MeV data sample (401 pb⁻¹), found negligible")
  .note(:ks_indistinguishable, "Two K_S^0 candidates are mixed randomly as they are indistinguishable")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)
  .execute_on([psip_data, psip_incMC, exMC_signal])