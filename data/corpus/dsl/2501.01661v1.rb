# Paper 2501.01661v1 — Search for η_c(2S)→p p̄ K⁺ K⁻ and measurement of χ_cJ→p p̄ K⁺ K⁻
# in ψ(3686) radiative decays with (2712.4±14.3)×10⁶ ψ(3686) events

### Dataset preparation ###
psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Decay card for the algorithm (shared final state: γ p p̄ K⁺ K⁻)
decay_card = <<~DECAYCARD
  Decay psi(2S)
  1.0000 gamma p+ anti-p- K+ K- PHSP;
  Enddecay
  End
DECAYCARD

# Exclusive MC — ψ(3686) → γ η_c(2S) → γ p p̄ K⁺ K⁻
exMC_etac2s = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_etac2s_ppKK"
  config.related_dataset = psip_data
  config.events = 500_000
  config.decay_card = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma eta_c(2S) VSP_PWAVE 0 0 0 0 0 0 0;
    Enddecay
    Decay eta_c(2S)
    1.0000 p+ anti-p- K+ K- PHSP;
    Enddecay
    End
  DECAYCARD
  config.cross_section = :default
end

# Exclusive MC — ψ(3686) → γ χ_c0 → γ p p̄ K⁺ K⁻
exMC_chic0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_chic0_ppKK"
  config.related_dataset = psip_data
  config.events = 500_000
  config.decay_card = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c0 VSP_PWAVE 0 0 0 0 0 0 0;
    Enddecay
    Decay chi_c0
    1.0000 p+ anti-p- K+ K- PHSP;
    Enddecay
    End
  DECAYCARD
  config.cross_section = :default
end

# Exclusive MC — ψ(3686) → γ χ_c1 → γ p p̄ K⁺ K⁻
exMC_chic1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_chic1_ppKK"
  config.related_dataset = psip_data
  config.events = 500_000
  config.decay_card = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c1 VSP_PWAVE 0 0 0 0 0 0 0;
    Enddecay
    Decay chi_c1
    1.0000 p+ anti-p- K+ K- PHSP;
    Enddecay
    End
  DECAYCARD
  config.cross_section = :default
end

# Exclusive MC — ψ(3686) → γ χ_c2 → γ p p̄ K⁺ K⁻
exMC_chic2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_chic2_ppKK"
  config.related_dataset = psip_data
  config.events = 500_000
  config.decay_card = <<~DECAYCARD
    Decay psi(2S)
    1.0000 gamma chi_c2 VSP_PWAVE 0 0 0 0 0 0 0;
    Enddecay
    Decay chi_c2
    1.0000 p+ anti-p- K+ K- PHSP;
    Enddecay
    End
  DECAYCARD
  config.cross_section = :default
end

# Exclusive MC — ψ(3686) → p p̄ K⁺ K⁻ background (no radiative photon)
exMC_bkg = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "bkg_psip_ppKK"
  config.related_dataset = psip_data
  config.events = 200_000
  config.decay_card = <<~DECAYCARD
    Decay psi(2S)
    1.0000 p+ anti-p- K+ K- PHSP;
    Enddecay
    End
  DECAYCARD
  config.cross_section = :default
end

### Event selection (BOSS) ###
algorithm = Algorithm.new("PsipRadiativePPKK")

algorithm.set_header(["PsipRadiativePPKKAlg/PsipRadiativePPKK.h"])
          .set_constant("ECMS" => [:double, 3.686])

event_selection = Selection.new
  # Charged track selection: 4 tracks (2 pos, 2 neg), net charge zero
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp "==2"
    nChrn "==2"
    nNet "==0"
  end
  # Photon selection: at least 1 photon, EMC barrel >25 MeV, endcap >50 MeV, TDC [0,700] ns
  .select_photon do
    tdc_emc_start 0
    tdc_emc_end 14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=1"
  end
  # Particle ID: first identify proton / anti-proton
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp "==1"
    nprm "==1"
  end
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  # Identify kaons from remaining tracks
  .pid(method: :probability) do
    prob_cut 0.001
    identify :kaon, against: [:pion]
    nkp "==1"
    nkm "==1"
  end
  # 3C kinematic fit (nominal) with vertex fit on charged tracks; photon energy floats
  .kinematic_fit([:prp, :prm, :kp, :km, :gamma]) do
    nominal
    constrain_three_momentum
    vertex_fit([0, 1, 2, 3])
    chi2_cut 200
  end

algorithm
  .note(:fake_photon_resampling, "fake photon (theta_gamma_pbar vs E_gamma) distribution re-sampled per 2D bin ratio r=nData/nMC; photons re-processed through 4C kinematic fit after resampling")
  .note(:fsr_correction, "FSR fraction ratio R_FSR = f_FSR_data / f_FSR_MC = 2.38±0.90 applied to psi(3686)→p pbar K+ K- exclusive MC for background PDF shape; determined from control sample psi(3686)→γ chi_c0→γ gamma_FSR p pbar K+ K-")
  .note(:helix_correction, "helix parameter correction applied to charged tracks before kinematic fit; efficiency difference between with/without helix correction taken as systematic")
  .note(:efficiency_curve, "detection efficiency curve obtained via Gaussian Process Regression on discrete MC efficiency points; 100 alternative GPR curves used for systematic uncertainty estimation")
  .note(:angular_distribution, "photon angular distribution 1+α cos²θ_γ applied in MC; α=1 for η_c(2S), α=1 for χ_c0, α=-1/3 for χ_c1 (E1 dominant), α≈1/12 for χ_c2; difference from M2/E3 mixing taken as systematic")
  .note(:damping_function, "M1/E1 damping function F(E_γ) used in signal PDF lineshape; alternative CLEO exponential damping exp(-E_γ²/8β²) with β=0.12 GeV tested for systematic")
  .note(:resolution_smearing, "MC-data resolution difference (δm₂, σ₂) for η_c(2S) linearly extrapolated from χ_cJ values; varied by ±1σ for systematic")
  .with_decay_card(decay_card)
  .apply(event_selection)

root_files = algorithm.execute_on([psip_data, psip_incMC, exMC_etac2s, exMC_chic0, exMC_chic1, exMC_chic2, exMC_bkg])