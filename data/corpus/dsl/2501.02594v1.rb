# Paper 2501.02594v1 — First observation of ψ(3686)→K⁻Λ(1520)Ξ̄⁺ + c.c.
# with (2712.4±14.3)×10⁶ ψ(3686) events
# Decay chain: ψ(3686)→K⁻ Λ(1520) Ξ̄⁺, Λ(1520)→p K⁻, Ξ̄⁺→Λ̄ π⁺, Λ̄→p̄ π⁺
# Final state: p p̄ K⁻ K⁻ π⁺ π⁺

### Dataset preparation ###
psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

decay_card = <<~DECAYCARD
  Decay psi(2S)
  1.0000 K- Lambda(1520) anti-Xi- PHSP;
  Enddecay

  Decay Lambda(1520)
  1.0000 p+ K- PHSP;
  Enddecay

  Decay anti-Xi-
  1.0000 anti-Lambda0 pi+ PHSP;
  Enddecay

  Decay anti-Lambda0
  1.0000 anti-p- pi+ HypWK;
  Enddecay

  End
DECAYCARD

# Singal exclusive MC
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_K_L1520_Xi"
  config.related_dataset = psip_data
  config.events = 500_000
  config.decay_card = decay_card
  config.cross_section = :default
end

### Event selection (BOSS) ###
algorithm = Algorithm.new("PsipKL1520XiBar")

algorithm.set_header(["PsipKL1520XiBarAlg/PsipKL1520XiBar.h"])
          .set_constant("ECMS" => [:double, 3.686])

event_selection = Selection.new
  # Charged track selection: 6 tracks (3 positive, 3 negative), net zero charge
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp "==3"
    nChrn "==3"
    nNet "==0"
  end
  # PID: identify proton, anti-proton, two K⁻, two π⁺
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp "==1"
    nprm "==1"
  end
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  # Identify kaons
  .pid(method: :probability) do
    prob_cut 0.001
    identify :kaon, against: [:pion]
    nkm "==2"
  end
  .remove([:km <= :chrgn])
  # Remaining tracks are pions; identify two π⁺
  .pid(method: :probability) do
    prob_cut 0.001
    identify :pion, against: [:kaon]
    npip "==2"
  end
  # Reconstruct Λ̄ → p̄ π⁺ via secondary vertex fit (long-lived V0)
  .secondary_vertex_fit([:prm, :pip]) do
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  # 4C kinematic fit over all 6 final-state charged particles
  .kinematic_fit([:prp, :prm, :km, :km, :pip, :pip]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end

algorithm
  .note(:recoil_mass_Xi_tag, "Ξ̄⁺ tagged via recoil mass RM(p K₁⁻ K₂⁻); signal region 1.31–1.34 GeV/c²; sidebands 1.242–1.272 and 1.372–1.402 GeV/c²; normalization factor f_sideband=0.47±0.01 from fitted background function")
  .note(:k1_k2_assignment, "K₁⁻ (from Λ(1520)) assigned by requiring M(p K₁⁻) closer to nominal Λ(1520) mass than M(p K₂⁻); two entries per event in combined M(p K⁻) fit")
  .note(:signal_region, "Λ(1520) signal region: 1.50–1.54 GeV/c²; sidebands: 1.43–1.47 and 1.57–1.61 GeV/c²")
  .note(:continuum_scale, "continuum background from √s=3.773 GeV data scaled by factor f_c = (L_3686/L_3773)×(σ_3686/σ_3773) with σ ∝ 1/s; f_c=1.386")
  .note(:signal_shape, "Λ(1520) signal shape: MC shape convolved with Gaussian for data-MC resolution difference; free parameters in unbinned maximum likelihood fit to combined M(p K⁻)")
  .note(:background_shape, "non-Λ(1520) background: ψ(3686)→p K⁻ K⁻ Ξ̄⁺ four-body MC shape with floating yield; non-Ξ background shape/yield fixed from Ξ sidebands after smoothing+normalization")
  .note(:track_ip_requirement, "at least 3 tracks required to originate from IP (bachelor K⁻ and Λ(1520)→p K⁻ tracks); standard V_r<1 cm, |V_z|<10 cm on all 6 tracks")
  .note(:pid_likelihood_sum, "best PID combination selected by maximizing sum of PID likelihoods for p p̄ K⁻ K⁻ π⁺ π⁺ combination")
  .with_decay_card(decay_card)
  .apply(event_selection)

root_files = algorithm.execute_on([psip_data, psip_incMC, exMC_signal])