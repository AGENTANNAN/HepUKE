### Dataset description ###
jpsi_data = DatasetManager.real_data.find("708_3097")       # J/ψ real data at 3.097 GeV
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")   # Corresponding inclusive MC sample

# Mode I (phase space): J/ψ → ω η′ π⁺π⁻ ; ω → π⁺π⁻π⁰, π⁰ → γγ ; η′ → η π⁺π⁻, η → γγ
decay_card_phsp = <<~DECAYCARD
    Decay J/psi
    1.0000 omega eta' pi+ pi-    PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0            OMEGA_DALITZ;
    Enddecay

    Decay eta'
    1.000 eta pi+ pi-            PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma            PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma            PHSP;
    Enddecay

    End
DECAYCARD

# Mode II: J/ψ → ω X(1835) ; X(1835) → η′ π⁺π⁻ ; ω → π⁺π⁻π⁰, π⁰ → γγ ; η′ → η π⁺π⁻, η → γγ
decay_card_x1835 = <<~DECAYCARD
    Decay J/psi
    1.0000 omega X(1835)         PHSP;
    Enddecay

    Decay X(1835)
    1.000 eta' pi+ pi-           PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0            OMEGA_DALITZ;
    Enddecay

    Decay eta'
    1.000 eta pi+ pi-            PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma            PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma            PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive MC samples for each of the two signal modes
exMC_phsp = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "jpsi_omega_etap_pipi_phsp"   # Mode I phase space
  config.related_dataset = jpsi_data
  config.events         = 100_000
  config.decay_card     = decay_card_phsp
  config.cross_section  = :default
end

exMC_x1835 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "jpsi_omega_x1835"            # Mode II X(1835)
  config.related_dataset = jpsi_data
  config.events         = 100_000
  config.decay_card     = decay_card_x1835
  config.cross_section  = :default
end

### Event selection (BOSS) ###
alg_name = "OmegaEtapX1835"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})     # CMS energy = 3.097 GeV (J/ψ)
            .set_alias({"std::vector<double>" => "Vdouble"})

# Both modes share the same 3π⁺3π⁻4γ final state and identical selection → one algorithm
event_selection = Selection.new
  .select_track {                       # Charged track selection
    cos_theta 0.93                      # |cosθ| < 0.93
    Vz 100.0                            # |Vz| < 100 cm
    Vr 10.0                             # Vr < 10 mm in transverse plane
    nChrp "==3"                         # exactly three positive tracks
    nChrn "==3"                         # exactly three negative tracks
    nNet  "==0"                         # net charge zero
  }
  .select_photon {                      # Photon selection
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 5.0                  # at least 5° from any charged track
    energyThreshold_b 0.025             # barrel energy threshold (25 MeV)
    energyThreshold_e 0.050             # endcap energy threshold (50 MeV)
    nGam ">=4"                          # at least four photons
  }
  .pid(method: :probability) {          # Particle identification
    prob_cut 0.001                      # PID probability > 0.001
    identify :pion, against: [:kaon]    # π⁺ and π⁻ identified against kaon
    npip "==3"
    npim "==3"
  }
  .remove([:pip <= :chrgp, :pim <= :chrgn])   # remove pion candidates from charged-track lists
  .kalman_kinematic_fit([:gamma, :gamma]) {   # Reconstruct π⁰ from γγ (1-C mass constraint)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"                                # at least one π⁰
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # Reconstruct η from γγ (1-C mass constraint)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
    chi2_cut 25
    neta ">=1"                                # at least one η
  }
  .kinematic_fit([:pip, :pip, :pip, :pim, :pim, :pim, :pi0, :eta]) {  # Nominal 4C fit
    nominal                                   # nominal fit — corrected four-momenta saved
    constrain_four_momentum                   # 4-momentum conservation to CMS
    chi2_cut 60                               # χ² < 60
    invariant_mass_of(:pip, :pim, :pi0).within(0.7607, 0.8047)  # ω window: |M(π⁺π⁻π⁰) − m_ω| < 22 MeV
    invariant_mass_of(:eta, :pip, :pim).within(0.9458, 0.9698)  # η′ window: |M(ηπ⁺π⁻) − m_η′| < 12 MeV
  }

# Best ω/η′ combination is a mass-difference ranking; the kinematic fit selects by smallest χ².
my_algorithm
  .note(:best_combination, "Best ω/η′ combination chosen by minimising sqrt((M(π+π-π0)-m_ω)^2 + (M(ηπ+π-)-m_η')^2); not directly expressible in DSL (the 4C fit auto-selects by smallest χ²), estimated/verified by re-running the selection")
  .with_decay_card(decay_card_phsp).apply(event_selection)

# Execute on J/ψ real data, inclusive MC, and both exclusive MC samples
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_phsp, exMC_x1835])