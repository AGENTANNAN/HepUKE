# =============================================================================
# ψ(3686) → Σ⁻ (→ n π⁻) Σ̄⁺ (→ n̄ π⁺)  —  measurement of α_{Σ⁻}
# BOSS part only: dataset preparation + event selection up to and including the
# final kinematic fit.  The anti-neutron is reconstructed from EMC showers,
# the neutron is left undetected (missing) and enters only through the fit.
# =============================================================================

### Dataset description ###
psip_data   = DatasetManager.real_data.find("709_3686")     # 448.1M ψ(3686) events at 3.686 GeV
psip_incMC  = DatasetManager.inclusive_mc.find("709_3686")  # matching inclusive MC at 3.686 GeV
offres_3650 = DatasetManager.real_data.find("709_3650")     # 3.65 GeV off-resonance data (non-ψ(3686) background)

# Decay card (EvtGen) for the signal process ψ(3686) → Σ⁻ Σ̄⁺ with the two-body Σ decays
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000 Sigma- anti-Sigma+      PHSP;
    Enddecay

    Decay Sigma-
    1.0000 n0 pi-                  PHSP;
    Enddecay

    Decay anti-Sigma+
    1.0000 anti-n0 pi+             PHSP;
    Enddecay

    End
DECAYCARD

# 500k exclusive signal MC events
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "psip_to_sigma_sigmabar_exclusive_mc"
  config.related_dataset = psip_data
  config.events          = 500_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "PsipToSigmaSigmaBar"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})          # √s = 3.686 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
  .select_track {              # exactly two charged tracks with zero net charge
    cos_theta 0.93             # |cosθ| < 0.93
    Vz        30.0             # |Vz| < 30 cm
    Vr        10.0             # Vr < 10 cm
    nChrp     "==1"            # one positive track  (π⁺)
    nChrn     "==1"            # one negative track  (π⁻)
    nNet      "==0"            # net charge zero
  }
  .select_photon {             # EMC showers kept as anti-neutron candidates
    tdc_emc_start     0        # EMC time window [0, 700] ns
    tdc_emc_end       14
    angle_to_track    10.0     # > 10° away from the nearest charged track
    energyThreshold_b 0.6      # E > 600 MeV (barrel)
    energyThreshold_e 0.6      # E > 600 MeV (endcap)
    nGam              ">=1"    # at least one anti-neutron candidate
  }
  .pid(method: :probability) {                     # pion identification (probability method)
    prob_cut 0.001                                 # L > 0.001
    identify :pion, against: [:kaon, :proton]      # L(π) > L(K) and L(π) > L(p)
    npip "==1"                                     # exactly one π⁺
    npim "==1"                                     # exactly one π⁻
  }
  # BOSS 4C kinematic fit of π⁺π⁻ with four-momentum conservation, χ² < 200
  .kinematic_fit([:pip, :pim]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }
  # Final fit: neutron treated as missing, M(n̄π⁺) constrained to the Σ⁺ mass, χ² < 50
  .kinematic_fit([:pip, :pim, :n_bar]) {
    miss_track_of(:n0)
    invariant_mass_of(:n_bar, :pip).constrain_to_nominal_mass_of(:"Sigma+")
    chi2_cut 50
  }

# BOSS-side procedures that cannot be expressed with the DSL primitives
my_algorithm
  .note(:background_veto, "ψ(3686) → π⁺π⁻ J/ψ with J/ψ → n n̄ is vetoed by requiring
    the π⁺π⁻ recoil mass (evaluated against the CMS four-momentum) to be below
    2.9 GeV/c²; the DSL has no recoil-mass primitive, so the veto is applied
    directly in the generated selection code.")
  .note(:nbar_shower_selection, "anti-neutron candidates are EMC showers with lateral
    moment > 20 and only the most energetic candidate per event is retained;
    neither the lateral-moment shower-shape cut nor the most-energetic-candidate
    retention is expressible with select_photon and both are applied in the
    generated code.")
  .note(:final_fit, "in the final kinematic fit the (missing) neutron three-momentum
    and the anti-neutron angles are left as free parameters in addition to the
    M(n̄π⁺) = M(Σ⁺) constraint; only the mass constraint is expressible in DSL
    kinematic_fit syntax, the remaining free-parameter parameterization is
    implemented in the generated code.")
  .note(:peaking_background, "peaking γχ_cJ (J = 0,1,2) / γη_c contributions and the
    non-ψ(3686) continuum are accounted for with the inclusive MC and with the
    3.65 GeV off-resonance data (709_3650) processed through the same selection.")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

# Execute on real data, inclusive MC, off-resonance data and the signal MC
root_files = my_algorithm.execute_on([psip_data, psip_incMC, offres_3650, exMC_signal])