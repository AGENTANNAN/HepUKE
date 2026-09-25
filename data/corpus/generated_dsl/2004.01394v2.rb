### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")      # ψ(3686) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # corresponding inclusive MC

# Decay card: ψ(3686) → γ χ_cJ (J = 0, 1, 2), χ_cJ → Σ− Σ̄+, Σ− → n π−, Σ̄+ → n̄ π+
# (all three χ_cJ modes share the same final state γ n̄ π+ π− with a missing neutron,
#  so a single merged card / single algorithm serves the whole χ_cJ family)
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    0.3333 gamma chi_c0 P2GC0;
    0.3333 gamma chi_c1 P2GC1;
    0.3334 gamma chi_c2 P2GC2;
    Enddecay

    Decay chi_c0
    1.0000 Sigma- anti-Sigma+ PHSP;
    Enddecay

    Decay chi_c1
    1.0000 Sigma- anti-Sigma+ PHSP;
    Enddecay

    Decay chi_c2
    1.0000 Sigma- anti-Sigma+ PHSP;
    Enddecay

    Decay Sigma-
    1.0000 n0 pi- PHSP;
    Enddecay

    Decay anti-Sigma+
    1.0000 anti-n0 pi+ PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC: 100k events of ψ(3686) → γ χ_cJ → γ Σ− Σ̄+
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gamchicJ_ssbar"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "PsipToGamChicJSigmaSigma"
psi_alg = Algorithm.new(alg_name)
psi_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
       .set_constant({"ECMS" => [:double, 3.686]})    # 3.686 GeV ψ(3686) centre-of-mass energy
       .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
  .select_track {                       # charged track (π+ π−) quality selection
      cos_theta 0.93                    # |cosθ| < 0.93
      Vz        30.0                    # |Vz| < 30 cm
      Vr        10.0                    # Vr < 10 cm
      nChrp     "==1"                   # exactly one positively charged track
      nChrn     "==1"                   # exactly one negatively charged track
  }
  .select_photon {                      # photon selection
      tdc_emc_start     0               # EMC TDC start
      tdc_emc_end       14              # EMC TDC end
      energyThreshold_b 0.025           # E > 25 MeV in the barrel
      energyThreshold_e 0.050           # E > 50 MeV in the endcap
      angle_to_track    10.0            # > 10° away from any charged track
      nGam              ">=1"           # at least one photon
  }
  .pid(method: :probability) {          # pion identification (probability method)
      prob_cut 0.001                    # PID probability > 0.001
      identify :pion, against: [:kaon, :proton]   # π+ and π−, separated from K and p
      npip     "==1"                    # one π+
      npim     "==1"                    # one π−
  }
  # Nominal kinematic fit: ψ(3686) → γ n̄ π+ π− with the neutron missing.
  # Σ̄+ is reconstructed from n̄ π+ by constraining the n̄π+ invariant mass to the
  # Σ̄+ nominal value; events are selected by the best (minimum-χ²) combination.
  .kinematic_fit([:gamma, :n_bar, :pip, :pim]) {
      nominal                                       # nominal fit — fitted four-momenta are stored
      miss_track_of :n                              # the neutron from Σ− → n π− is not detected
      constrain_four_momentum                       # 4C energy-momentum constraint to ECMS
      invariant_mass_of(:n_bar, :pip).constrain_to_nominal_mass_of(:Sigma_bar_plus)
      # K_S0 veto: reject π+π− combinations within 10 MeV of the K_S0 mass
      invariant_mass_of(:pip, :pim).out_of(0.4876, 0.5076)
      chi2_cut 20                                   # χ² < 20 (default best combination by min χ²)
  }
  # Competing-hypothesis fit (Σ− Σ̄+ without the radiative photon) — no chi2_cut and no
  # nominal: only the competing χ² is stored, and the veto
  # (chi2_4c_γΣΣ < chi2_4c_ΣΣ) is applied at the ROOT level.
  .kinematic_fit([:n_bar, :pip, :pim]) {
      miss_track_of :n
      constrain_four_momentum
  }

# BOSS-side procedures that cannot be expressed with the current DSL
psi_alg
  .note(:nbar_identification, "antineutron (n̄) taken as the most energetic EMC shower with \
    deposited energy 0.2-2.0 GeV, second moment > 20 cm^2 and more than 20 EMC hits inside a \
    40 degree cone; the neutron (n) from Sigma- -> n pi- is not reconstructed and is treated \
    as missing in the kinematic fit")
  .note(:background_veto, "before the kinematic fit: photon pairs with invariant mass within \
    12 MeV of the pi0 mass are vetoed, and events whose pi+pi- recoil mass is within 10 MeV of \
    the J/psi mass are vetoed (recoil_mass_of has no DSL equivalent, applied in the generated code)")
  .note(:radiative_photon_selection, "the radiative photon is required to have energy > 80 MeV, \
    to be more than 40 degrees away from the antineutron shower and more than 10 degrees from \
    any charged track")

# Generate the algorithm for the decay card and run on data, inclusive MC and signal MC
psi_alg.with_decay_card(decay_card_signal).apply(event_selection)
root_files = psi_alg.execute_on([psip_data, psip_incMC, exMC_signal])