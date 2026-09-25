# ============================================================================
# ψ(3686) → γ χ_cJ (J = 0, 1, 2) → γ Σ+ anti-p- K_S0
#   Σ+ → p+ π0,  π0 → γγ,  K_S0 → π+ π−
# All three χ_cJ states share the final state: γ p pbar π+ π− γγ
# ============================================================================

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")      # 3.686 GeV real data (ψ(3686))
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # corresponding inclusive MC sample

# Decay card for the combined χ_c0 / χ_c1 / χ_c2 signal chain (EvtGen format)
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    0.3333 gamma chi_c0 P2GC0;
    0.3333 gamma chi_c1 P2GC1;
    0.3334 gamma chi_c2 P2GC2;
    Enddecay

    Decay chi_c0
    1.0000 gamma Sigma+ anti-p- K_S0 PHSP;
    Enddecay

    Decay chi_c1
    1.0000 gamma Sigma+ anti-p- K_S0 PHSP;
    Enddecay

    Decay chi_c2
    1.0000 gamma Sigma+ anti-p- K_S0 PHSP;
    Enddecay

    Decay Sigma+
    1.0000 p+ pi0 PHSP;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# One 100k-event exclusive MC sample covering the combined χ_c0, χ_c1, χ_c2 chain
# (the three states share the same final state and the same event selection)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gammachicJ_sigma_pbar_ks"
  config.related_dataset = psip_data          # associated real dataset
  config.events          = 100_000            # 100k events
  config.decay_card      = decay_card_signal  # decay card of the combined signal chain
  config.cross_section   = :default           # default cross section
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')  # keep the MC configuration

### Event selection (BOSS) ###
alg_name     = "GamChiCJSigmaPbarKS"
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})          # ECMS = 3.686 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
event_selection.select_track {               # Charged track selection
                  cos_theta   0.93           # |cos(theta)| < 0.93
                  Vz          100.0          # |Vz| < 100 cm
                  Vr          10.0           # Vr < 10 cm in the transverse plane
                  nChrp       ">=2"          # At least two positive tracks
                  nChrn       ">=2"          # At least two negative tracks
                }
               .select_photon {              # Photon selection
                  tdc_emc_start     0        # TDC start time
                  tdc_emc_end       14       # TDC end time
                  angle_to_track    10.0     # Min angle to nearest charged track (deg)
                  energyThreshold_b 0.025    # Min energy in the EMC barrel (25 MeV)
                  energyThreshold_e 0.050    # Min energy in the EMC endcap (50 MeV)
                  nGam              ">=3"    # At least three photons
                }
               .pid(method: :probability) {  # Particle identification (no leptons expected)
                  prob_cut 0.001                                                   # PID probability > 0.001
                  identify :proton, against: [:kaon, :pion]                        # p+ and anti-p- vs K, pi
                  nprp ">=1"                                                       # At least one proton
                  nprm ">=1"                                                       # At least one anti-proton
                }
               .remove([:prp <= :chrgp, :prm <= :chrgn])   # Remove identified protons (both charges)
               .assign({:chrgp => :pip, :chrgn => :pim})   # Remaining charged tracks are pi+ / pi-
               .secondary_vertex_fit([:pip, :pim]) {       # Build K_S0 from pi+ pi-
                  build_virtual_particle(:K_S0).by_minimizing_mass_difference  # pick pair closest to K_S0 mass
                  remove_used_particle_from_candidate_list                     # remove used pions from lists
                }
               .kalman_kinematic_fit([:gamma, :gamma]) {   # Kalman (1C) fit: gamma gamma -> pi0
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # constrain to nominal pi0 mass
                  chi2_cut 25                              # chi2 < 25
                  npi0 ">=1"                               # At least one pi0 candidate
                  # NOTE (ROOT level): among the photon pairs the one whose pi0 mass after the
                  # 4C fit is closest to the nominal pi0 mass is chosen; the remaining photon
                  # is the radiative photon.
                }
               # ------------------------------------------------------------------
               # Nominal 4C kinematic fit of the p pbar K_S0 pi0 gamma final state
               # (i.e. gamma p pbar pi+ pi- gamma gamma, built from the reconstructed
               #  K_S0, pi0 and one radiative photon). Only this fit provides the
               #  corrected four-momenta.
               # ------------------------------------------------------------------
               .kinematic_fit([:prp, :prm, :K_S0, :pi0, :gamma]) {
                  nominal                  # nominal fit: corrected four-momenta are saved
                  constrain_four_momentum  # 4C energy-momentum constraint
                  chi2_cut 50              # chi2 < 50
                }
               # Competing 4C hypothesis: no radiative photon (pi0 only) — chi2 stored
               # for the ROOT-level background veto.
               .kinematic_fit([:prp, :prm, :K_S0, :pi0]) {
                  constrain_four_momentum
                }
               # Competing 4C hypothesis: one extra photon (two photons in the fit) — chi2
               # stored for the ROOT-level background veto.
               .kinematic_fit([:prp, :prm, :K_S0, :pi0, :gamma, :gamma]) {
                  constrain_four_momentum
                }

# Anti-Lambda background veto and photon-multiplicity chi2 comparison are applied on the
# fitted quantities (ROOT level) and cannot be expressed in the BOSS selection — kept as a note.
my_Algorithm
  .note(:background_veto, "anti-Lambda background veto |M(pbar pi+) - m_Lambda| > 6 MeV " \
                          "(m_Lambda = 1.11568 GeV) applied on the nominal 4C fit result; " \
                          "additional suppression from comparing the chi2 of the competing 4C " \
                          "fits with different photon multiplicities (no radiative photon / " \
                          "extra photon), both evaluated in the ROOT analysis on the stored chi2 values")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

# Execute on real data, inclusive MC and the signal exclusive MC
root_files = my_Algorithm.execute_on([psip_data, psip_incMC, exMC_signal])