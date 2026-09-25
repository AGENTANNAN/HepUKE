# =============================================================================
# e+e- -> Ds*+ Ds-, D*0 anti-D0, D*+ D-   at sqrt(s) = 4.178 GeV
# Partial reconstruction: reconstruct exactly ONE D*(s) per event; the other
# D is recovered from the recoil four-momentum.
# =============================================================================

### ---------------------------- Dataset preparation ---------------------------- ###
data_4180  = DatasetManager.real_data.find("703_4180")     # real data, 3.19 fb^-1 @ 4.178 GeV
incMC_4180 = DatasetManager.inclusive_mc.find("703_4180")  # corresponding inclusive MC

# Decay cards (EvtGen syntax). The D on the "other side" is present in the card
# but is NOT reconstructed - it is the recoil.

# Mode I: Ds*+ Ds-  (Ds*+ -> gamma Ds+, Ds+ -> K_S0 K+ ; Ds- present but not reconstructed)
decay_card_ds = <<~DECAYCARD
    Decay psi(4260)
    1.000 D_s*+ D_s- PHSP;
    Enddecay

    Decay D_s*+
    1.000 gamma D_s+ PHSP;
    Enddecay

    Decay D_s+
    1.000 K_S0 K+ PHSP;
    Enddecay

    Decay D_s-
    1.000 K+ K- pi- PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Mode II: D*0 anti-D0  (D*0 -> D0 pi0, D0 -> K- pi+ ; anti-D0 -> K+ pi- present but not reconstructed)
decay_card_d0 = <<~DECAYCARD
    Decay psi(4260)
    1.000 D*0 anti-D0 PHSP;
    Enddecay

    Decay D*0
    1.000 D0 pi0 PHSP;
    Enddecay

    Decay D0
    1.000 K- pi+ PHSP;
    Enddecay

    Decay anti-D0
    1.000 K+ pi- PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Mode III: D*+ D-  (D*+ -> D+ pi0, D+ -> K- pi+ pi+ ; D- -> K+ pi- pi- present but not reconstructed)
decay_card_dp = <<~DECAYCARD
    Decay psi(4260)
    1.000 D*+ D- PHSP;
    Enddecay

    Decay D*+
    1.000 D+ pi0 PHSP;
    Enddecay

    Decay D+
    1.000 K- pi+ pi+ PHSP;
    Enddecay

    Decay D-
    1.000 K+ pi- pi- PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC samples: 500k events each
exMC_ds = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4180_dsstar_ds"
  config.related_dataset = data_4180
  config.events          = 500000
  config.decay_card      = decay_card_ds
  config.cross_section   = :default
end

exMC_d0 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4180_dstar0_d0bar"
  config.related_dataset = data_4180
  config.events          = 500000
  config.decay_card      = decay_card_d0
  config.cross_section   = :default
end

exMC_dp = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4180_dstar_d"
  config.related_dataset = data_4180
  config.events          = 500000
  config.decay_card      = decay_card_dp
  config.cross_section   = :default
end

### ---------------------------- Event selection (BOSS) ---------------------------- ###

# ---------------------------------------------------------------------------
# Mode I : e+e- -> Ds*+ Ds- , reconstruct Ds*+ -> gamma Ds+ -> gamma K_S0 K+
# ---------------------------------------------------------------------------
alg_ds = Algorithm.new("DsStarDs")
alg_ds.set_header(["DsStarDsAlg/DsStarDs.h"])
      .set_constant({"ECMS" => [:double, 4.178]})
      .note(:kinematic_fit, "final fit constrains m(Ds+) and the Ds+gamma (Ds*+) mass / recoil mass to their known values with a chi^2 cut; the corrected four-momenta feed the helicity (theta0, theta1, phi1, m12) unbinned ML analysis")
      .note(:mass_window, "K_S0 mass window 0.487-0.511 GeV applied on the pi+pi- secondary-vertex candidates")
      .note(:background_veto, "DeltaE and Ds-mass cuts (Table 1) with M(Ds+ gamma) vs recoil-mass M(Ds-) bands; overlapping modes rejected so the background stays below 8%")

sel_ds = Selection.new
sel_ds.select_track {
            cos_theta 0.93          # |cos(theta)| < 0.93
            Vz 10.0                 # |Vz| < 10 cm
            Vr 1.0                  # Vr < 1 cm
          }
       .select_photon {
            energyThreshold_b 0.025 # E > 25 MeV (barrel)
            energyThreshold_e 0.050 # E > 50 MeV (endcap)
            tdc_emc_start 0         # EMC time 0 - 700 ns
            tdc_emc_end 14
          }
       .pid(method: :probability) {
            prob_cut 0.001
            identify :kaon, against: [:pion, :proton]   # K+ / K-
            identify :pion, against: [:kaon, :proton]   # pi+ / pi-
          }
       .secondary_vertex_fit([:pip, :pim]) {            # K_S0 -> pi+ pi-
            build_virtual_particle(:K_S0).by_minimizing_mass_difference
            remove_used_particle_from_candidate_list
          }
       .partial_rec([1]) {                              # reconstruct the Ds*+ subtree (recID 1)
            best_combination_by_mass :"D_s+", 1.96835   # m(Ds+)
            require_recoil_mass 1.90, 2.04              # recoil-mass band around m(Ds-)
          }

alg_ds.with_decay_card(decay_card_ds).apply(sel_ds)

# ---------------------------------------------------------------------------
# Mode II : e+e- -> D*0 anti-D0 , reconstruct D*0 -> D0 pi0 -> K- pi+ pi0
# ---------------------------------------------------------------------------
alg_d0 = Algorithm.new("DStar0D0bar")
alg_d0.set_header(["DStar0D0barAlg/DStar0D0bar.h"])
      .set_constant({"ECMS" => [:double, 4.178]})
      .note(:kinematic_fit, "final fit constrains m(D0) and the D0 pi0 (D*0) mass / recoil mass to their known values with a chi^2 cut; corrected four-momenta feed the helicity unbinned ML analysis")
      .note(:mass_window, "pi0 mass window 0.115-0.150 GeV applied on the gamma gamma pair")
      .note(:background_veto, "DeltaE and D0-mass cuts (Table 1) with M(D0 pi0) vs recoil-mass M(anti-D0) bands; overlaps rejected, background kept below 8%")

sel_d0 = Selection.new
sel_d0.select_track {
            cos_theta 0.93
            Vz 10.0
            Vr 1.0
          }
      .select_photon {
            energyThreshold_b 0.025
            energyThreshold_e 0.050
            tdc_emc_start 0
            tdc_emc_end 14
          }
      .pid(method: :probability) {
            prob_cut 0.001
            identify :kaon, against: [:pion, :proton]
            identify :pion, against: [:kaon, :proton]
          }
      .kalman_kinematic_fit([:gamma, :gamma]) {          # pi0 -> gamma gamma
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
            chi2_cut 25
            npi0 ">=1"
          }
      .partial_rec([1]) {                                # reconstruct the D*0 subtree (recID 1)
            best_combination_by_mass :D0, 1.86483        # m(D0)
            require_recoil_mass 1.80, 1.93               # recoil-mass band around m(anti-D0)
          }

alg_d0.with_decay_card(decay_card_d0).apply(sel_d0)

# ---------------------------------------------------------------------------
# Mode III : e+e- -> D*+ D- , reconstruct D*+ -> D+ pi0 -> K- pi+ pi+ pi0
# ---------------------------------------------------------------------------
alg_dp = Algorithm.new("DStarD")
alg_dp.set_header(["DStarDAlg/DStarD.h"])
      .set_constant({"ECMS" => [:double, 4.178]})
      .note(:kinematic_fit, "final fit constrains m(D+) and the D+ pi0 (D*+) mass / recoil mass to their known values with a chi^2 cut; corrected four-momenta feed the helicity unbinned ML analysis")
      .note(:mass_window, "pi0 mass window 0.115-0.150 GeV applied on the gamma gamma pair")
      .note(:background_veto, "DeltaE and D+-mass cuts (Table 1) with M(D+ pi0) vs recoil-mass M(D-) bands; overlaps rejected, background kept below 8%")

sel_dp = Selection.new
sel_dp.select_track {
            cos_theta 0.93
            Vz 10.0
            Vr 1.0
          }
      .select_photon {
            energyThreshold_b 0.025
            energyThreshold_e 0.050
            tdc_emc_start 0
            tdc_emc_end 14
          }
      .pid(method: :probability) {
            prob_cut 0.001
            identify :kaon, against: [:pion, :proton]
            identify :pion, against: [:kaon, :proton]
          }
      .kalman_kinematic_fit([:gamma, :gamma]) {          # pi0 -> gamma gamma
            invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
            chi2_cut 25
            npi0 ">=1"
          }
      .partial_rec([1]) {                                # reconstruct the D*+ subtree (recID 1)
            best_combination_by_mass :"D+", 1.86966      # m(D+)
            require_recoil_mass 1.80, 1.94               # recoil-mass band around m(D-)
          }

alg_dp.with_decay_card(decay_card_dp).apply(sel_dp)

### ---------------------------- Execute on datasets ---------------------------- ###
root_files_ds = alg_ds.execute_on([data_4180, incMC_4180, exMC_ds])
root_files_d0 = alg_d0.execute_on([data_4180, incMC_4180, exMC_d0])
root_files_dp = alg_dp.execute_on([data_4180, incMC_4180, exMC_dp])