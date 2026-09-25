### Dataset description ###
# Multi-energy scan: 19 c.m. energies from 4395.38 to 4950.93 MeV
# Using representative energy point 4681.92 MeV (highest significance 4.7sigma)
# For full scan, use DatasetManager.real_data.where(...) to collect all energy points

rscan_data_4682 = DatasetManager.real_data.find("703_4682")
rscan_incMC_4682 = DatasetManager.inclusive_mc.find("703_4682")

# Decay card for signal: e+e- -> K_S0 K- pi+ J/psi + c.c., J/psi -> l+l-
# Uses KKMC + psi(4260) as top mother (energy-scan convention)
decay_card_signal_ks_ee = <<~DECAYCARD
    Decay psi(4260)
    1.000 K_S0 K- pi+ J/psi PHSP;
    Enddecay

    Decay J/psi
    1.000 e+ e- VLL;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

decay_card_signal_ks_mumu = <<~DECAYCARD
    Decay psi(4260)
    1.000 K_S0 K- pi+ J/psi PHSP;
    Enddecay

    Decay J/psi
    1.000 mu+ mu- VLL;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

exMC_ks_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exclusive_K0SKpiJpsi_ee"
  config.related_dataset = rscan_data_4682
  config.events = 100000
  config.decay_card = decay_card_signal_ks_ee
  config.cross_section = :default
end

exMC_ks_mumu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exclusive_K0SKpiJpsi_mumu"
  config.related_dataset = rscan_data_4682
  config.events = 100000
  config.decay_card = decay_card_signal_ks_mumu
  config.cross_section = :default
end

### Event selection (BOSS) ###

# ============================================================
# Method I: Full reconstruction of e+e- -> K_S0 K- pi+ J/psi, J/psi -> e+e-
# Final state: pi+, pi- (from K_S0), K-, pi+, e+, e-
# ============================================================

alg_full_ee = Algorithm.new("K0SKpiJpsiFullEE")
alg_full_ee.set_header(["K0SKpiJpsiFullEEAlg/K0SKpiJpsiFullEE.h"])

selection_full_ee = Selection.new
  .select_track {
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp "==3"
    nChrn "==3"
  }
  .pid(method: :probability) {
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 1.0
    identify :kaon, against: [:pion, :proton]
    nkp "==1"
    nkm "==0"
    nlp "==1"
    nlm "==1"
  }
  .remove([:kp <= :chrgp, :lp <= :chrgp, :lm <= :chrgn])
  .assign({chrgp: :pip, chrgn: :pim})
  # Reconstruct K_S0 -> pi+ pi- via secondary vertex fit
  .secondary_vertex_fit([:pip, :pim]) {
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # Main 4C kinematic fit
  .kinematic_fit([:K_S0, :kp, :pip, :lp, :lm]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_full_ee
  .note(:helix_correction, "helix parameters of charged tracks corrected in MC to improve data-MC agreement of kinematic fit chi2 distributions")
  .note(:multi_energy, "analysis performed at 19 c.m. energies from 4395 to 4951 MeV; this DSL represents the full-reconstruction method (Method I) with KS0 and J/psi->e+e- at a single representative energy point")
  .note(:other_methods, "two additional reconstruction methods used: Method II (missing a K± or pi±, 1C fit) and Method III (missing a K0, 1C fit); both J/psi->e+e- and J/psi->mu+mu- final states analyzed; KL0 mode also reconstructed via missing-KL0 method")
  .with_decay_card(decay_card_signal_ks_ee)
  .apply(selection_full_ee)

# ============================================================
# Method I: Full reconstruction with J/psi -> mu+mu- (MUC depth cut)
# ============================================================

alg_full_mumu = Algorithm.new("K0SKpiJpsiFullMuMu")
alg_full_mumu.set_header(["K0SKpiJpsiFullMuMuAlg/K0SKpiJpsiFullMuMu.h"])

selection_full_mumu = Selection.new
  .select_track {
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp "==3"
    nChrn "==3"
  }
  .pid(method: :probability) {
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.4
    identify :kaon, against: [:pion, :proton]
    nkp "==1"
    nkm "==0"
    nlp "==1"
    nlm "==1"
  }
  .remove([:kp <= :chrgp, :lp <= :chrgp, :lm <= :chrgn])
  .assign({chrgp: :pip, chrgn: :pim})
  .secondary_vertex_fit([:pip, :pim]) {
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:K_S0, :kp, :pip, :lp, :lm]) {
    nominal
    constrain_four_momentum
    chi2_cut 200
  }

alg_full_mumu
  .note(:helix_correction, "helix parameters of charged tracks corrected in MC to improve data-MC agreement of kinematic fit chi2 distributions")
  .note(:muc_depth_cut, "MUC hit depth sum (Depth_mu+ + Depth_mu-) > 75 cm applied to suppress mu/pi misidentification background in the mu+mu- channel; applied in ROOT analysis stage")
  .with_decay_card(decay_card_signal_ks_mumu)
  .apply(selection_full_mumu)

# Execute the algorithms
root_files_full_ee = alg_full_ee.execute_on([rscan_data_4682, rscan_incMC_4682, exMC_ks_ee])
root_files_full_mumu = alg_full_mumu.execute_on([rscan_data_4682, rscan_incMC_4682, exMC_ks_mumu])