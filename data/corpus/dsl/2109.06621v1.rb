# BESIII DSL: First observation of psi(3686) -> Xi(1530)0 anti-Xi(1530)0 and Xi(1530)0 anti-Xi0
# arXiv: 2109.06621v1
# Single-baryon tagging technique with psi(3686) data at 3.686 GeV

### Dataset preparation ###
psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Decay card for psi(3686) -> Xi(1530)0 anti-Xi(1530)0
decay_card_mode1 = <<~DECAYCARD
    Decay psi(3686)
    1.000 Xi(1530)0 anti-Xi(1530)0  PHSP;
    Enddecay

    Decay anti-Xi(1530)0
    1.000 anti-Xi- pi-              PHSP;
    Enddecay

    Decay Xi(1530)0
    1.000 Xi- pi+                   PHSP;
    Enddecay

    Decay anti-Xi-
    1.000 anti-Lambda pi-           PHSP;
    Enddecay

    Decay Xi-
    1.000 Lambda pi-                PHSP;
    Enddecay

    Decay Lambda
    1.000 p+ pi-                    HypWK;
    Enddecay

    Decay anti-Lambda
    1.000 anti-p- pi+               HypWK;
    Enddecay

    End
DECAYCARD

# Decay card for psi(3686) -> Xi(1530)0 anti-Xi0
decay_card_mode2 = <<~DECAYCARD
    Decay psi(3686)
    1.000 Xi(1530)0 anti-Xi0        PHSP;
    Enddecay

    Decay anti-Xi0
    1.000 anti-Lambda pi0           PHSP;
    Enddecay

    Decay Xi(1530)0
    1.000 Xi- pi+                   PHSP;
    Enddecay

    Decay Xi-
    1.000 Lambda pi-                PHSP;
    Enddecay

    Decay Lambda
    1.000 p+ pi-                    HypWK;
    Enddecay

    Decay anti-Lambda
    1.000 anti-p- pi+               HypWK;
    Enddecay

    Decay pi0
    1.000 gamma gamma               PHSP;
    Enddecay

    End
DECAYCARD

exMC_mode1 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "psi3686_Xi1530_Xi1530_exclusive_mc"
  config.related_dataset = psip_data
  config.events = 1_000_000
  config.decay_card = decay_card_mode1
  config.cross_section = :default
end

exMC_mode2 = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "psi3686_Xi1530_Xi0_exclusive_mc"
  config.related_dataset = psip_data
  config.events = 1_000_000
  config.decay_card = decay_card_mode2
  config.cross_section = :default
end

### Event selection for mode I: psi(3686) -> Xi(1530)0 anti-Xi(1530)0 ###
alg_mode1 = Algorithm.new("Psi3686Xi1530Xi1530")
alg_mode1.set_header(["Psi3686Xi1530Xi1530Alg/Psi3686Xi1530Xi1530.h"])
alg_mode1.set_constant({"ECMS" => [:double, 3.686]})

sel_mode1 = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 100.0
    Vr 10.0
    nChrp ">=3"
    nChrn ">=3"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp ">=1"
    nprm ">=1"
  end
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  .pid(method: :probability) do
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
  end
  # Reconstruct Lambda -> p pi- (tag side)
  .secondary_vertex_fit([:prp, :pim]) do
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  # Reconstruct anti-Lambda -> anti-p pi+ (tag side)
  .secondary_vertex_fit([:prm, :pip]) do
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end

sel_mode1.note(:cascade_reconstruction, "Xi- reconstructed as Lambda + pi- cascade secondary vertex; decay length > 0 required; mass windows: |M(p pi-)-M(Lambda)| < 5 MeV, |M(Lambda pi-)-M(Xi-)| < 8 MeV")
  .note(:xi1530_reconstruction, "Xi(1530)0 reconstructed from Xi- pi+ combination minimizing |M(Xi pi) - M(Xi(1530)0)|")
  # ...but we need to mark the kinematic fit. Let me just go to the final approach.

sel_mode1.note(:single_baryon_tag, "Single-baryon tagging: tag side reconstructs anti-Xi(1530)0 -> anti-Xi+ pi- (or Xi(1530)0 -> Xi- pi+); signal side inferred from recoil mass M_recoil(Xi- pi+)")
  .note(:recoil_mass_analysis, "Signal yields extracted from simultaneous unbinned ML fit to M_recoil(Xi- pi+) spectra; signal shapes from MC convolved with Gaussian; backgrounds: pi+pi-J/psi with J/psi->Xi-anti-Xi+, wrong-combination background, Chebychev polynomial for other backgrounds")
  .note(:angular_analysis, "Angular distribution parameter alpha for psi(3686)->Xi(1530)0 anti-Xi(1530)0 determined from efficiency-corrected cos(theta_B) fit; alpha = 0.32 +/- 0.19 +/- 0.07")
  .note(:branching_fractions, "B(psi(3686)->Xi(1530)0 anti-Xi(1530)0) = (6.77 +/- 0.14 +/- 0.39)e-5; B(psi(3686)->Xi(1530)0 anti-Xi0) = (0.53 +/- 0.04 +/- 0.03)e-5")
  .with_decay_card(decay_card_mode1)
  .apply(sel_mode1)

root_files_mode1 = alg_mode1.execute_on([psip_data, psip_incMC, exMC_mode1])

### Event selection for mode II: psi(3686) -> Xi(1530)0 anti-Xi0 ###
alg_mode2 = Algorithm.new("Psi3686Xi1530Xi0")
alg_mode2.set_header(["Psi3686Xi1530Xi0Alg/Psi3686Xi1530Xi0.h"])
alg_mode2.set_constant({"ECMS" => [:double, 3.686]})

sel_mode2 = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 100.0
    Vr 10.0
    nChrp ">=3"
    nChrn ">=3"
  end
  .select_photon do
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=2"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp ">=1"
    nprm ">=1"
  end
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  .pid(method: :probability) do
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
  end
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  end
  .secondary_vertex_fit([:prp, :pim]) do
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .secondary_vertex_fit([:prm, :pip]) do
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end

sel_mode2.note(:cascade_reconstruction, "Xi- reconstructed as Lambda + pi- cascade secondary vertex with positive decay length; mass windows: |M(p pi-)-M(Lambda)| < 5 MeV, |M(Lambda pi-)-M(Xi-)| < 8 MeV")
  .note(:single_baryon_tag, "Single-baryon tagging: tag side reconstructs anti-Xi(1530)0 -> anti-Xi+ pi-; signal side (anti-Xi0) inferred from recoil mass M_recoil(Xi- pi+)")
  .note(:pippim_jpsi_veto, "Veto psi(3686)->pi+pi-J/psi background: |M_recoil(pi+pi-) - M(J/psi)| > 4 MeV/c^2")
  .note(:recoil_mass_analysis, "Signal yields extracted from simultaneous unbinned ML fit to M_recoil spectra; signal shapes from MC convolved with Gaussian; WCB and pi+pi-J/psi background fixed")
  .with_decay_card(decay_card_mode2)
  .apply(sel_mode2)

root_files_mode2 = alg_mode2.execute_on([psip_data, psip_incMC, exMC_mode2])