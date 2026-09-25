# DSL auto-generated from 2505.14988v1
# Paper: Test of local realism via entangled Lambda Lambda_bar at BESIII
# Analysis: J/ψ → η_c γ, η_c → Λ Λ_bar, Λ → p π-, Λ_bar → pbar π+
# Single energy: 708_3097 (J/ψ)

### Dataset preparation ###
jpsi_data = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

decay_card = <<~DECAYCARD
    Decay J/psi
    1.000 gamma eta_c HELAMP 1 0 1 0 1 0 0 0 1 0 1 0 0 0;
    Enddecay

    Decay eta_c
    1.000 Lambda0 anti-Lambda0 PHSP;
    Enddecay

    Decay Lambda0
    1.000 p+ pi- HypWK;
    Enddecay

    Decay anti-Lambda0
    1.000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_etac_LambdaLambdabar"
  config.related_dataset = jpsi_data
  config.events = 100_000
  config.decay_card = decay_card
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg = Algorithm.new("LambdaLambdabarEntanglement")
alg.set_header(["LambdaLambdabarEntanglementAlg/LambdaLambdabarEntanglement.h"])
    .set_constant({ "ECMS" => [:double, 3.097] })
    .note(:sigma0_veto, "|M(γ Λ) - M_Σ0| > 9 MeV veto; applied in ROOT via invariant mass check on fit-corrected momenta")
    .note(:space_like_separation, "1/S ≤ L1/L2 ≤ S with S = (1+β_Λ)/(1-β_Λ) space-like separation requirement; applied in ROOT via decay-length ratio")
    .note(:lambda_mass_window, "Λ mass window [1.008, 1.124] GeV applied in ROOT after secondary vertex fit")
    .note(:proton_pion_separation, "Proton/pion separated by momentum: momentum > 0.4 GeV/c → proton, else pion; applied in ROOT")
    .note(:chi2_cut, "Paper applies chi2_4C < 30; using loose chi2_cut 200 per Rule T3 (tight cut is ROOT-level)")
    .note(:photon_angle_track, "Photon-proton/antiproton angle > 20°, photon-pion angle > 30° to reject hadronic shower photons; inexpressible in DSL — use photon isolation note")

sel = Selection.new
sel.select_track do
      nChrp ">=2"
      nChrn ">=2"
      nNet "==0"
      cos_theta 0.93
      Vz 100.0
      Vr 10.0
    end
    .select_photon do
      tdc_emc_start 0
      tdc_emc_end 14
      angle_to_track 10.0
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      nGam ">=1"
    end
    .pid(method: :probability) do
      prob_cut 0.001
      identify :proton, against: [:kaon, :pion]
      nprp ">=1"
      nprm ">=1"
    end
    .remove([:prp <= :chrgp])
    .remove([:prm <= :chrgn])
    .assign({ chrgp: :pip, chrgn: :pim })
    # Reconstruct Λ → p π- via secondary vertex fit
    .secondary_vertex_fit([:prp, :pim]) do
      build_virtual_particle(:Lambda).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
    end
    # Reconstruct Λ_bar → pbar π+ via secondary vertex fit
    .secondary_vertex_fit([:prm, :pip]) do
      build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
      remove_used_particle_from_candidate_list
    end
    # 4C kinematic fit: e+e- → γ Λ Λ_bar
    .kinematic_fit([:gamma, :Lambda, :Lambda_bar]) do
      nominal
      constrain_four_momentum
      chi2_cut 200
    end

alg.with_decay_card(decay_card).apply(sel)
alg.execute_on([jpsi_data, jpsi_incMC, exMC_signal])