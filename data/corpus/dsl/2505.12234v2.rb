# DSL auto-generated from 2505.12234v2
# Paper: Observation of χ_cJ → p pbar η η at BESIII
# Analysis: ψ(3686) → γ χ_cJ, χ_cJ → p pbar η η, η → γγ each
# Single energy: 709_3686 (ψ(2S))

### Dataset preparation ###
psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

decay_card = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_cJ PHSP;
    Enddecay

    Decay chi_cJ
    1.000 p+ anti-p- eta eta PHSP;
    Enddecay

    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_chicJ_ppbar_etaeta"
  config.related_dataset = psip_data
  config.events = 100_000
  config.decay_card = decay_card
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg = Algorithm.new("ChicJPPbarEtaEta")
alg.set_header(["ChicJPPbarEtaEtaAlg/ChicJPPbarEtaEta.h"])
    .set_constant({ "ECMS" => [:double, 3.686] })
    .note(:pi0_veto, "|M(γγ) - m_pi0| > 15 MeV pi0 veto; applied in ROOT via invariant mass window")
    .note(:eta_pair_selection, "ηη pair selected by minimizing Δ = sqrt((M(γ1γ2)-m_eta)^2 + (M(γ3γ4)-m_eta)^2) in ROOT")
    .note(:eta_signal_window, "|M(γγ) - m_eta| < 20 MeV applied in ROOT")
    .note(:chi2_cut, "Paper applies chi2_4C < 35; using loose chi2_cut 200 per Rule T3 (tight cut is ROOT-level)")
    .note(:chi2_comparison, "chi2(5γ hypothesis) < chi2(4γ hypothesis) and chi2(6γ hypothesis); competing hypothesis comparison in ROOT")

sel = Selection.new
sel.select_track do
      nChrp "==1"
      nChrn "==1"
      nNet "==0"
      cos_theta 0.93
      Vz 100.0
      Vr 10.0
    end
    .select_photon do
      tdc_emc_start 0
      tdc_emc_end 14          # EMC time [0, 700] ns (14 × 50 ns)
      angle_to_track 10.0
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      nGam ">=5"              # 4 for two η → γγ + 1 transition photon from ψ(2S)→γ χ_cJ
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
    # Reconstruct two η from photon pairs via Kalman fit
    .kalman_kinematic_fit([:gamma, :gamma]) do
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
      chi2_cut 25
      neta ">=2"
    end
    # Competing 5γ + p pbar hypothesis (no chi2_cut, no nominal — for ROOT-level veto)
    .kinematic_fit([:prp, :prm, :gamma, :gamma, :gamma, :gamma, :gamma]) do
      constrain_four_momentum
    end
    # 4C kinematic fit: e+e- → p pbar η η (2 composite eta tokens)
    .kinematic_fit([:prp, :prm, :eta, :eta]) do
      nominal
      constrain_four_momentum
      chi2_cut 200
    end

alg.with_decay_card(decay_card).apply(sel)
alg.execute_on([psip_data, psip_incMC, exMC_signal])