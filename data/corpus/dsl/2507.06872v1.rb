# Paper 2507.06872v1: Search for J/psi -> K+ K+ e- e- + c.c. (LNV)
# BESIII, (10087±44)×10^6 J/psi events

psip_data = DatasetManager.real_data.find("708_3097")
psip_incMC = DatasetManager.inclusive_mc.find("708_3097")

decay_card = <<~DECAYCARD
    Decay J/psi
    1.000 K+ K+ e- e- PHSP;
    Enddecay
    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "jpsi_kkee_lnv_signal_mc"
  config.related_dataset = psip_data
  config.events = 100000
  config.decay_card = decay_card
  config.cross_section = :default
end

alg = Algorithm.new("Jpsi2KKeeLNV")
alg.set_header(["Jpsi2KKeeLNVAlg/Jpsi2KKeeLNV.h"])
  .set_constant({ "ECMS" => [:double, 3.097] })

event_selection = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nTot "==4"
    nNet "==0"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.8
    identify :kaon, against: [:pion]
    nkp "==2"
    nkm "==0"
    nlp "==0"
    nlm "==2"
  end
  # Nominal 4C kinematic fit for signal: K+ K+ e- e-
  .kinematic_fit([:kp, :kp, :lm, :lm]) do
    nominal
    constrain_four_momentum
    chi2_cut 20
  end
  # Competing hypothesis: K+ K- pi+ pi-
  .kinematic_fit([:kp, :km, :lp, :lm]) do
    constrain_four_momentum
  end
  # Competing hypothesis: K+ K- K+ K-
  .kinematic_fit([:kp, :km, :kp, :km]) do
    constrain_four_momentum
  end
  # Competing hypothesis: pi+ pi- pi+ pi-
  .kinematic_fit([:lp, :lm, :lp, :lm]) do
    constrain_four_momentum
  end

alg.note(:electron_pid, "CL_e > 0.001, CL_e/(CL_e+CL_K+CL_pi) > 0.8, 0.8 < E/p < 1.2; additional E/p upper bound of 1.2 applied")
  .note(:kaon_pid, "CL_K > 0, CL_K > CL_pi, E/p < 0.8 for kaon anti-electron veto")
  .note(:charge_conjugate, "K- K- e+ e+ charge-conjugate mode handled analogously with counters swapped")
  .note(:gamma_conversion_veto, "K+K-e+e- and pi+pi-e+e- combinations rejected if e+e- from gamma conversions")
  .note(:signal_region, "M(K+K+e-e-) in [3.07, 3.12] GeV/c^2 (J/psi mass window, 3sigma)")
  .note(:competing_chi2_veto, "ROOT-level: chi2_4c_signal is smallest among the 4 kinematic hypotheses")
  .with_decay_card(decay_card)
  .apply(event_selection)

alg.execute_on([psip_data, psip_incMC, exMC_signal])