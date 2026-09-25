# Paper 2507.10331v2: Search for CLFV decay psi(3686) -> e+/- mu-/+
# BESIII, (2367.0±11.1)×10^6 psi(3686) events

psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

decay_card = <<~DECAYCARD
    Decay psi(2S)
    1.000 e+ mu- PHSP;
    Enddecay
    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "psip_emu_clfv_signal_mc"
  config.related_dataset = psip_data
  config.events = 100000
  config.decay_card = decay_card
  config.cross_section = :default
end

alg = Algorithm.new("PsiP2eMuCLFV")
alg.set_header(["PsiP2eMuCLFVAlg/PsiP2eMuCLFV.h"])
  .set_constant({ "ECMS" => [:double, 3.686] })

event_selection = Selection.new
  .select_track do
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nTot "==2"
    nNet "==0"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                   treat_as_electron_if_energy_above: 0.94
    nlp "==1"
    nlm "==1"
    nep "==1"
    nem "==1"
  end
  .kinematic_fit([:ep, :mum]) do
    nominal
    constrain_four_momentum
    chi2_cut 200
  end

alg.note(:muon_pid, "Muon EMC energy [0.1, 0.3] GeV, MUC depth > 20 cm, MUC layers > 3, MUC chi2 < 100, chi2_dE/dx_e < -1.1")
  .note(:electron_pid, "E/p > 0.94, |chi2_dE/dx_e| < 1.2, MUC hits == 0")
  .note(:TOF_cosmic_veto, "|TOF(e) - TOF(mu)| < 1.0 ns")
  .note(:photon_veto, "N_gamma == 0 (photon veto: E>25 MeV barrel, E>50 MeV endcap, TDC [0,700] ns, angle to track > 10 deg)")
  .note(:back_to_back, "|Delta_theta| < 1.4 deg, |Delta_phi| < 1.7 deg (back-to-back e-mu)")
  .note(:bhabha_veto, "Muon veto: 0.82<|cos(theta_mu)|<0.86 EMC gap; MUC-MDC polar angle diff < 20 deg, azimuthal diff < 30 deg")
  .note(:signal_region, "|sum(p_i)|/sqrt(s) <= 0.03, 0.97 <= sum(E_i)/sqrt(s) <= 1.04")
  .note(:charge_conjugate, "e- mu+ charge-conjugate mode handled analogously")
  .with_decay_card(decay_card)
  .apply(event_selection)

alg.execute_on([psip_data, psip_incMC, exMC_signal])