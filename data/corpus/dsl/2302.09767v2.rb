# 2302.09767v2: Xi0 asymmetry parameters from psi(3686) -> Xi0 Xi0_bar
# psi(3686) -> Xi0 Xi0_bar -> pi0 pi0 Lambda Lambda_bar -> 4gamma + 2p2pi- + 2pbar2pi+
# 6C kinematic fit (energy-momentum conservation + 2 pi0 mass constraints)

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

psip_data  = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

decay_card_signal = <<~DECAYCARD
  Decay psi(3686)
  1.0000 Xi0 anti-Xi0  PHSP;
  Enddecay
  Decay Xi0
  1.0000 pi0 Lambda0  PHSP;
  Enddecay
  Decay anti-Xi0
  1.0000 pi0 anti-Lambda0  PHSP;
  Enddecay
  Decay Lambda0
  1.0000 p+ pi-  PHSP;
  Enddecay
  Decay anti-Lambda0
  1.0000 anti-p- pi+  PHSP;
  Enddecay
  Decay pi0
  1.0000 gamma gamma  PHSP;
  Enddecay
  End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "sig_Xi0_psip"
  config.related_dataset = psip_data
  config.events          = 5_000_000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

algorithm = Algorithm.new("Xi0AsymmetryPsip")
algorithm
  .set_header(["Xi0AsymmetryPsipAlg/Xi0AsymmetryPsip.h"])
  .set_constant({ "ECMS" => [:double, 3.686] })
  .with_decay_card(decay_card_signal)

event_selection = Selection.new
  .select_track do
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nChrp     ">=2"
    nChrn     ">=2"
  end
  .pid(method: :probability) do
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
  end
  .for_each(:prp) { where { p < 0.5 }; remove }
  .for_each(:prm) { where { p < 0.5 }; remove }
  .pid(method: :probability) do
    prob_cut 0.001
    identify :pion, against: [:kaon, :proton]
  end
  .for_each(:pip) { where { p > 0.5 }; remove }
  .for_each(:pim) { where { p > 0.5 }; remove }
  .select_photon do
    tdc_emc_start     0
    tdc_emc_end       14
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track    10.0
    nGam              ">=4"
  end
  .secondary_vertex_fit(:Lambda, daughters: [:prp, :pim]) do
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .secondary_vertex_fit(:Lambda_bar, daughters: [:prm, :pip]) do
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :Lambda, :Lambda_bar]) do
    constrain_four_momentum
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 200
    nominal
  end

algorithm
  .note(:background_veto, "psi(3686) -> pi0 pi0 J/psi veto: reject events if recoil mass of pi0 pi0 system is within 20 MeV/c^2 of nominal J/psi mass")
  .note(:xi0_selection, "Best Xi0/Xi0bar selected by minimizing sqrt((M_pi0Lambda - m_Xi0)^2 + (M_pi0Lambda_bar - m_Xi0)^2); Xi0 mass window |M - m_Xi0| < 15 MeV/c^2 applied in ROOT")
  .note(:lambda_decay_length, "Lambda decay length > 0 required to suppress non-Lambda background")
  .note(:helix_correction, "Helix parameter corrections applied to charged tracks for 6C kinematic fit momentum resolution differences between data and MC")
  .note(:tracking_efficiency, "Tracking and pi0 reconstruction efficiency corrections evaluated from control sample of psi(3686) -> Xi0 Xi0_bar single-tag events")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

algorithm.execute_on([psip_data, psip_incMC, exMC_signal])