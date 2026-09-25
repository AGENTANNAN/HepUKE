# BESIII DSL: Search for Xi- -> Xi0 e- nu_e in J/psi decays
# arXiv: 2108.09948v2
# Uses single-tag technique with Xi+ tag

### Dataset preparation ###
jpsi_data = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Decay card for the signal process J/psi -> Xi- anti-Xi+, anti-Xi+ -> anti-Lambda pi+, Xi- -> Xi0 e- anti-nu_e
decay_card = <<~DECAYCARD
    Decay J/psi
    1.0000 Xi- anti-Xi+         PHSP;
    Enddecay

    Decay anti-Xi+
    1.0000 anti-Lambda pi+      PHSP;
    Enddecay

    Decay Xi-
    1.0000 Xi0 e- anti-nu_e     PHSP;
    Enddecay

    Decay Xi0
    1.0000 Lambda pi0           PHSP;
    Enddecay

    Decay Lambda
    1.0000 p+ pi-               HypWK;
    Enddecay

    Decay anti-Lambda
    1.0000 anti-p- pi+          HypWK;
    Enddecay

    Decay pi0
    1.0000 gamma gamma          PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "xine_semileptonic_signal_exclusive_mc"
  config.related_dataset = jpsi_data
  config.events = 1_000_000
  config.decay_card = decay_card
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg = Algorithm.new("XiSemileptonic")
alg.set_header(["XiSemileptonicAlg/XiSemileptonic.h"])
alg.set_constant({"ECMS" => [:double, 3.097]})

event_selection = Selection.new
  # Charged track selection
  .select_track do
    cos_theta 0.93
    Vz 100.0
    Vr 10.0
    nChrp ">=3"
    nChrn ">=2"
    nTot "==5"   # total of 5 charged tracks (4 from ST + DT sides)
  end
  # Photon selection
  .select_photon do
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=2"
  end
  # Particle identification: proton and pion hypotheses
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
  # Reconstruct pi0 from photon pairs via Kalman fit
  .kalman_kinematic_fit([:gamma, :gamma]) do
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  end
  # Reconstruct Lambda -> p pi- with secondary vertex fit (signal side)
  .secondary_vertex_fit([:prp, :pim]) do
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  # Reconstruct anti-Lambda -> anti-p pi+ with secondary vertex fit (tag side)
  .secondary_vertex_fit([:prm, :pip]) do
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  # Kinematic fit: Xi- -> Xi0 e- nu_e (missing nu_e)
  # Note: the anti-Xi+ is tagged, Xi0 = Lambda pi0, electron identified
  .kinematic_fit([:Lambda_bar, :pip, :Lambda, :pi0, :em]) do
    nominal
    constrain_four_momentum
    miss_track_of :nu_e    # missing neutrino
    chi2_cut 200
  end

# Note: electron identification and tag-side cuts (mass windows, recoil mass)
# are captured via notes since they involve procedures not yet in DSL
alg.note(:electron_pid, "Electron identified via PID confidence level (CL_e) highest among e/pi/K hypotheses; PID method uses dE/dx and TOF")
  .note(:tag_mass_window, "Anti-Xi+ tagged via anti-Lambda pi+; mass windows: |M(pbar pi+) - M(Lambda)| < 5 MeV/c^2, |M(Lambda_bar pi+) - M(Xi+)| < 5 MeV/c^2; recoil mass cut [1.290, 1.342] GeV/c^2")
  .note(:signal_mass_windows, "Xi0 mass window: |M(Lambda pi0) - M(Xi0)| < 14.5 MeV/c^2; Xi0 momentum (0.79, 0.84) GeV/c")
  .note(:q2_analysis, "Signal yield extracted from fit to q^2 distribution; q^2 = missing mass squared of e-nu system")
  .note(:st_yield_fit, "ST yield extracted from binned ML fit to M(Lambda_bar pi+) distribution; ST yield = 1,780,070 +/- 1366")
  .note(:upper_limit, "No signal observed; upper limit B(Xi- -> Xi0 e- nu_e) < 2.59e-4 at 90% CL using Bayesian method")
  .with_decay_card(decay_card)
  .apply(event_selection)

root_files = alg.execute_on([jpsi_data, jpsi_incMC, exMC_signal])