# Search for K_S0 invisible decays via J/ψ → φ K_S0 K_S0
# J/ψ, 1.0087×10^10 events, arXiv:2501.06426v1

jpsi_data = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

decay_card = <<~DECAYCARD
  Decay J/psi
  1.000 phi K_S0 K_S0 PHSP;
  Enddecay

  Decay phi
  1.000 K+ K- VSS;
  Enddecay

  Decay K_S0
  1.000 pi+ pi- PHSP;
  Enddecay

  End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "jpsi_phi_KsKs_invisible"
  config.related_dataset = jpsi_data
  config.events = 100_000
  config.decay_card = decay_card
  config.cross_section = :default
end

# Algorithm: J/ψ → φ K_S0 K_S0, tag one K_S0 → π⁺π⁻ + φ → K⁺K⁻, signal K_S0 invisible
# φ reconstructed from K⁺K⁻, tag K_S0 via secondary vertex fit, invisible K_S0 via missing mass
alg = Algorithm.new("JpsiPhiKsKsInvisible")
alg.set_header(["JpsiPhiKsKsInvisibleAlg/JpsiPhiKsKsInvisible.h"])
   .set_constant({ "ECMS" => [:double, 3.097] })

# ── Selection ──
event_selection = Selection.new

event_selection
  .select_track {
    cos_theta 0.93
    Vz 100.0
    Vr 10.0
    nChrg "==4"
    nNet "==0"
  }
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    nGam ">=0"
  }
  # PID: kaons from φ (L(K) > L(π)), pions from K_S0 tag (L(π) > L(K))
  .pid(method: :probability) {
    identify :kaon, against: [:pion]
    identify :pion, against: [:kaon]
    nkp "==1"
    nkm "==1"
    npip "==1"
    npim "==1"
  }
  # Secondary vertex fit for K_S0 tag → π⁺π⁻
  .secondary_vertex_fit([:pip, :pim]) do
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  end
  # Kinematic fit with missing-particle treatment for invisible K_S0
  .kinematic_fit([:kp, :km, :K_S0]) do
    nominal
    miss_track_of :K_S0
    constrain_four_momentum
    chi2_cut 200
  end

alg.with_decay_card(decay_card).apply(event_selection)

# Additional cuts applied in ROOT:
# - φ mass window: [1.00, 1.04] GeV/c² for K⁺K⁻
# - K_S0 mass window: [0.486, 0.510] GeV/c² for π⁺π⁻
# - Recoil mass of φ > 1.08 GeV/c²
# - cosθ(φK_S0) in [-0.80, 0.80]
# - |M_recoil(φK_S0) - M(K_S0)| < 40 MeV/c²
# - K_S0 decay length significance > 2
# - Signal side: E_EMC (sum of EMC shower energies) for invisible search
alg.note(:root_cuts, "ROOT-level: φ mass [1.00,1.04] GeV/c², K_S0 mass [0.486,0.510] GeV/c², recoil mass of φ > 1.08 GeV/c², cosθ(φK_S0) in [-0.80,0.80], |M_recoil(φK_S0)-M(K_S0)|<40 MeV/c², K_S0 decay length significance > 2, E_EMC for invisible search")
alg.note(:signal_extraction, "Signal yield from fit to E_EMC distribution; upper limits set via Bayesian method")

alg.execute_on([jpsi_data, jpsi_incMC, exMC_signal])