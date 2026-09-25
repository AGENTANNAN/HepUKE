# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")      # ψ(3686) data at √s = 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")   # corresponding inclusive MC sample

# Decay card — Mode A: ψ(2S) → γ χc1, χc1 → K+ K− K+ K− ω, ω → π+π−π0, π0 → γγ
decay_card_modeA = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c1 P2GC1;
    Enddecay

    Decay chi_c1
    1.000 K+ K+ K- K- omega PHSP;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Decay card — Mode B: ψ(2S) → γ χc1, χc1 → φ K+ K− ω, φ → K+ K−, ω → π+π−π0, π0 → γγ
# (identical final state to Mode A, but through the φ intermediate resonance)
decay_card_modeB = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c1 P2GC1;
    Enddecay

    Decay chi_c1
    1.000 phi K+ K- omega PHSP;
    Enddecay

    Decay phi
    1.000 K+ K- VSS;
    Enddecay

    Decay omega
    1.000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC samples (500k events each), both associated with the ψ(3686) data
exMC_modeA = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gammachic1_2K2Komega"
  config.related_dataset = psip_data
  config.events          = 500_000
  config.decay_card      = decay_card_modeA
  config.cross_section   = :default
end

exMC_modeB = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_gammachic1_phiKKomega"
  config.related_dataset = psip_data
  config.events          = 500_000
  config.decay_card      = decay_card_modeB
  config.cross_section   = :default
end

### Event selection (BOSS) — shared by both decay modes ###
alg_name_A = "ChicJTo2K2KOmega"
alg_A = Algorithm.new(alg_name_A)
alg_A.set_header(["#{alg_name_A}Alg/#{alg_name_A}.h"])
     .set_constant({"ECMS" => [:double, 3.686]})    # ECMS = 3.686 GeV

alg_name_B = "ChicJToPhiKKOmega"
alg_B = Algorithm.new(alg_name_B)
alg_B.set_header(["#{alg_name_B}Alg/#{alg_name_B}.h"])
     .set_constant({"ECMS" => [:double, 3.686]})

# Common selection chain for the final state γ K+ K+ K− K− π+ π− π0
event_selection = Selection.new
    .select_track {                 # charged-track selection
        cos_theta 0.93              # |cos θ| < 0.93
        Vz        10.0             # |Vz| < 10 cm
        Vr        1.0              # Vr < 1 cm
        nChrp     "==3"            # exactly three positive tracks
        nChrn     "==3"            # exactly three negative tracks
        nNet      "==0"            # net charge zero
    }
    .select_photon {                # photon selection
        tdc_emc_start     0        # EMC timing 0–700 ns
        tdc_emc_end       14
        angle_to_track    10.0     # at least 10° from any charged track
        energyThreshold_b 0.025    # > 25 MeV in the barrel
        energyThreshold_e 0.050    # > 50 MeV in the endcap
        nGam              ">=3"    # at least three photons
    }
    .pid(method: :probability) {    # kaon / pion identification (probability method)
        prob_cut 0.001             # PID probability > 0.001
        identify :kaon, against: [:pion, :proton]   # K+ and K−
        identify :pion, against: [:kaon, :proton]   # π+ and π−
        nkp  ">=2"                 # at least two K+
        nkm  ">=2"                 # at least two K−
        npip ">=1"                 # at least one π+
        npim ">=1"                 # at least one π−
    }
    .kalman_kinematic_fit([:gamma, :gamma]) {   # 1C mass-constrained fit for π0
        invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
        chi2_cut 25               # χ² < 25
        npi0 ">=1"                # at least one π0
    }
    # 5C kinematic fit (4C four-momentum conservation + the 1C π0 mass from the Kalman fit)
    .kinematic_fit([:gamma, :kp, :kp, :km, :km, :pip, :pim, :pi0]) {
        nominal                   # nominal fit — corrected four-momenta are saved
        constrain_four_momentum
        chi2_cut 40               # χ² < 40
    }

# Attach the shared selection to each chain's algorithm
alg_A.with_decay_card(decay_card_modeA).apply(event_selection.dup)
alg_B.with_decay_card(decay_card_modeB).apply(event_selection.dup)

# Execute on data, inclusive MC and the corresponding signal exclusive MC
root_files_A = alg_A.execute_on([psip_data, psip_incMC, exMC_modeA])
root_files_B = alg_B.execute_on([psip_data, psip_incMC, exMC_modeB])