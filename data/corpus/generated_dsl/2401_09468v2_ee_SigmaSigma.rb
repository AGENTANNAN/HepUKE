# =====================================================================================
# e+e- -> Sigma+ Sigma-   (Sigma+ -> p pi0, Sigma- -> anti-p pi0, pi0 -> gamma gamma)
# 41-point R scan, 3.510 - 4.951 GeV ; observed final state: p anti-p gamma gamma gamma gamma
# =====================================================================================

### Dataset preparation ###
# 41 R-scan energy points (BOSS 713); sample naming convention is [BOSS]_[ECMS(MeV)]
scan_energy_points_mev = (3510..4951).step(36).to_a   # 3510, 3546, ... , 4950 -> 41 points

# Real scan data and the matching inclusive MC, one sample per energy point
r_scan_data  = scan_energy_points_mev.map { |e| DatasetManager.real_data.find("713_#{e}") }
r_scan_incMC = scan_energy_points_mev.map { |e| DatasetManager.inclusive_mc.find("713_#{e}") }

# ConExc decay card: continuum (R-scan) signal with ISR modelling for the Born cross section.
# The literal token "ConExc" switches the DSL to the no-KKMC simulation template and injects
# "Particle vpho <ECMS> 0.0" per energy point -> NO explicit Particle vpho line for a scan.
decay_card_signal = <<~DECAYCARD
    Decay vpho
    1.000 ConExc 44;
    Enddecay

    Decay Sigma+
    1.000 p+ pi0 PHSP;
    Enddecay

    Decay anti-Sigma-
    1.000 anti-p- pi0 PHSP;
    Enddecay

    Decay pi0
    1.000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# One and the same 50k-event signal MC sample per scan point, from the single ConExc card
exMC_signal = DatasetManager.create_exclusive_mc_for(r_scan_data) do |config|
  config.sample_name   = "sigma_sigma_conexc"   # auto-suffixed per energy point
  config.events        = 50_000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "SigmaPairScan"
my_algorithm = Algorithm.new(alg_name)
my_algorithm
  .set_header(["#{alg_name}Alg/#{alg_name}.h"])
  .set_constant({"ECMS" => [:double, 4.230]})   # representative value of the scan range
  # The CM energy is not a fixed constant in this analysis: it is read per run from the
  # measured beam energy. Not expressible in the DSL selection grammar.
  .note(:measured_beam_energy,
        "CMS energy is taken per run from the measured beam energy (conditions DB) instead of
         a fixed ECMS constant; the ECMS constant declared above is only a representative
         scan value used for template bookkeeping")

event_selection = Selection.new
event_selection
  .select_track {                 # Charged track selection
    cos_theta 0.93                # |cos(theta)| < 0.93
    Vz        10.0                # |Vz| < 10 cm along the beam axis
    Vr        1.0                 # Vr < 1 cm in the transverse plane
    nChrp     ">=2"               # at least 2 positively charged tracks (p, and possibly p-bar)
    nChrn     ">=2"               # at least 2 negatively charged tracks
    nNet      "==0"               # net charge zero
  }
  .select_photon {                # Photon selection
    tdc_emc_start     0           # EMC timing window start
    tdc_emc_end       14          # EMC timing window end
    energyThreshold_b 0.025       # barrel energy threshold 25 MeV
    energyThreshold_e 0.050       # endcap energy threshold 50 MeV
    angle_to_track    10.0        # at least 10 degrees away from any charged track
    nGam              ">=4"       # at least 4 photons (4 gamma from the two pi0)
  }
  .pid(method: :probability) {    # Probability PID (no leptons expected)
    prob_cut 0.001                # PID probability > 0.001
    identify :proton, against: [:kaon, :pion]   # p+ and anti-p- vs K and pi
    nprp ">=1"                    # at least one proton
    nprm ">=1"                    # at least one anti-proton
  }
  # Additional proton momentum requirement: p > 0.5 GeV/c (drop softer protons)
  .for_each(:prp) {
    where { three_momentum_of(:prp) < 0.5 }   # active subset: protons below 0.5 GeV/c
    remove                                    # erase them from the proton candidate list
  }
  .for_each(:prm) {
    where { three_momentum_of(:prm) < 0.5 }   # active subset: anti-protons below 0.5 GeV/c
    remove                                    # erase them from the anti-proton candidate list
  }
  # Step 1: reconstruct pi0 from each gamma-gamma pair (1C mass-constrained Kalman fit)
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # m(gamma gamma) -> m(pi0)
    chi2_cut 25                 # chi2 < 25 for the pi0 mass constraint
    npi0 ">=2"                  # at least two pi0 candidates (one per Sigma)
  }
  # Step 2: 6C kinematic fit on p anti-p pi0 pi0 = 4C energy-momentum conservation
  # (constrain_four_momentum) plus the two pi0 mass constraints already imposed by the
  # Kalman step above.
  .kinematic_fit([:prp, :prm, :pi0, :pi0]) {
    nominal                     # nominal fit: corrected four-momenta are the ones saved
    constrain_four_momentum     # total four-momentum constrained to the CMS energy
    chi2_cut 200                # loose chi2 cut; optimal cut applied later in ROOT
  }

# Generate the complete algorithm for the process described by the decay card
my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)

# Execute on all scan points (real data + inclusive MC) and the signal MC samples
root_files = my_algorithm.execute_on(r_scan_data + r_scan_incMC + exMC_signal)