### Dataset description ###
psip_data = DatasetManager.real_data.find("709_3686")        # ψ(2S) real data at 3.686 GeV
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")    # Corresponding inclusive MC sample

# Decay card for the signal process: ψ(2S) → γχ_c0 / γχ_c1 / γχ_c2 / γη_c(2S), each → p pbar
# (equal 25% rates among the four radiative channels, as specified for the signal MC)
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    0.25  gamma  chi_c0      PHSP;
    0.25  gamma  chi_c1      PHSP;
    0.25  gamma  chi_c2      PHSP;
    0.25  gamma  eta_c(2S)   PHSP;
    Enddecay

    Decay chi_c0
    1.000  p+  anti-p-       PHSP;
    Enddecay

    Decay chi_c1
    1.000  p+  anti-p-       PHSP;
    Enddecay

    Decay chi_c2
    1.000  p+  anti-p-       PHSP;
    Enddecay

    Decay eta_c(2S)
    1.000  p+  anti-p-       PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC sample (100k events) covering all four radiative channels
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3686_gamma_chic_ppbar"
  config.related_dataset = psip_data   # Associated real dataset
  config.events          = 100000      # Number of events to generate
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# All four signal channels share the same final state γ p pbar and identical selection,
# so a single Algorithm instance and a single decay card are used.
alg_name = "GammaChiCToPPbar"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})   # CMS energy 3.686 GeV

event_selection = Selection.new
event_selection.select_track {          # Charged track selection
                  cos_theta   0.93      # |cosθ| < 0.93
                  Vz          10.0      # |Vz| < 10 cm
                  Vr          1.0       # Vr < 1 cm in the transverse plane
                  nChrp       "==1"     # Exactly one positively charged track
                  nChrn       "==1"     # Exactly one negatively charged track
                  nNet        "==0"     # Net charge zero
                }
               .select_photon {         # Photon selection
                  tdc_emc_start     0       # TDC start time
                  tdc_emc_end       14      # TDC end time
                  energyThreshold_b 0.025   # E > 25 MeV in the EMC barrel
                  energyThreshold_e 0.050   # E > 50 MeV in the EMC endcap
                  nGam              ">=1"   # At least one photon
                }
               .pid(method: :probability) {   # Particle identification (probability method)
                  prob_cut   0.001               # PID probability > 0.001
                  identify :proton, against: [:kaon, :pion]  # p+ and anti-p- simultaneously
                  nprp   ">=1"                   # At least one proton
                  nprm   ">=1"                   # At least one anti-proton
                }
               .select_isolated_photon {         # Isolated photon against proton/anti-proton showers
                  angle_to_prp_track   20.0      # Opening angle to the positive (proton) track > 20°
                  angle_to_prm_track   30.0      # Opening angle to the negative (anti-proton) track > 30°
                  nGam   ">=1"                   # At least one photon survives isolation
                }
                # Nominal 4C kinematic fit to γ p pbar with a common vertex constraint on the two tracks
               .kinematic_fit([:gamma, :prp, :prm]) {
                  nominal                    # Nominal fit: corrected four-momenta are saved
                  vertex_fit([1, 2])         # Vertex fit on the prp (1) and prm (2) tracks
                  constrain_four_momentum    # 4C energy-momentum conservation
                  chi2_cut 60                # χ² < 60
                }
                # Competing 4C hypothesis p pbar (no nominal, no chi2 cut):
                # its χ² is stored so that the downstream ROOT analysis can require
                # χ²_4C(γ p pbar) < χ²_4C(p pbar), suppressing ψ(2S) → p pbar with a fake/FSR photon.
               .kinematic_fit([:prp, :prm]) {
                  constrain_four_momentum
                }

my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
# Same selection applied to real data, inclusive MC and the exclusive signal MC
root_files = my_algorithm.execute_on([psip_data, psip_incMC, exMC_signal])