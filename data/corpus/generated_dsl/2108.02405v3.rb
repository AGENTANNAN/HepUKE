# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
data_3773 = DatasetManager.real_data.find("712_3773")        # ψ(3770) real data at √s = 3.773 GeV
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")    # Corresponding inclusive MC sample

# Decay card for the signal process ψ(3770) → D0 anti-D0, D0 → ω φ,
# with the full chain ω → π+π−π0, φ → K+K−, π0 → γγ.
# The anti-D0 side is left to decay generically (EvtGen default decay table).
decay_card_for_signal = <<~DECAYCARD
    Decay psi(3770)
    1.0000 D0 anti-D0 PHSP;
    Enddecay

    Decay D0
    1.0000 omega phi PHSP;
    Enddecay

    Decay omega
    1.0000 pi+ pi- pi0 OMEGA_DALITZ;
    Enddecay

    Decay phi
    1.0000 K+ K- VSS;
    Enddecay

    Decay pi0
    1.0000 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive MC for the signal D0 → ω φ (full ω/φ/π0 decay chain)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name    = "exmc_3773_D0_omega_phi"
  config.related_dataset = data_3773          # associated real dataset
  config.events         = 100000              # 100k events
  config.decay_card     = decay_card_for_signal
  config.cross_section  = :default
end

### Event selection (BOSS) ###
alg_name = "D0ToOmegaPhi"
my_Algorithm = Algorithm.new(alg_name)
my_Algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.773]})   # √s = 3.773 GeV

event_selection = Selection.new
event_selection.select_track {            # Charged track selection
                  cos_theta 0.93    # |cosθ| < 0.93
                  Vz        10.0    # |Vz| < 10 cm
                  Vr        1.0     # Vr < 1 cm
                  nChrp     ">=2"   # at least 2 positive tracks
                  nChrn     ">=2"   # at least 2 negative tracks
                }
               .select_photon {           # Photon selection
                  energyThreshold_b 0.025  # 25 MeV in the barrel
                  energyThreshold_e 0.050  # 50 MeV in the endcap
                  angle_to_track    10.0   # at least 10° from any charged track
                  nGam              ">=2"  # at least 2 photons (π0 → γγ)
                }
               .pid(method: :probability) {  # PID by the probability method, 0.001 cut
                  prob_cut 0.001
                  identify :pion, against: [:kaon]  # π+ / π− separated from kaons
                  identify :kaon, against: [:pion]  # K+ / K− separated from pions
                  npip ">=1"
                  npim ">=1"
                  nkp  ">=1"
                  nkm  ">=1"
                }
               .kalman_kinematic_fit([:gamma, :gamma]) {  # Reconstruct π0 → γγ
                  invariant_mass_of(:gamma, :gamma).within(0.115, 0.150)  # pre-fit M(γγ) window
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)  # constrain to nominal π0 mass
                  chi2_cut 25
                  npi0 ">=1"
                }
               # Main 4C kinematic fit to π+π−π0K+K−
               .kinematic_fit([:pip, :pim, :pi0, :kp, :km]) {
                  nominal                  # nominal fit — its corrected four-momenta are used
                  constrain_four_momentum  # constrain total four-momentum to the CMS energy
                  chi2_cut 200             # loose χ² < 200 (tight cut applied in ROOT)
                }
# NOTE: the downstream D0 candidate selection (M_BC > 1.84 GeV/c², −0.03 < ΔE < 0.02 GeV,
# best-candidate choice by minimum |ΔE|, K_S0 veto on M(π+π−) ∈ [0.490, 0.503] GeV/c²),
# the 2D unbinned ML fit to M(π+π−π0) vs M(K+K−), and the f_L extraction from cosθω / cosθK
# are post-kinematic-fit procedures and belong to the ROOT analysis stage, not to BOSS.

my_Algorithm.with_decay_card(decay_card_for_signal).apply(event_selection)
root_files = my_Algorithm.execute_on([data_3773, incMC_3773, exMC_signal])