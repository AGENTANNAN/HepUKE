# Paper: 2305.14631v1
# Title: Determination of spin and parity of D*(s) mesons
# Energy: 4.178 GeV (single energy point)
# Partial reconstruction: one D*(s) reconstructed, other D(s) undetected
# Processes: e+e- -> Ds*+ Ds-, D*0 D0bar, D*+ D-
# Ds*+ -> gamma Ds+, D*0 -> D0 pi0, D*+ -> D+ pi0
# Ds+ -> K_S0 K+, D0 -> K- pi+, D+ -> K- pi+ pi+
# Helicity amplitude analysis to determine J^P

### Dataset preparation ###
data_703_4178 = DatasetManager.load_real_data.find("703_4178")
incMC_703_4178 = DatasetManager.load_inclusive_mc.find("703_4178")

all_data = [data_703_4178]
all_incMC = [incMC_703_4178]

# Decay card: e+e- -> Ds*+ Ds- (partial reconstruction)
# Ds*+ -> gamma Ds+; Ds+ -> K_S0 K+
# Ds- undetected (inclusive)
decay_card_DsStar = <<~DECAYCARD
    Decay e+ e-
    1.000  D_s*+  D_s-                            VSP_PWAVE 1.0 0.0 0.0;
    Enddecay

    Decay D_s*+
    1.000  gamma  D_s+                             VSP_PWAVE 1.0 0.0 0.0;
    Enddecay

    Decay D_s+
    1.000  K_S0  K+                                PHSP;
    Enddecay

    Decay D_s-
    1.000  K+  K-  pi-                             PHSP;
    Enddecay

    Decay K_S0
    1.000  pi+  pi-                                PHSP;
    Enddecay
End
DECAYCARD

# Decay card: e+e- -> D*0 D0bar
# D*0 -> D0 pi0; D0 -> K- pi+
decay_card_D0Star = <<~DECAYCARD
    Decay e+ e-
    1.000  D*0  anti-D0                            VSP_PWAVE 1.0 0.0 0.0;
    Enddecay

    Decay D*0
    1.000  D0  pi0                                 VSP_PWAVE 1.0 0.0 0.0;
    Enddecay

    Decay D0
    1.000  K-  pi+                                 PHSP;
    Enddecay

    Decay anti-D0
    1.000  K+  pi-                                 PHSP;
    Enddecay

    Decay pi0
    1.000  gamma  gamma                            PHSP;
    Enddecay
End
DECAYCARD

# Decay card: e+e- -> D*+ D-
# D*+ -> D+ pi0; D+ -> K- pi+ pi+
decay_card_DStar = <<~DECAYCARD
    Decay e+ e-
    1.000  D*+  D-                                 VSP_PWAVE 1.0 0.0 0.0;
    Enddecay

    Decay D*+
    1.000  D+  pi0                                 VSP_PWAVE 1.0 0.0 0.0;
    Enddecay

    Decay D+
    1.000  K-  pi+  pi+                            PHSP;
    Enddecay

    Decay D-
    1.000  K+  pi-  pi-                            PHSP;
    Enddecay

    Decay pi0
    1.000  gamma  gamma                            PHSP;
    Enddecay
End
DECAYCARD

exMC_DsStar = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "ee_to_DsStar_Ds"
  config.events         = 500000
  config.decay_card     = decay_card_DsStar
  config.cross_section  = :default
end

exMC_D0Star = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "ee_to_D0Star_D0bar"
  config.events         = 500000
  config.decay_card     = decay_card_D0Star
  config.cross_section  = :default
end

exMC_DStar = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "ee_to_DStar_D"
  config.events         = 500000
  config.decay_card     = decay_card_DStar
  config.cross_section  = :default
end

### Common event selection for all three processes ###
alg = Algorithm.new("DsStar_SpinParity")
alg.set_header(["DsStarSpinParityAlg/DsStarSpinParity.h"])

event_selection = Selection.new

# Charged tracks: track quality + PID for kaons and pions
event_selection.select_track {
                 cos_theta 0.93
                 Vz   10.0
                 Vr   1.0
               }
               # PID: kaon vs pion
               .pid(method: :probability) {
                 prob_cut 0.001
                 identify :kaon, against: [:pion]
                 identify :pion, against: [:kaon]
               }
               .assign({:kp => :kp, :km => :km, :pip => :pip, :pim => :pim})

# Photon selection for pi0 and radiative photon
event_selection.select_photon {
                 energy 0.025
                 emc_time [0, 700]
               }

# Build K_S0 from pi+ pi-
event_selection.build_virtual_particle(:K_S0, [:pip, :pim]) {
                 secondary_vertex_fit
                 mass_window_lo 0.487
                 mass_window_hi 0.511
               }

# Build pi0 from gamma gamma
event_selection.build_virtual_particle(:pi0, [:gamma, :gamma]) {
                 mass_window_lo 0.115
                 mass_window_hi 0.150
               }

alg.with_decay_card(decay_card_DsStar).apply(event_selection)

# Partial reconstruction technique:
# - Reconstruct only one D*(s) per event
# - Ds+: K_S0 K+, D0: K- pi+, D+: K- pi+ pi+
# - Ds* -> gamma Ds, D* -> pi0 D
# - Use DeltaE and M(D(s)) cuts, plus M(D(s) pi0/gamma) and RM(D(s))
# - Tag/recoil regions defined in M(D(s)pi0/gamma) vs RM(D(s)) 2D plane
alg.note(:partial_reconstruction,
  "Partial reconstruction: one D*(s) reconstructed per event. Ds+: K_S0 K+, D0: K- pi+, D+: K- pi+ pi+. Ds* -> gamma Ds, D* -> pi0 D. DeltaE and M(D(s)) cuts per Table 1. Signal region: M(D(s)pi0/gamma) and RM(D(s)) bands. Horizontal/vertical band regions used; overlapping region rejected. Background <8%. Tag and recoil samples used for helicity amplitude analysis. Applied in ROOT.")

# Kinematic fit: constrain D(s) mass + D(s)* recoil mass
# Updated four-momenta used for helicity analysis
alg.note(:kinematic_fit,
  "Kinematic fit: D(s) mass constrained to known mass, recoil mass of D(s)pi0(gamma) constrained to known D(s) mass. Chi2 cut applied. Updated four-momenta used for helicity amplitude analysis. Applied in ROOT.")

# Helicity amplitude analysis:
# - Unbinned ML fit to angular distributions (theta0, theta1, phi1, m12)
# - Test J^P = 1-, 2+, 3- hypotheses
# - Background subtracted using inclusive MC
# - Significance >10sigma for J^P=1- over 2+ and 3-
alg.note(:helicity_analysis,
  "Helicity amplitude analysis: unbinned ML fit to joint angular distribution (theta0, theta1, phi1, m12). D*(s) spin-parity tested: 1-, 2+, 3-. Background subtracted with inclusive MC. Normalization from PHSP MC with efficiency. Results: J^P=1- confirmed with >10sigma significance for Ds*+, D*0, D*+. Moments <sin^2 theta1> vs phi1 also shown. Applied in ROOT.")

# 3.19 fb^-1 at 4.178 GeV
alg.note(:luminosity,
  "3.19 fb^-1 at 4.178 GeV.")

all_datasets = all_data + all_incMC + exMC_DsStar + exMC_D0Star + exMC_DStar
root_files = alg.execute_on(all_datasets)