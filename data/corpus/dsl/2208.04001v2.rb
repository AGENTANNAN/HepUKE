# Paper: 2208.04001v2
# Title: Measurement of the Branching Fraction of the Singly Cabibbo-Suppressed
#        Decay Lambda_c+ -> Lambda K+
# Energy: 4.599-4.950 GeV (13 energy points, 6.44 fb^-1)
# Single-tag (ST) method: Lambda_c+ fully reconstructed from Lambda + bachelor K+/pi+
# Lambda -> p pi-
# Signal mode: Lambda_c+ -> Lambda K+
# Reference mode: Lambda_c+ -> Lambda pi+
# Relative BF measurement via simultaneous M_BC fit

### Dataset preparation ###
data_703_4600 = DatasetManager.load_real_data.find("703_4600")

all_data = [data_703_4600]
all_incMC = DatasetManager.load_inclusive_mc

# Decay card: e+e- -> Lambda_c+ anti-Lambda_c-
# Signal: Lambda_c+ -> Lambda K+, Lambda -> p pi-
# anti-Lambda_c- inclusive
decay_card_LcK = <<~DECAYCARD
    Decay e+ e-
    1.000  Lambda_c+  anti-Lambda_c-               VSP_PWAVE 1.0 0.0 0.0;
    Enddecay

    Decay Lambda_c+
    1.000  Lambda  K+                              PHSP;
    Enddecay

    Decay Lambda
    1.000  p+  pi-                                 PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.000  anti-p-  K+  pi-                        PHSP;
    Enddecay
End
DECAYCARD

# Reference: Lambda_c+ -> Lambda pi+
decay_card_LcPi = <<~DECAYCARD
    Decay e+ e-
    1.000  Lambda_c+  anti-Lambda_c-               VSP_PWAVE 1.0 0.0 0.0;
    Enddecay

    Decay Lambda_c+
    1.000  Lambda  pi+                             PHSP;
    Enddecay

    Decay Lambda
    1.000  p+  pi-                                 PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.000  anti-p-  K+  pi-                        PHSP;
    Enddecay
End
DECAYCARD

exMC_LcK = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "Lc_to_Lambda_K"
  config.events         = 500000
  config.decay_card     = decay_card_LcK
  config.cross_section  = :default
end

exMC_LcPi = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "Lc_to_Lambda_pi"
  config.events         = 500000
  config.decay_card     = decay_card_LcPi
  config.cross_section  = :default
end

### Event selection for Lambda_c+ -> Lambda K+ / Lambda pi+ ###
alg_LcK = Algorithm.new("Lc_to_Lambda_K")
alg_LcK.set_header(["LcLambdaKAlg/LcLambdaK.h"])

event_selection_LcK = Selection.new

# Charged track selection
event_selection_LcK.select_track {
                 cos_theta 0.93
                 Vz   10.0   # tight for prompt bachelor track
                 Vr   1.0
               }
               .pid(method: :probability) {
                 prob_cut 0.001
                 identify :proton, against: [:kaon, :pion]
                 identify :pion, against: [:kaon]
                 identify :kaon, against: [:pion]
               }
               .assign({:prp => :prp, :pim => :pim, :kp => :kp})

# Build Lambda from p pi- with secondary vertex fit
event_selection_LcK.build_virtual_particle(:Lambda, [:prp, :pim]) {
                 secondary_vertex_fit
                 mass_window_lo 1.111
                 mass_window_hi 1.121
               }

# Lambda_c+ identified by M_BC and DeltaE
# DeltaE: -0.009 < DeltaE < 0.012 GeV
# (applied in ROOT via algorithm.note)

alg_LcK.with_decay_card(decay_card_LcK).apply(event_selection_LcK)

alg_LcK.note(:lambdac_reconstruction,
  "Lambda_c+ reconstruction: Lambda -> p pi- (secondary vertex fit, decay length > 2*resolution). Bachelor K+/pi+ from IP (tight track). DeltaE = E(Lambda_c+) - E_beam, M_BC = sqrt(E_beam^2 - |p(Lambda_c+)|^2). DeltaE window: (-0.009, 0.012) GeV. Multiple candidates: keep minimum |DeltaE|. E/p < 0.9 for kaon candidate to suppress Lambda_c+ -> Lambda e+ nu_e background. Applied in ROOT.")

# Signal extraction: simultaneous M_BC fit
alg_LcK.note(:signal_extraction,
  "Simultaneous unbinned ML fit to M_BC distributions of Lambda_c+ -> Lambda K+ and Lambda_c+ -> Lambda pi+ across 13 c.m. energies. Signal: MC shape convolved with Gaussian. Background: ARGUS function. Relative BF R = B(Lambda_c+ -> Lambda K+) / B(Lambda_c+ -> Lambda pi+) = (4.78+/-0.34+/-0.20)%. Absolute BF using world-average B(Lambda_c+ -> Lambda pi+): B(Lambda_c+ -> Lambda K+) = (6.21+/-0.44+/-0.26+/-0.34)e-4. Applied in ROOT.")

# Energy points and luminosity
alg_LcK.note(:energy_points,
  "13 c.m. energies: 4.599, 4.612, 4.628, 4.641, 4.661, 4.682, 4.698, and 6 points above 4.700 GeV (merged). Total 6.44 fb^-1.")

### Reference mode ###
alg_LcPi = Algorithm.new("Lc_to_Lambda_pi")
alg_LcPi.set_header(["LcLambdaPiAlg/LcLambdaPi.h"])

event_selection_LcPi = Selection.new

event_selection_LcPi.select_track {
                  cos_theta 0.93
                  Vz   10.0
                  Vr   1.0
                }
                .pid(method: :probability) {
                  prob_cut 0.001
                  identify :proton, against: [:kaon, :pion]
                  identify :pion, against: [:kaon]
                }
                .assign({:prp => :prp, :pim => :pim, :pip => :pip})

event_selection_LcPi.build_virtual_particle(:Lambda, [:prp, :pim]) {
                  secondary_vertex_fit
                  mass_window_lo 1.111
                  mass_window_hi 1.121
                }

alg_LcPi.with_decay_card(decay_card_LcPi).apply(event_selection_LcPi)

alg_LcPi.note(:reference_mode,
  "Reference mode Lambda_c+ -> Lambda pi+ uses identical selection except bachelor pion PID: L(pi) > L(K). Track quality: |cos theta|<0.93, tight track for bachelor pion. Applied in ROOT.")

all_datasets = all_data + all_incMC + exMC_LcK + exMC_LcPi
root_files_LcK = alg_LcK.execute_on(all_datasets)
root_files_LcPi = alg_LcPi.execute_on(all_datasets)