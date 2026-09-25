# Paper: 2306.10319v3
# Title: Precise measurement of J/psi -> Lambda_bar pi+ Sigma- + c.c. and
#        J/psi -> Lambda_bar pi- Sigma+ + c.c.
# Energy: 3.097 GeV (J/psi, single energy point)
# Final state: J/psi -> Lambda_bar pi+ Sigma- (and charge conjugates)
# Lambda -> p pi-; Sigma not reconstructed (inferred from recoil mass)
# Four separate but related channels

### Dataset preparation ###
data_708_3097 = DatasetManager.load_real_data.find("708_3097")
incMC_708_3097 = DatasetManager.load_inclusive_mc.find("708_3097")

all_data = [data_708_3097]
all_incMC = [incMC_708_3097]

# Decay card: J/psi -> Lambda_bar pi+ Sigma-
# Lambda -> p pi-; Sigma- -> n pi- (not reconstructed, only recoil mass)
# Only reconstruct Lambda and prompt pion; Sigma inferred via M_recoil
decay_card = <<~DECAYCARD
    Decay J/psi
    1.000  anti-Lambda  pi+  Sigma-              HELAMP 1.0 0.0;
    Enddecay

    Decay anti-Lambda
    1.000  anti-p-  pi+                           PHSP;
    Enddecay

    Decay Sigma-
    1.000  n0  pi-                                PHSP;
    Enddecay
End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "Jpsi_Lambdabar_pi_Sigma"
  config.events         = 500000
  config.decay_card     = decay_card
  config.cross_section  = :default
end

### Event selection ###
alg = Algorithm.new("Jpsi_LambdabarPiSigma")
alg.set_header(["LambdabarPiSigmaAlg/LambdabarPiSigma.h"])

event_selection = Selection.new

# At least 3 charged tracks: proton, pion from Lambda + prompt pion from J/psi
# PID: identify proton and pion
event_selection.select_track {
                 cos_theta 0.93
                 nTot ">=3"
               }
               .pid(method: :probability) {
                 prob_cut 0.001
                 identify :proton, against: [:pion, :kaon]
                 identify :pion, against: [:kaon, :proton]
               }
               .assign({:prp => :prp, :pim => :pim})

# Build Lambda from p pi- with vertex fit
event_selection.build_virtual_particle(:Lambda, [:prp, :pim]) {
                 secondary_vertex_fit
                 mass_window 0.005   # |M(p pi-) - M(Lambda)| < 5 MeV
               }

# Prompt pion from J/psi: distance of closest approach cuts
# Sigma not reconstructed; inferred from M_recoil(Lambda_bar pi)
# M_recoil in [1.14, 1.24] GeV/c^2
alg.with_decay_card(decay_card).apply(event_selection)

# Lambda reconstruction:
# - Secondary vertex fit on p pi-
# - Mass window |M(p pi-) - M(Lambda)| < 5 MeV
# - Decay length > 2*sigma
# - Best Lambda candidate: min [M(p pi-) - M(Lambda)]^2
alg.note(:lambda_reconstruction,
  "Lambda: vertex-fit p pi-, |M-M(Lambda)| < 5 MeV. Decay length/error > 2. Best candidate selected by min [M-M(Lambda)]^2. Lambda reconstruction efficiency corrected via J/psi -> p K- Lambda_bar control sample. Applied in ROOT.")

# Prompt pion selection: distance of closest approach < 10 cm (Vz), < 1 cm (Vxy)
# PID: L(pi) > 0.001, L(pi) > L(p) and L(pi) > L(K)
# Sigma inferred from M_recoil(Lambda_bar pi) in [1.14, 1.24] GeV/c^2
alg.note(:sigma_reconstruction,
  "Sigma not reconstructed. Inferred from recoil mass of Lambda_bar pi system. M_recoil window [1.14,1.24] GeV/c^2. Binned ML fit to M_recoil distribution: signal=RooKeysPdf conv. Gaussian, background=3rd-order Chebyshev. Four separate channels measured. Applied in ROOT.")

# PWA: full reconstruction including Sigma -> n pi for efficiency determination
# 1C kinematic fit with missing neutron
# Detection efficiency from PWA MC weighted by data/MC efficiency corrections
alg.note(:pwa_efficiency,
  "PWA for efficiency: full reconstruction including Sigma -> n pi. 1C kinematic fit missing neutron (chi2_1C < 30). Lambda vertex chi2 < 30. 11 intermediate resonances in PWA. Detection efficiency from PWA MC weighted by Lambda reconstruction + pi tracking/PID data/MC corrections. Efficiency ~31-35%. Applied in ROOT.")

# Branching fractions:
# B(J/psi -> Lambda_bar pi+ Sigma- + c.c.) = (1.221+/-0.002+/-0.038)e-3
# B(J/psi -> Lambda_bar pi- Sigma+ + c.c.) = (1.244+/-0.002+/-0.045)e-3
# Isospin symmetry test and 12% rule comparison
alg.note(:branching_fractions,
  "Branching fractions: B(J/psi->Lambda_bar pi+ Sigma- +c.c.)=(1.221+/-0.002+/-0.038)e-3, B(J/psi->Lambda_bar pi- Sigma+ +c.c.)=(1.244+/-0.002+/-0.045)e-3. Isospin ratio = 0.98+/-0.04. 12% rule ratio with psi(3686) = (11.5+/-1.4)% and (12.4+/-1.4)%. QED background from 3.080 GeV data. Applied in ROOT.")

# (1.0087 +/- 0.0044) x 10^10 J/psi events
alg.note(:total_jpsi_events,
  "Total of (1.0087 +/- 0.0044) x 10^10 J/psi events used. Continuum data at 3.080 GeV (166 pb^-1) for QED background estimation.")

all_datasets = all_data + all_incMC + exMC_signal
root_files = alg.execute_on(all_datasets)