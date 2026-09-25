# Paper: 2209.08464v3
# Title: Partial wave analysis of the charmed baryon hadronic decay
#        Lambda_c+ -> Lambda pi+ pi0
# Energy: 4.600-4.699 GeV (7 energy points)
# Single-tag (ST) method: Lambda_c+ fully reconstructed
# Lambda -> p pi-, pi0 -> gamma gamma
# First PWA of Lambda_c+ -> Lambda pi+ pi0
# 4.4 fb^-1 total

### Dataset preparation ###
data_703_4600 = DatasetManager.load_real_data.find("703_4600")

all_data = [data_703_4600]
all_incMC = DatasetManager.load_inclusive_mc

# Decay card: e+e- -> Lambda_c+ anti-Lambda_c-
# Lambda_c+ -> Lambda pi+ pi0
# anti-Lambda_c- inclusive
decay_card = <<~DECAYCARD
    Decay e+ e-
    1.000  Lambda_c+  anti-Lambda_c-               VSP_PWAVE 1.0 0.0 0.0;
    Enddecay

    Decay Lambda_c+
    1.000  Lambda  pi+  pi0                       PHSP;
    Enddecay

    Decay Lambda
    1.000  p+  pi-                                 PHSP;
    Enddecay

    Decay pi0
    1.000  gamma  gamma                            PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.000  anti-p-  K+  pi-                        PHSP;
    Enddecay
End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "Lc_to_Lambda_pipi0"
  config.events         = 500000
  config.decay_card     = decay_card
  config.cross_section  = :default
end

### Event selection ###
alg = Algorithm.new("Lc_to_Lambda_pipi0")
alg.set_header(["LcLambdaPiPi0Alg/LcLambdaPiPi0.h"])

event_selection = Selection.new

# Charged track selection: p, pi-, pi+ (bachelor from Lambda_c+), and potentially pi- from Lambda
event_selection.select_track {
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

# Photon selection: at least 2 photons for pi0
event_selection.select_photon {
                 energy 0.025
                 emc_time [0, 700]
                 min_angle_to_charged 10   # degrees
               }

# Build pi0 from gamma gamma with 1C kinematic fit
event_selection.build_virtual_particle(:pi0, [:gamma, :gamma]) {
                 mass_window_lo 0.115
                 mass_window_hi 0.150
               }

# Build Lambda from p pi- with secondary vertex fit
event_selection.build_virtual_particle(:Lambda, [:prp, :pim]) {
                 secondary_vertex_fit
                 mass_window_lo 1.111
                 mass_window_hi 1.121
               }

alg.with_decay_card(decay_card).apply(event_selection)

# Lambda and pi0 reconstruction details
alg.note(:lambda_reconstruction,
  "Lambda: p pi- secondary vertex fit (chi2<100). Decay length > 2*sigma. Proton PID: L(p)>L(K) and L(p)>L(pi). Pion from Lambda: no PID. pi0: gamma gamma with 1C kinematic fit constraining to pi0 mass. Sigma0 veto: M(Lambda gamma) in (1.179,1.203) GeV/c2 rejected.")

# Lambda_c+ reconstruction via M_BC and DeltaE
alg.note(:lambdac_reconstruction,
  "Lambda_c+ -> Lambda pi+ pi0. Single-tag reconstruction with M_BC and DeltaE. DeltaE window: (-0.03, 0.02) GeV. Best candidate: minimum |DeltaE|. M_BC signal region per energy point (Table 3). Signal yields from unbinned ML fits to M_BC (signal: MC shape convolved with Gaussian; background: ARGUS). Sum of signal yields ~10k events. 3C kinematic fit applied for PWA: constrain Lambda mass, Lambda_c+ mass, and recoil mass to Lambda_c+ mass. Applied in ROOT.")

# PWA results
alg.note(:pwa_results,
  "PWA with helicity amplitude formalism (TFPWA framework). Nominal components: Lambda rho(770)+, Sigma(1385)+ pi0, Sigma(1385)0 pi+, Sigma(1670)+ pi0, Sigma(1670)0 pi+, Sigma(1750)+ pi0, Sigma(1750)0 pi+, Lambda NR(1-). First measurement of B(Lambda_c+ -> Lambda rho(770)+) = (4.06+/-0.30+/-0.35+/-0.23)e-2, B(Lambda_c+ -> Sigma(1385)+ pi0) = (5.86+/-0.49+/-0.52+/-0.35)e-3, B(Lambda_c+ -> Sigma(1385)0 pi+) = (6.47+/-0.59+/-0.66+/-0.38)e-3. Decay asymmetry parameters also measured. Applied in ROOT.")

# 7 energy points
alg.note(:energy_points,
  "7 c.m. energies: 4.600, 4.612, 4.628, 4.641, 4.661, 4.682, 4.699 GeV. Total 4.4 fb^-1.")

all_datasets = all_data + all_incMC + exMC_signal
root_files = alg.execute_on(all_datasets)