# Paper: 2305.07231v1
# Title: Search for baryon and lepton number violating decays of Xi0 hyperons
# Energy: 3.097 GeV (J/psi, single energy point)
# Final state: J/psi -> Xi0 anti-Xi0; anti-Xi0 -> anti-Lambda pi0 (tag);
#             Xi0 -> K± e∓ (signal, BNV)
# Double-tag technique: ST anti-Xi0 -> anti-Lambda pi0, DT Xi0 -> K± e∓
# BNV search with Delta(B-L)=0 and |Delta(B-L)|=2 channels

### Dataset preparation ###
data_708_3097 = DatasetManager.load_real_data.find("708_3097")
incMC_708_3097 = DatasetManager.load_inclusive_mc.find("708_3097")

all_data = [data_708_3097]
all_incMC = [incMC_708_3097]

# Decay card: J/psi -> Xi0 anti-Xi0
# Tag side: anti-Xi0 -> anti-Lambda pi0; anti-Lambda -> anti-p pi+
# Signal side: Xi0 -> K- e+ or Xi0 -> K+ e-
# pi0 -> gamma gamma
decay_card = <<~DECAYCARD
    Decay J/psi
    1.000  Xi0  anti-Xi0                         HELAMP 1.0 0.0;
    Enddecay

    Decay anti-Xi0
    1.000  anti-Lambda  pi0                       PHSP;
    Enddecay

    Decay Xi0
    1.000  K-  e+                                 PHSP;
    Enddecay

    Decay anti-Lambda
    1.000  anti-p-  pi+                           PHSP;
    Enddecay

    Decay pi0
    1.000  gamma  gamma                           PHSP;
    Enddecay
End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "Jpsi_Xi0_Xi0bar_BNV"
  config.events         = 500000
  config.decay_card     = decay_card
  config.cross_section  = :default
end

### Event selection ###
alg = Algorithm.new("Jpsi_Xi0Xi0bar_BNV")
alg.set_header(["Xi0BNVAlg/Xi0Xi0barBNV.h"])

event_selection = Selection.new

# ST selection: anti-Xi0 -> anti-Lambda pi0
# Charged tracks: anti-p, pi+ from anti-Lambda decay + K±, e∓ from signal
# Photon candidates for pi0 -> gamma gamma
event_selection.select_track {
                 cos_theta 0.93
                 Vz   10.0
                 Vr   1.0
               }
               # PID: proton and pion for anti-Lambda
               .pid(method: :probability) {
                 prob_cut 0.001
                 identify :proton, against: [:pion, :kaon]
                 identify :pion, against: [:kaon, :proton]
               }
               .select_photon {
                 energy 0.025         # E > 25 MeV in barrel
                 cos_theta_min 0.80   # barrel region
               }
               # Build anti-Lambda from anti-p pi+
               .build_virtual_particle(:anti_Lambda, [:prm, :pip]) {
                 secondary_vertex_fit
                 mass_window 0.005
               }
               # Build pi0 from gamma gamma
               .build_virtual_particle(:pi0, [:gamma, :gamma]) {
                 mass_window_lo 0.115
                 mass_window_hi 0.150
               }
               # Build anti-Xi0 from anti-Lambda pi0
               .build_virtual_particle(:anti_Xi0, [:anti_Lambda, :pi0]) {
                 mass_window 0.020
               }
               # 4C kinematic fit
               .kinematic_fit([:prm, :pip, :gamma, :gamma]) {
                 nominal
                 constrain_four_momentum
                 invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                 chi2_cut 200
               }

alg.with_decay_card(decay_card).apply(event_selection)

# ST selection: tagged anti-Xi0 reconstructed via anti-Lambda pi0
# anti-Lambda -> anti-p pi+ with vertex fit and mass window |M(pbar pi+) - M_Lambda| < 5 MeV
# pi0 -> gamma gamma, mass window (115, 150) MeV, mass-constrained kinematic fit
# ST yield from fit to M_BC distribution
alg.note(:st_selection,
  "ST: anti-Xi0 -> anti-Lambda pi0. anti-Lambda: vertex-fit pbar pi+, |M-M(Lambda)|<5 MeV. pi0: gamma gamma in (115,150) MeV/c^2, mass-constrained fit. ST M_BC fit with signal MC shape conv. Gaussian + 3rd-order Chebyshev. ST yield = 2,538,372+/-2593. ST efficiency = 16.26%. Applied in ROOT.")

# DT selection: Xi0 -> K± e∓ from remaining tracks recoiling against tagged anti-Xi0
# Electron PID: CL_e > 0.1% and CL_e/(CL_e+CL_pi+CL_K) > 0.8
# Kaon PID: kaon hypothesis has highest confidence level
# Signal region: |M(Ke)-M(Xi0)| < 20 MeV/c^2 and |M(Lambdabar pi0)-M(Xi0)| < 20 MeV/c^2
# Opening angle cut: theta_Xi0_Xibar0 > 178.5 deg
alg.note(:dt_selection,
  "DT: Xi0 -> K± e∓ recoiling vs anti-Xi0 tag. e PID: CL_e>0.1%, CL_e/(CL_e+CL_pi+CL_K)>0.8. K PID: highest CL. Opening angle > 178.5 deg. Signal region: |M(Ke)-M(Xi0)|<20 MeV, |M(Lambdabar pi0)-M(Xi0)|<20 MeV. 1 event in K-e+, 0 events in K+e-. Upper limit set via profile likelihood. Applied in ROOT.")

# BNV search: Delta(B-L) = 0 (K- e+) and |Delta(B-L)| = 2 (K+ e-)
# Upper limits: B(Xi0->K-e+) < 3.6x10^-6, B(Xi0->K+e-) < 1.9x10^-6 at 90% CL
alg.note(:bnv_upper_limits,
  "Upper limits at 90% CL: B(Xi0->K-e+) < 3.6e-6, B(Xi0->K+e-) < 1.9e-6. Frequentist profile likelihood method with systematic uncertainties. Applied in ROOT.")

# 1.0087 x 10^10 J/psi events
alg.note(:total_jpsi_events,
  "Total of (1.0087 +/- 0.0044) x 10^10 J/psi events used.")

all_datasets = all_data + all_incMC + exMC_signal
root_files = alg.execute_on(all_datasets)