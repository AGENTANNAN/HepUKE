# Paper: 2207.14350v4
# Title: Evidence of the isospin-violating decay psi(3686) -> Lambda Lambdabar pi0
#        and study of psi(3686) -> Lambda Lambdabar eta
# Energy: 3.686 GeV (psi(3686), single energy point)
# Final state: psi(3686) -> Lambda Lambdabar pi0/eta
# Lambda -> p pi-, Lambdabar -> pbar pi+, pi0/eta -> gamma gamma
# Common final state: p pbar pi+ pi- gamma gamma
# 4C kinematic fit (not constraining gamma gamma mass to distinguish pi0 vs eta)

### Dataset preparation ###
data_709_3686 = DatasetManager.load_real_data.find("709_3686")
incMC_709_3686 = DatasetManager.load_inclusive_mc.find("709_3686")

all_data = [data_709_3686]
all_incMC = [incMC_709_3686]

# Decay card: psi(3686) -> Lambda Lambdabar pi0
# Lambda -> p pi-, Lambdabar -> pbar pi+, pi0 -> gamma gamma
# (Same final state for eta mode; generated separately)
decay_card = <<~DECAYCARD
    Decay psi(3686)
    1.000  Lambda  anti-Lambda  pi0              HELAMP 1.0 0.0;
    Enddecay

    Decay Lambda
    1.000  p+  pi-                                 PHSP;
    Enddecay

    Decay anti-Lambda
    1.000  anti-p-  pi+                            PHSP;
    Enddecay

    Decay pi0
    1.000  gamma  gamma                            PHSP;
    Enddecay
End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "psi3686_to_LL_pi0_eta"
  config.events         = 500000
  config.decay_card     = decay_card
  config.cross_section  = :default
end

### Event selection ###
alg = Algorithm.new("psi3686_LL_pi0_eta")
alg.set_header(["psi3686LLPi0EtaAlg/psi3686LLPi0Eta.h"])

event_selection = Selection.new

# 4 charged tracks with net charge zero
event_selection.select_track {
                 cos_theta 0.93
                 Vz   10.0
                 Vr   1.0
               }
               .pid(method: :probability) {
                 prob_cut 0.001
                 identify :proton, against: [:pion, :kaon]
                 identify :anti_proton, against: [:pion, :kaon]
                 identify :pion, against: [:kaon, :proton]
               }
               .assign({:prp => :prp, :prm => :prm, :pip => :pip, :pim => :pim})

# Photon selection: at least 2 photons
event_selection.select_photon {
                 energy 0.025
                 emc_time [0, 700]
                 min_angle_to_charged 10   # degrees
               }

# Build Lambda from p pi- with secondary vertex fit
event_selection.build_virtual_particle(:Lambda, [:prp, :pim]) {
                 secondary_vertex_fit
                 mass_window_lo 1.111
                 mass_window_hi 1.121
               }

# Build Lambdabar from pbar pi+ with secondary vertex fit
event_selection.build_virtual_particle(:Lambda_bar, [:prm, :pip]) {
                 secondary_vertex_fit
                 mass_window_lo 1.111
                 mass_window_hi 1.121
               }

# 4C kinematic fit: energy-momentum conservation under Lambda Lambdabar gamma gamma
event_selection.kinematic_fit {
                 constrain_four_momentum
                 nominal
                 chi2_cut 200
               }

alg.with_decay_card(decay_card).apply(event_selection)

# Lambda/Lambdabar reconstruction
alg.note(:lambda_reconstruction,
  "Lambda/Lambdabar: secondary vertex fit p pi- / pbar pi+. Best Lambda-Lambdabar pair by min chi2_svtx(Lambda) + chi2_svtx(Lambdabar). 4C kinematic fit (chi2<40): energy-momentum conservation under Lambda Lambdabar gamma gamma hypothesis. J/psi veto: 3.087 < M(Lambda Lambdabar) < 3.107 GeV/c2 rejected. pi+pi- recoil mass J/psi veto: 3.087 < M_recoil(pi+pi-) < 3.107 GeV/c2 rejected. Applied in ROOT.")

# pi0 signal extraction
alg.note(:pi0_analysis,
  "psi(3686) -> Lambda Lambdabar pi0: additional chi2_4C<15, M(Lambda Lambdabar)<3.4 GeV/c2 veto, M(p pi0/pbar pi0) vetoes for background suppression. Signal yield from fit to M(gamma gamma): 23.0+/-6.3 events, significance 3.7sigma. BF = (1.42+/-0.39+/-0.59)e-6. Applied in ROOT.")

# eta signal extraction
alg.note(:eta_analysis,
  "psi(3686) -> Lambda Lambdabar eta: eta mass window (0.525,0.560) GeV/c2. N_obs = 218+/-17. Dalitz plot analysis for Lambda(1670) resonance. PWA yields: M(Lambda(1670)) = 1672+/-5 MeV/c2, Gamma = 38+/-10 MeV. BF(psi(3686)->Lambda Lambdabar eta) = (2.34+/-0.18+/-0.52)e-5. Applied in ROOT.")

# 448.1M psi(3686) events
alg.note(:total_psi3686_events,
  "Total of (448.1+/-2.9) x 10^6 psi(3686) events. Continuum data at 3.773 GeV (2.92 fb^-1) for background estimation. Applied in ROOT.")

all_datasets = all_data + all_incMC + exMC_signal
root_files = alg.execute_on(all_datasets)