# Paper: 2207.11666v2
# Title: First observation of psi(3686) -> Lambda Lambdabar omega and
#        search for excited Lambda states
# Energy: 3.686 GeV (psi(3686), single energy point)
# Final state: psi(3686) -> Lambda Lambdabar omega
# Lambda -> p pi-, Lambdabar -> pbar pi+, omega -> pi+ pi- pi0, pi0 -> gamma gamma
# 5C kinematic fit with pi0 mass constraint

### Dataset preparation ###
data_709_3686 = DatasetManager.load_real_data.find("709_3686")
incMC_709_3686 = DatasetManager.load_inclusive_mc.find("709_3686")

all_data = [data_709_3686]
all_incMC = [incMC_709_3686]

# Decay card: psi(3686) -> Lambda Lambdabar omega
# Lambda -> p pi-, Lambdabar -> pbar pi+, omega -> pi+ pi- pi0, pi0 -> gamma gamma
decay_card = <<~DECAYCARD
    Decay psi(3686)
    1.000  Lambda  anti-Lambda  omega             HELAMP 1.0 0.0;
    Enddecay

    Decay Lambda
    1.000  p+  pi-                                 PHSP;
    Enddecay

    Decay anti-Lambda
    1.000  anti-p-  pi+                            PHSP;
    Enddecay

    Decay omega
    1.000  pi+  pi-  pi0                           PHSP;
    Enddecay

    Decay pi0
    1.000  gamma  gamma                            PHSP;
    Enddecay
End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "psi3686_to_LL_omega"
  config.events         = 500000
  config.decay_card     = decay_card
  config.cross_section  = :default
end

### Event selection ###
alg = Algorithm.new("psi3686_LLomega")
alg.set_header(["psi3686LLomegaAlg/psi3686LLomega.h"])

event_selection = Selection.new

# At least 6 charged tracks: proton, antiproton from Lambda/Lambdabar + 2 pions from omega
event_selection.select_track {
                 cos_theta 0.93
                 Vz   10.0
                 Vr   1.0
               }
               .pid(method: :probability) {
                 prob_cut 0.001
                 identify :proton, against: [:pion, :kaon]
                 identify :pion, against: [:kaon, :proton]
                 identify :anti_proton, against: [:pion, :kaon]
               }
               .assign({:prp => :prp, :prm => :prm, :pip => :pip, :pim => :pim})

# Photon selection for pi0 -> gamma gamma
event_selection.select_photon {
                 energy 0.025
                 emc_time [0, 700]
                 min_angle_to_charged 10   # degrees
               }

# Build Lambda from p pi- with secondary vertex fit
event_selection.build_virtual_particle(:Lambda, [:prp, :pim]) {
                 secondary_vertex_fit
                 mass_window 0.005   # |M-M(Lambda)| < 5 MeV
               }

# Build Lambdabar from pbar pi+ with secondary vertex fit
event_selection.build_virtual_particle(:Lambda_bar, [:prm, :pip]) {
                 secondary_vertex_fit
                 mass_window 0.005
               }

# Build pi0 from gamma gamma with 1C kinematic fit
event_selection.build_virtual_particle(:pi0, [:gamma, :gamma]) {
                 mass_window_lo 0.115
                 mass_window_hi 0.150
               }

# 5C kinematic fit: energy-momentum conservation + pi0 mass constraint
event_selection.kinematic_fit {
                 constrain_four_momentum
                 invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                 nominal
                 chi2_cut 200
               }

alg.with_decay_card(decay_card).apply(event_selection)

# Lambda reconstruction: secondary vertex fit, mass window
alg.note(:lambda_reconstruction,
  "Lambda/Lambdabar: secondary vertex fit p pi- (chi2<50). |M-M(Lambda)| < 5 MeV. Best Lambda-Lambdabar pair by min [M(p pi-)-M(Lambda)]^2 + [M(pbar pi+)-M(Lambdabar)]^2. 5C kinematic fit (chi2<40): energy-momentum + pi0 mass constrained. Applied in ROOT.")

# omega reconstruction and J/psi veto
alg.note(:omega_reconstruction,
  "omega -> pi+ pi- pi0. pi0 -> gamma gamma with 1C kinematic fit (chi2<25). J/psi veto: |M_recoil(pi+pi-) - M(J/psi)| > 30 MeV. omega signal region: [0.753,0.813] GeV. Background from omega sidebands. Applied in ROOT.")

# omega signal extraction and BF
alg.note(:bf_measurement,
  "Unbinned ML fit to M(pi+pi-pi0): omega signal from MC shape, background = 1st order Chebyshev. N_obs = 207+/-21. Detection efficiency = 3.89%. B(psi(3686)->Lambda Lambdabar omega) = (3.30+/-0.34+/-0.29)e-5. (448.1+/-2.9)e6 psi(3686) events. Applied in ROOT.")

# Lambda* resonance search via Dalitz plot fit
alg.note(:lambdastar_search,
  "Lambda* search: 2D unbinned ML fit to M^2(Lambda omega) vs M^2(Lambdabar omega) Dalitz plot. Lambda* signal: S-wave Breit-Wigner. Fitted mass = 2.001+/-0.007 GeV, width = 0.036+/-0.014 GeV. Significance 3.0sigma. Upper limit B(psi(3686)->Lambda Lambdabar* + c.c. -> Lambda Lambdabar omega) < 1.4e-5 at 90% CL. Applied in ROOT.")

# 448.1M psi(3686) events
alg.note(:total_psi3686_events,
  "Total of (448.1+/-2.9) x 10^6 psi(3686) events used. Continuum data at 3.650 GeV (44.49 pb^-1) and 3.773 GeV (2917 pb^-1) for background estimation.")

all_datasets = all_data + all_incMC + exMC_signal
root_files = alg.execute_on(all_datasets)