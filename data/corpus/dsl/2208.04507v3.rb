# Paper: 2208.04507v3
# Title: Measurement of e+e- -> omega pi+ pi- cross section at
#        sqrt(s) = 2.000 to 3.080 GeV
# Energy: 2.000-3.080 GeV (19 energy points, 647 pb^-1 total)
# Final state: omega -> pi+ pi- pi0, pi0 -> gamma gamma
# Cross section measurement with PWA for intermediate subprocesses
# X(2230) resonance search

### Dataset preparation ###
# 19 energy points from 2.000 to 3.080 GeV
# Using a representative dataset; actual multi-point handling in ROOT
data_scan = DatasetManager.load_real_data.find("703_4600")

all_data = [data_scan]
all_incMC = DatasetManager.load_inclusive_mc

# Decay card: e+e- -> omega pi+ pi- with ISR
# omega -> pi+ pi- pi0, pi0 -> gamma gamma
decay_card = <<~DECAYCARD
    Decay e+ e-
    1.000  omega  pi+  pi-                        PHSP;
    Enddecay

    Decay omega
    1.000  pi+  pi-  pi0                          PHSP;
    Enddecay

    Decay pi0
    1.000  gamma  gamma                           PHSP;
    Enddecay
End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "ee_to_omega_pipi"
  config.events         = 500000
  config.decay_card     = decay_card
  config.cross_section  = :default
end

### Event selection ###
alg = Algorithm.new("ee_to_omega_pipi")
alg.set_header(["eeOmegaPiPiAlg/eeOmegaPiPi.h"])

event_selection = Selection.new

# 4 charged tracks with net charge zero
event_selection.select_track {
                 cos_theta 0.93
                 Vz   10.0
                 Vr   1.0
               }
               .pid(method: :probability) {
                 prob_cut 0.001
                 identify :pion, against: [:kaon]
               }
               .assign({:pip1 => :pip, :pim1 => :pim, :pip2 => :pip, :pim2 => :pim})

# Photon selection: at least 2 photons
event_selection.select_photon {
                 energy 0.025
                 emc_time [0, 700]
                 min_angle_to_charged 10   # degrees
               }

# Build pi0 from gamma gamma
event_selection.build_virtual_particle(:pi0, [:gamma, :gamma]) {
                 mass_window_lo 0.120
                 mass_window_hi 0.150
               }

# 4C kinematic fit imposing four-momentum conservation
# under e+e- -> 2(pi+pi-) gamma gamma hypothesis
event_selection.kinematic_fit {
                 constrain_four_momentum
                 nominal
                 chi2_cut 200
               }

alg.with_decay_card(decay_card).apply(event_selection)

# Event selection details
alg.note(:event_selection,
  "4 charged tracks (net charge zero), at least 2 photons. Track: |cos theta|<0.93, |Vz|<10 cm, |Vxy|<1 cm. Photon: E>25 MeV barrel, E>50 MeV endcap. pi0: M(gamma gamma) in (0.120,0.150) GeV/c2. 4C kinematic fit (chi2<40): e+e- -> 2(pi+pi-) gamma gamma. Best photon combination by min chi2_4C. omega: combination of pi+pi-pi0 with minimal mass difference from known omega mass, signal region (0.758, 0.808) GeV/c2. Sidebands for background estimation. Applied in ROOT.")

# Cross section measurement
alg.note(:cross_section,
  "Cross sections via fits to M(pi+pi-pi0). Signal: MC shape convolved with Gaussian. Background: 2nd-order Chebychev (dominated by e+e- -> 2(pi+pi-)pi0). ISR correction via ConExc, iterated to convergence. X(2230) observed with 10.3sigma significance. M=2250+/-25+/-27 MeV/c2, Gamma=125+/-43+/-15 MeV. Applied in ROOT.")

# PWA for intermediate subprocesses
alg.note(:pwa,
  "PWA with helicity amplitude formalism. Intermediate states: omega f0(500), omega f0(980), omega f0(1370), omega f2(1270), b1(1235)pi. f0(980) with Flatte formula. Simultaneous fit across c.m. energies. Group A: 2.000-2.232 GeV, Group B: 2.309-2.900 GeV. Non-resonant amplitude also included. Applied in ROOT.")

# 19 energy points
alg.note(:energy_points,
  "19 c.m. energies: 2.0000, 2.0500, 2.1000, 2.1250, 2.1500, 2.1750, 2.2000, 2.2324, 2.3094, 2.3864, 2.3960, 2.6444, 2.6464, 2.9000, 2.9500, 2.9810, 3.0000, 3.0200, 3.0800 GeV. Total 647 pb^-1.")

all_datasets = all_data + all_incMC + exMC_signal
root_files = alg.execute_on(all_datasets)