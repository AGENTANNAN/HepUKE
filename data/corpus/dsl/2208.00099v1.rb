# Paper: 2208.00099v1
# Title: Measurement of e+e- -> pi+ pi- D+ D- cross sections at
#        center-of-mass energies from 4.190 to 4.946 GeV
# Energy: 4.190-4.946 GeV (37 energy points)
# Partial reconstruction method:
#   D+ -> K- pi+ pi+ (tag side)
#   pi+ pi- from remaining tracks
#   D- identified via recoil mass of D+ pi+_d pi-_d
# Cross section measurement; resonance search (Y(4390), new state at 4.706 GeV)
# 17.4 fb^-1 total

### Dataset preparation ###
# 37 energy points from 4.190 to 4.946 GeV
data_703_4600 = DatasetManager.load_real_data.find("703_4600")

all_data = [data_703_4600]
all_incMC = DatasetManager.load_inclusive_mc

# Decay card: e+e- -> pi+ pi- D+ D- (phase space)
# D+ -> K- pi+ pi+, D- treated via recoil (partial reconstruction)
decay_card = <<~DECAYCARD
    Decay e+ e-
    1.000  pi+  pi-  D+  D-                       PHSP;
    Enddecay

    Decay D+
    1.000  K-  pi+  pi+                           PHSP;
    Enddecay

    Decay D-
    1.000  K+  pi-  pi-                           PHSP;
    Enddecay
End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "ee_to_pipi_DpDm"
  config.events         = 500000
  config.decay_card     = decay_card
  config.cross_section  = :default
end

### Event selection ###
alg = Algorithm.new("ee_to_pipi_DpDm")
alg.set_header(["eePiPiDpDmAlg/eePiPiDpDm.h"])

event_selection = Selection.new

# Charged track selection
event_selection.select_track {
                 cos_theta 0.93
                 Vz   10.0
                 Vr   1.0
               }
               .pid(method: :probability) {
                 prob_cut 0.001
                 identify :kaon, against: [:pion]
                 identify :pion, against: [:kaon]
               }
               .assign({:km => :km, :pip1 => :pip, :pip2 => :pip, :pip_d1 => :pip, :pim_d1 => :pim})

# Build D+ candidate from K- pi+ pi+ with vertex fit
event_selection.build_virtual_particle(:D_plus, [:km, :pip1, :pip2]) {
                 secondary_vertex_fit
                 mass_window 0.011   # |M-M(D+)| < 11 MeV/c2
               }

# Partial reconstruction: D- inferred from recoil mass of D+ pi+_d pi-_d
event_selection.partial_rec {
                 set :D_minus, method: :recoil_mass, of: [:D_plus, :pip_d1, :pim_d1]
                 mass_window 0.009   # recoil mass window ~9 MeV/c2
               }

alg.with_decay_card(decay_card).apply(event_selection)

# Track selection details
alg.note(:track_selection,
  "Charged tracks: |cos theta|<0.93, |Vz|<10 cm (tight for prompt tracks), |Vxy|<1 cm. Tighter Vxy<0.55 cm and Vz<3 cm after optimization. Proton/anti-proton veto to suppress Lambda_c backgrounds. K/pi PID via likelihood comparison. Applied in ROOT.")

# D+ reconstruction and D- identification
alg.note(:d_reconstruction,
  "D+: K- pi+ pi+ combination with vertex fit chi2<100. D-: recoil mass of D+ pi+_d pi-_d system, mass window d_RM=6-9 MeV/c2 depending on sqrt(s). D+ mass constrained kinematic fit applied. D0 veto: |M(K- pi+_d pi-_d pi+)-M(D0)|>0.01 GeV/c2. K_S0 veto: secondary vertex fit on pi+_d pi-_d pair, |L/Delta_L|<2. Applied in ROOT.")

# Cross section measurement
alg.note(:cross_section,
  "Cross section measurement via RM(D+ pi+_d pi-_d) fits in D+ signal and sideband regions. Signal: MC shape convolved with Gaussian. Background: 2nd-order Chebychev. ISR correction via kkmc. Two resonance structures observed: Y(4390) with M=4373.1+/-4.0 MeV/c2, Gamma=146.5+/-7.4 MeV; new state at M=4706+/-11 MeV/c2, Gamma=45+/-28 MeV (4.1sigma). Applied in ROOT.")

# X(3842) search
alg.note(:x3842_search,
  "X(3842) search via RM(pi+_d pi-_d) distributions (equivalent to M(D+ D-)). Evidence with 4.2sigma found for sqrt(s)=4.600-4.700 GeV. N_sig=155+/-38. Applied in ROOT.")

# 37 energy points
alg.note(:energy_points,
  "37 c.m. energies: 4.190, 4.200, 4.210, 4.220, 4.230, 4.237, 4.245, 4.246, 4.260, 4.270, 4.280, 4.290, 4.310, 4.315, 4.340, 4.360, 4.380, 4.390, 4.400, 4.420, 4.440, 4.470, 4.530, 4.575, 4.600, 4.612, 4.620, 4.640, 4.660, 4.680, 4.700, 4.740, 4.750, 4.780, 4.840, 4.914, 4.946 GeV. Total 17.4 fb^-1.")

all_datasets = all_data + all_incMC + exMC_signal
root_files = alg.execute_on(all_datasets)