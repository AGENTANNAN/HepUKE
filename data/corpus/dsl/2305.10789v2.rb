# Paper: 2305.10789v2
# Title: Precise measurement of e+e- -> D_s*+ D_s*- cross sections
# Energy: 76 energy points from 4.226 to 4.95 GeV
# Semi-inclusive method: reconstruct only one D_s*±
# D_s* -> gamma D_s, D_s -> K+ K- pi±
# ConExc generator for continuum/R-scan cross section measurement

### Dataset preparation ###
# 76 energy points from threshold to 4.95 GeV
# Representative sample at 4.226 GeV area
data_703_4226 = DatasetManager.load_real_data.find("703_4226")

all_data = [data_703_4226]
all_incMC = DatasetManager.load_inclusive_mc  # multiple energy points

# Decay card: e+e- -> D_s*+ D_s*- (ConExc: continuum production)
# D_s*+ -> gamma D_s+; D_s+ -> K+ K- pi+
# Semi-inclusive: reconstruct one D_s* side
decay_card = <<~DECAYCARD
    Decay e+ e-
    1.000  D_s*+  D_s*-                           VSP_PWAVE 1.0 0.0 0.0;
    Enddecay

    Decay D_s*+
    1.000  gamma  D_s+                            VSP_PWAVE 1.0 0.0 0.0;
    Enddecay

    Decay D_s*-
    1.000  gamma  D_s-                            VSP_PWAVE 1.0 0.0 0.0;
    Enddecay

    Decay D_s+
    1.000  K+  K-  pi+                            PHSP;
    Enddecay

    Decay D_s-
    1.000  K+  K-  pi-                            PHSP;
    Enddecay
End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "DsStarDsStar_semiinc"
  config.events         = 500000
  config.decay_card     = decay_card
  config.cross_section  = :default
end

### Event selection ###
alg = Algorithm.new("DsStarDsStar_CrossSection")
alg.set_header(["DsStarDsStarCSAlg/DsStarDsStarCS.h"])

event_selection = Selection.new

# Charged tracks: at least 2 kaons + 1 pion for D_s -> K+ K- pi±
# + 1 radiative photon
event_selection.select_track {
                 cos_theta 0.93
                 Vz   10.0
                 Vr   1.0
               }
               # PID: kaon and pion identification
               .pid(method: :probability) {
                 prob_cut 0.001
                 identify :kaon, against: [:pion, :proton]
                 identify :pion, against: [:kaon, :proton]
               }
               .assign({:kp => :kp, :km => :km, :pip => :pip, :pim => :pim})

# Photon selection for D_s* -> gamma D_s
event_selection.select_photon {
                 energy 0.025
               }

# Build D_s candidates: K+ K- pi± within |M(KKpi)-M(D_s)| < 15 MeV
# Then D_s* candidates: gamma + D_s
alg.with_decay_card(decay_card).apply(event_selection)

# Semi-inclusive reconstruction:
# - Reconstruct only one D_s*± per event
# - D_s -> K+ K- pi±, mass window |M(KKpi)-M(D_s)| < 15 MeV
# - D_s* -> gamma D_s, modified missing mass: M_miss = m_miss + m(gamma KKpi) - m(D_s*)
# - Signal region: |M_miss - m(D_s*)| < 5 * sigma_M_miss(E_CM)
# - Signal yield from fit to M(gamma KKpi) distribution
alg.note(:semi_inclusive_method,
  "Semi-inclusive: reconstruct one D_s*±. D_s mass window 15 MeV. Modified missing mass cut at 5sigma. Signal yield from M(gamma KKpi) fit: MC signal shape conv. Gaussian + random combination shape + 2nd-order Chebyshev background. Applied in ROOT.")

# Born cross section calculation:
# sigma_Born = N_sig / [2 * B(D_s->KKpi) * epsilon * (1+delta) * (1/|1+Pi|^2) * L_int]
# ISR correction via iterative weighting with fitted line shape
# VP correction applied
alg.note(:cross_section_calculation,
  "Born cross sections measured at 76 energy points. ISR correction via 4 iterations. Line shape: coherent sum of 3 BWs + PHSP amplitude. Two resonances: M1=4186.5 MeV (Gamma=55 MeV), M2=4414.5 MeV (Gamma=122.6 MeV). Third BW at ~4.79 GeV, significance >6.1sigma. Applied in ROOT.")

# 76 energy points, 15.67 fb^-1 total
alg.note(:energy_points,
  "76 c.m. energies from 4.226 to 4.95 GeV. Total integrated luminosity 15.67 fb^-1.")

all_datasets = all_data + all_incMC + exMC_signal
root_files = alg.execute_on(all_datasets)