# Paper: 2305.12166v2
# Title: Production of doubly-charged Delta baryon in e+e- annihilation
# Energy: 6 energies from 2.3094 to 2.6464 GeV
# Final state: e+e- -> Delta++ Delta-- or Delta++ pbar pi- + c.c.
# Delta++ -> p pi+, Delta-- -> pbar pi-
# Reconstructed final state: p pbar pi+ pi-
# ConExc generator for continuum production at low energies

### Dataset preparation ###
# 6 energy points: 2.3094, 2.3864, 2.3960, 2.5000, 2.6444, 2.6464 GeV
# Use low-energy datasets
data_703_2309 = DatasetManager.load_real_data.find("703_2309")
data_703_2386 = DatasetManager.load_real_data.find("703_2386")
data_703_2396 = DatasetManager.load_real_data.find("703_2396")
data_703_2500 = DatasetManager.load_real_data.find("703_2500")
data_703_2644 = DatasetManager.load_real_data.find("703_2644")
data_703_2646 = DatasetManager.load_real_data.find("703_2646")

all_data = [data_703_2309, data_703_2386, data_703_2396,
            data_703_2500, data_703_2644, data_703_2646]
all_incMC = DatasetManager.load_inclusive_mc  # multiple low-energy points

# Decay card: e+e- -> Delta++ Delta-- (ConExc: continuum)
# Delta++ -> p pi+; Delta-- -> pbar pi-
decay_card = <<~DECAYCARD
    Decay e+ e-
    1.000  Delta++  anti-Delta--                 HELAMP 1.0 0.0;
    Enddecay

    Decay Delta++
    1.000  p+  pi+                                PHSP;
    Enddecay

    Decay anti-Delta--
    1.000  anti-p-  pi-                           PHSP;
    Enddecay
End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc_for(all_data) do |config|
  config.sample_name   = "ee_to_DeltaPP_DeltaMM"
  config.events         = 500000
  config.decay_card     = decay_card
  config.cross_section  = :default
end

### Event selection ###
alg = Algorithm.new("DeltaPPDeltaMM")
alg.set_header(["DeltaPPDeltaMMAlg/DeltaPPDeltaMM.h"])

event_selection = Selection.new

# Exactly 4 charged tracks, zero net charge
# Proton and pion PID
event_selection.select_track {
                 cos_theta 0.93
                 Vz   10.0
                 Vr   1.0
                 nTot "==4"
               }
               # PID: identify protons and pions
               .pid(method: :probability) {
                 prob_cut 0.001
                 identify :proton, against: [:kaon, :pion]
                 identify :pion, against: [:kaon]
                 nproton "==2"
                 npion   "==2"
               }
               .assign({:prp => :prp, :prm => :prm, :pip => :pip, :pim => :pim})
               # 4C kinematic fit for p pbar pi+ pi- hypothesis
               .kinematic_fit([:prp, :prm, :pip, :pim]) {
                 nominal
                 constrain_four_momentum
                 chi2_cut 50
               }

alg.with_decay_card(decay_card).apply(event_selection)

# Delta++ -> p pi+ and Delta-- -> pbar pi- reconstruction
# Signal extraction via 2D simultaneous unbinned ML fit:
# m(p pi+) vs m(pbar pi-) and m(p pi-) vs m(pbar pi+)
# Five components in fit: Delta++ Delta-- (signal), PHSP (p pbar pi+ pi-),
# semi-Delta (Delta++ pbar pi- and p pi+ Delta--), Lambda anti-Lambda
alg.note(:signal_extraction,
  "Signal extraction: simultaneous 2D unbinned ML fit to m(p pi+) vs m(pbar pi-) and m(p pi-) vs m(pbar pi+). 5 components: signal Delta++Delta--, PHSP, semi-Delta (2 modes), Lambda anti-Lambda. Broad resonance fits with MC shapes. No significant Delta++Delta-- signal; significant semi-Delta above 2.6 GeV. Applied in ROOT.")

# Born cross sections:
# sigma_Born = N_sig / [L_int * epsilon * (1+delta) * Br]
# ConExc generator with ISR and VP corrections
# Upper limits via Bayesian method for non-significant signals
alg.note(:cross_section_results,
  "Born cross sections / upper limits at 6 energy points. Delta++Delta--: no significant signal, upper limits at 90% CL. Semi-Delta (Delta++ pbar pi- + c.c.): significant at 2.6444 and 2.6464 GeV. Combined 2.6454 GeV: sigma_Born = 58.2+/-5.0+/-5.9 pb. Applied in ROOT.")

# 6 energy points, 179 pb^-1 total
alg.note(:energy_points,
  "6 c.m. energies: 2.3094, 2.3864, 2.3960, 2.5000, 2.6444, 2.6464 GeV. Total 179 pb^-1.")

all_datasets = all_data + all_incMC + exMC_signal
root_files = alg.execute_on(all_datasets)