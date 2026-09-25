# DSL auto-generated from 2505.05888v2
# Paper: Measurement of the relative phase between strong and electromagnetic amplitudes in J/ψ → φη
# Analysis: 26 energy points 3000-3120 MeV, ConExc signal MC (mode 23: φ η)
# Selection: K+K- + γγ, 4C kinematic fit

### Dataset preparation ###
# 26 scan points: identifiable subset from dataset tables; remaining points to be added
data_scan = [
  DatasetManager.real_data.find("713_Rscan_2950"),
  DatasetManager.real_data.find("713_Rscan_2981"),
  DatasetManager.real_data.find("713_Rscan_3000"),
  DatasetManager.real_data.find("713_Rscan_3020"),
  DatasetManager.real_data.find("713_Rscan_3080"),
]
incMC_scan = data_scan.map { |d| DatasetManager.inclusive_mc.find(d.sample_name) }

# ConExc decay card for e+e- → φη (mode 23), φ→K+K-, η→γγ
# Particle vpho omitted — DSL auto-injects per energy point for multi-energy scan
decay_card_phi_eta = <<~DECAYCARD
    Decay vpho
    1 ConExc 23;
    Enddecay
    Decay vhdr
    1 phi eta PHSP;
    Enddecay
    Decay phi
    1.000 K+ K- VSS;
    Enddecay
    Decay eta
    1.000 gamma gamma PHSP;
    Enddecay
    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc_for(data_scan) do |config|
  config.sample_name = "sig_conexc_phi_eta"
  config.events = 50_000
  config.decay_card = decay_card_phi_eta
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg = Algorithm.new("JpsiPhiEtaScan")
alg.set_header(["JpsiPhiEtaScanAlg/JpsiPhiEtaScan.h"])
    .set_constant({ "ECMS" => [:double, 3.097] })
    .note(:scan_points_partial, "26 scan points 3000-3120 MeV; only 5 identifiable from dataset tables listed. Remaining points require lookup.")
    .note(:eta_mass_window, "|M(γγ) - M_eta| < 30 MeV applied as post-fit window in ROOT analysis")
    .note(:chi2_rejection, "Paper applies chi2_4C > 85 event rejection; using loose chi2_cut 200 per Rule T3 (tight cut is ROOT-level)")
    .note(:phi_window, "φ mass window applied in ROOT analysis on fit-corrected M(K+K-)")

sel = Selection.new
sel.select_track do
      nChrp ">=1"
      nChrn ">=1"
      nNet "==0"
      cos_theta 0.93
      Vz 100.0
      Vr 10.0
    end
    .select_photon do
      tdc_emc_start 0
      tdc_emc_end 14          # EMC time [0, 700] ns (14 × 50 ns units)
      angle_to_track 10.0
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      nGam ">=2"
    end
    .pid(method: :probability) do
      prob_cut 0.001
      identify :kaon, against: [:pion, :proton]
      nkp ">=1"
      nkm ">=1"
    end
    .remove([:kp <= :chrgp])
    .remove([:km <= :chrgn])
    .assign({ chrgp: :pip, chrgn: :pim })
    .kinematic_fit([:kp, :km, :gamma, :gamma]) do
      nominal
      constrain_four_momentum
      chi2_cut 200
    end

alg.with_decay_card(decay_card_phi_eta).apply(sel)
alg.execute_on(data_scan + incMC_scan + exMC_signal)