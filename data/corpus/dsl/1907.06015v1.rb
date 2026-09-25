# 1907.06015v1: Cross sections of e+e- → K+K-K+K- and φK+K- at √s = 2.100-3.080 GeV

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

# Representative energy point (R-scan data at 20 c.m. energies; use as template)
scan_data = DatasetManager.real_data.find("713_Rscan_2125")
scan_incMC = DatasetManager.inclusive_mc.find("713_Rscan_2125")

# ConExc decay card: mode 16 = 2(K+K-)
decay_card_4k = <<~DECAYCARD
  Decay vpho
  1 ConExc 16;
  Enddecay
  Decay vhdr
  1 K+ K- K+ K- PHSP;
  Enddecay
  End
DECAYCARD

exMC_4k = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_conexc_4k"
  config.related_dataset = scan_data
  config.events = 100_000
  config.decay_card = decay_card_4k
  config.cross_section = :default
end

# ====== Algorithm 1: e+e- → K+K-K+K- ======
alg_4k = Algorithm.new("FourKaon4K")
alg_4k.set_header(["FourKaon4KAlg/FourKaon4K.h"])
       .set_constant({"ECMS" => [:double, 2.125]})

sel_4k = Selection.new
sel_4k.select_track {
  cos_theta 0.93
  Vz 100.0
  Vr 10.0
  nChrp ">=2"
  nChrn ">=2"
  nNet "==0"
}
.pid(method: :probability) {
  prob_cut 0.001
  identify :kaon, against: [:pion, :proton]
  nkp "==2"
  nkm "==2"
}
.kinematic_fit([:kp, :kp, :km, :km]) {
  nominal
  constrain_four_momentum
  chi2_cut 200
}

alg_4k.with_decay_card(decay_card_4k).apply(sel_4k)
alg_4k.note(:signal_extraction, "Signal yield from unbinned ML fit to M_recoil(K+K-K) in ROOT; signal shape from MC convolved with Gaussian; bg: 2nd-order Chebyshev polynomial")
       .note(:qed_suppression, "Identified kaon momenta required to be < 0.8*p_beam to suppress e+e- and mu+mu- QED backgrounds")
       .note(:isr_correction, "ISR correction factor (1+deltar) obtained by QED calculation iterating until Born cross section converges")
       .note(:weighted_efficiency, "Detection efficiency combined from phiK+K- weighted PHSP MC and 4K PHSP MC weighted by signal yields")

alg_4k.execute_on([scan_data, scan_incMC, exMC_4k])

# ====== Algorithm 2: e+e- → φK+K- (φ → K+K-) ======
alg_phikk = Algorithm.new("PhiKK")
alg_phikk.set_header(["PhiKKAlg/PhiKK.h"])
         .set_constant({"ECMS" => [:double, 2.125]})

sel_phikk = Selection.new
sel_phikk.select_track {
  cos_theta 0.93
  Vz 100.0
  Vr 10.0
  nChrp ">=1"
  nChrn ">=1"
}
.pid(method: :probability) {
  prob_cut 0.001
  identify :kaon, against: [:pion, :proton]
  nkp ">=1"
  nkm ">=1"
}
# 1C kinematic fit: 3 identified kaons + 1 missing kaon (Kmiss = kaon mass)
.kinematic_fit([:kp, :km, :kp]) {
  nominal
  miss_track_of :km
  constrain_four_momentum
  chi2_cut 20
}

alg_phikk.with_decay_card(decay_card_4k).apply(sel_phikk)
alg_phikk.note(:phi_signal, "phi->K+K- signal from P-wave Breit-Wigner fit to M(K+K-) in ROOT; non-phi bg described by ARGUS function")
         .note(:qed_suppression, "Identified kaon momenta < 0.8*p_beam")
         .note(:branching_fraction_correction, "B(phi->K+K-) = 49.2% applied in cross section calculation")

alg_phikk.execute_on([scan_data, scan_incMC, exMC_4k])