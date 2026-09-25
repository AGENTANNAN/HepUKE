# 1902.09756v1: psi(3686) -> p pbar phi, phi -> K+ K-
# Data: psi(2S), 4.48 x 10^8 events, 3.686 GeV
# Final state: p pbar K+ K- (4 charged tracks)
# Missing kaon: 1C kinematic fit with missing K+ or K-

DatasetManager.load_real_data('config/BES3_dataset.md')
DatasetManager.load_inclusive_mc('config/BES3_incMC.md')

data_psip = DatasetManager.real_data.find("709_3686")
incMC_psip = DatasetManager.inclusive_mc.find("709_3686")

decay_card = <<~DECAYCARD
  Decay psi(2S)
  1.000 anti-p- p+ phi PHSP;
  Enddecay

  Decay phi
  1.000 K+ K- VSS;
  Enddecay

  End
DECAYCARD

exMC = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_psi2StoPPbarPhi"
  config.related_dataset = data_psip
  config.events = 100000
  config.decay_card = decay_card
  config.cross_section = :default
end

alg = Algorithm.new("Psi2SPPbarPhi")
alg.set_header(["Psi2SPPbarPhiAlg/Psi2SPPbarPhi.h"])
    .set_constant({"ECMS" => [:double, 3.686]})

sel = Selection.new
sel.select_track {
      cos_theta 0.93
      Vz 100.0
      Vr 10.0
      nChrp ">=2"
      nChrn ">=2"
      nNet "==0"
    }
    .pid(method: :probability) {
      prob_cut 0.001
      identify :proton, against: [:kaon, :pion]
      identify :kaon, against: [:pion, :proton]
      nprp ">=1"
      nprm ">=1"
      # At least one kaon identified (events with 3 or 4 charged tracks)
      nkp ">=0"
      nkm ">=0"
    }
    .remove([:prp <= :chrgp, :prm <= :chrgn])
    .remove([:kp <= :chrgp, :km <= :chrgn])
    .assign({:chrgp => :pip, :chrgn => :pim})
    .kinematic_fit([:prp, :prm, :kp, :km]) {
      nominal
      miss_track_of :km
      constrain_four_momentum
    }
alg.with_decay_card(decay_card).apply(sel)

alg.note(:Lambda_veto,
  "Lambda/Lambda-bar veto: events excluded if any p pi- or pbar pi+ combination " \
  "has invariant mass within |M_{p pi- (pbar pi+)} - M_{Lambda(Lambda-bar)}| < 3 MeV/c^2; " \
  "remaining tracks assumed as pions for veto purpose")
alg.note(:one_c_kinematic_fit,
  "1C kinematic fit under psi(2S)->p pbar K+ K- with missing K+ or K-; " \
  "for events with both kaons detected, both fits performed and best chi2_1c kept")
alg.note(:continuum_background,
  "Continuum background e+e- -> p pbar phi evaluated using off-resonance " \
  "sample at 3.773 GeV, scaled by luminosity and cross-section ratios")
alg.note(:phi_signal_extraction,
  "Phi signal extracted by unbinned maximum likelihood fit to K+K- invariant " \
  "mass in [0.985, 1.115] GeV/c^2; signal shape from MC convoluted with Gaussian; " \
  "continuum background shape and yield fixed; non-peaking background from ARGUS function")
alg.note(:Xppbar_upper_limit,
  "Upper limit on X(ppbar): X(ppbar) parameterized by S-wave BW with FSI factor; " \
  "B(psi(2S)->X(ppbar)phi->p pbar phi) < 1.82 x 10^-7 at 90% CL; " \
  "no significant near-threshold enhancement observed in p pbar mass spectrum")

alg.execute_on([data_psip, incMC_psip, exMC])