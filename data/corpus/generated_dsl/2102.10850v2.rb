# ============================================================
# ψ(3770) → D0 D̄0  (double-tag / semileptonic tag, √s = 3.773 GeV)
#   tag side : D̄0 → K+π− , K+π−π0 , K+π−π−π+   (charm = −1)
#   signal   : D0 → K1(1270)− e+νe , K1(1270)− → K−π+π−
# ============================================================

### Dataset preparation ###
data_3773  = DatasetManager.real_data.find("712_3773")      # 3.773 GeV real data (ψ(3770))
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")   # corresponding inclusive MC

# Full decay card: three tag modes in equal fractions + semileptonic signal decay
decay_card = <<~DECAYCARD
  Decay psi(3770)
  1.0000 D0 anti-D0 PHSP;
  Enddecay

  Decay anti-D0
  0.3333 K+ pi- PHSP;
  0.3333 K+ pi- pi0 PHSP;
  0.3334 K+ pi- pi- pi+ PHSP;
  Enddecay

  Decay D0
  1.0000 K1(1270)- e+ nu_e PHSP;
  Enddecay

  Decay K1(1270)-
  1.0000 K- pi+ pi- PHSP;
  Enddecay

  Decay pi0
  1.0000 gamma gamma PHSP;
  Enddecay

  End
DECAYCARD

# 400k exclusive MC events for the full decay card (three equal tag modes + signal)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_D0D0bar_K1enu"
  config.related_dataset = data_3773
  config.events          = 400_000
  config.decay_card      = decay_card
  config.cross_section   = :default
end

### Tag-based event selection (TagAnalysis) ###
alg_name = "D0D0barDTK1enu"
alg = TagAnalysis.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 3.773]})
   .with_decay_card(decay_card)

# Tag side: anti-D0 (charm = −1) in the three hadronic modes
alg.tag_side(:D0) do |t|
  t.modes :D0toKPi, :D0toKPiPi0, :D0toKPiPiPi
  t.charm -1
  t.window :mBC, min: 1.858, max: 1.874   # explicitly requested M_BC window
end

# Signal side: D0 → K1(1270)− e+νe → K−π+π− e+ νe
alg.signal_side do |s|
  s.charged(km: 1, pip: 1, pim: 1, ep: 1)   # one K−, one π+, one π−, one e+
  s.require_charge 0                        # net charge zero
  s.missing :nu_e                           # one undetected νe
end

# 4C kinematic fit on the tag D plus the K−π+π−e+νe system
alg.fit do |f|
  f.constrain_four_momentum
  f.chi2_cut 200
end

# Inexpressible BOSS-side tag configuration captured as notes
alg.note(:tag_delta_e_selection, "tag-side ΔE windows are mode dependent and are applied in ROOT on the stored ΔE: K+π- in [-0.029, 0.027] GeV, K+π-π0 in [-0.069, 0.038] GeV, K+π-π-π+ in [-0.031, 0.028] GeV; within each tag mode the accepted candidate is the one with the smallest |ΔE|")
   .note(:tag_reconstruction_selection, "DTagAlg local re-reconstruction applies |cosθ|<0.93, |Vz|<10 cm and Vr<1 cm to charged tracks (tracks from K_S0 decays are not used); photons must satisfy E>25 MeV in the barrel (|cosθ|<0.80) or E>50 MeV in the endcap (0.86<|cosθ|<0.92), be separated by >10° from the nearest extrapolated track, and have an EMC shower time within 700 ns; π0 candidates are formed from two photons with 0.115<M(γγ)<0.150 GeV/c2 and a 1C mass-constrained fit to the nominal π0 mass")
   .note(:pid_definition, "K/π separation combines dE/dx and TOF, requiring L(K)>L(π) for kaons and L(π)>L(K) for pions; electron identification combines dE/dx, TOF and EMC, requiring E/p>0.8 and L'(e)/(L'(e)+L'(π)+L'(K))>0.8")

alg.apply
root_files = alg.execute_on([data_3773, incMC_3773, exMC_signal])